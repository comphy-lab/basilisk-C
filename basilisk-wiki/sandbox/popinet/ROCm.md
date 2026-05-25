# Porting Basilisk to ROCm according to ChatGPT

The key thing to understand is that `grid/gpu/grid.h` is *not* “a CUDA/OpenCL backend”. It is really:

1. a **runtime kernel-generation system**,
2. tightly coupled to Basilisk’s `foreach` IR/preprocessor,
3. using GLSL compute shaders as a portable execution backend,
4. with automatic CPU fallback and synchronization.

A ROCm/HIP port can work, but the design choices depend on how compatible you want to remain with the current Basilisk GPU model.

The original file is here: [Basilisk GPU grid.h](https://basilisk.fr/src/grid/gpu/grid.h)

A few important observations from the current implementation:

* kernels are generated as strings (`build_shader()`)
* compiled dynamically at runtime
* fields live in SSBOs
* reductions are custom GPU kernels
* synchronization is automatic through `s.gpu.stored`
* unsupported loops transparently fall back to CPU
* execution model assumes GLSL/OpenGL compute dispatch (`glDispatchCompute`) ([Basilisk][1])

A ROCm version should preserve the *semantic contract* more than the implementation details.

# The biggest architectural decision

You have two fundamentally different approaches.

## Option A — HIP runtime codegen (closest to existing design)

Replace:

* GLSL codegen
* OpenGL shader compilation
* SSBO management

with:

* HIP C++ kernel generation
* `hiprtc` runtime compilation
* HIP device buffers

This preserves Basilisk’s dynamic execution model almost exactly.

You would still:

* generate kernels from `foreach`
* compile at runtime
* cache compiled kernels
* dispatch dynamically

This is probably the cleanest conceptual mapping.

## Option B — static AOT compilation (better performance, harder integration)

Generate HIP kernels ahead-of-time during `qcc`.

This is more “native ROCm style”, but much more invasive because Basilisk currently assumes:

* loops can be generated dynamically
* loop variants depend on runtime metadata
* externals/reductions are assembled dynamically

You would basically be redesigning the GPU backend.

Unless you want a long-term rewrite, I strongly recommend Option A first.

# Mapping OpenGL concepts to ROCm/HIP

Here’s the practical translation table.

| Current GLSL/OpenGL | ROCm/HIP equivalent              |
| ------------------- | -------------------------------- |
| SSBO                | `hipMalloc` buffers              |
| shader              | HIP kernel                       |
| GLSL source string  | HIP C++ source string            |
| GLSL compiler       | `hiprtc`                         |
| `glDispatchCompute` | `hipModuleLaunchKernel`          |
| uniforms            | kernel parameters                |
| `glMemoryBarrier`   | `hipDeviceSynchronize` / streams |
| reductions          | HIP reduction kernels            |
| shader cache        | HIP module cache                 |

# What I would preserve unchanged

These parts are actually very good abstractions:

## 1. `External`

Keep this almost exactly.

It is already basically a kernel ABI descriptor.

This part:

```c
typedef struct {
  ...
} External;
```

is backend-independent.

Good design.

## 2. `ForeachData`

Also backend-independent.

The GPU backend should remain a consumer of loop metadata.

## 3. CPU/GPU hybrid mode

This is one of the best features of the current implementation.

Keep:

```c
s.gpu.stored
```

and the synchronization logic.

You will need:

```c
hipMemcpyAsync(...)
```

instead of buffer mapping.

# The hardest part

The hardest part is NOT kernel launch.

It is:

# translating generated GLSL semantics into HIP C++

The generated GLSL currently assumes:

* GLSL vector types
* GLSL intrinsics
* GLSL indexing semantics
* shader-global variables
* compute-grid IDs
* OpenGL memory layout

You will need a compatibility layer.

# I would introduce a GPU backend abstraction

Right now `grid.h` is very OpenGL-specific.

You want:

```c
typedef struct GPUBackend {
    void * (*compile)(const char * source);
    void (*launch)(...);
    void * (*alloc)(size_t);
    void (*copy_to_gpu)(...);
    void (*copy_from_gpu)(...);
} GPUBackend;
```

Then:

* GLSL backend
* HIP backend
* maybe CUDA backend later

become implementations.

This prevents another monolithic file.

# Runtime compilation with HIPRTC

This is probably the core mechanism you want.

ROCm provides:

* `hiprtcCreateProgram`
* `hiprtcCompileProgram`
* `hiprtcGetCode`

Very analogous to shader compilation.

Your current:

```c
Shader * compile_shader(...)
```

becomes:

```c
HipKernel * compile_kernel(...)
```

with:

* generated HIP source
* runtime compile
* cached module/function handles

# Replace GLSL thread indexing

Current GLSL:

```glsl
gl_GlobalInvocationID
```

HIP:

```cpp
int i = blockIdx.x * blockDim.x + threadIdx.x;
int j = blockIdx.y * blockDim.y + threadIdx.y;
```

I would define compatibility macros:

```cpp
#define GPU_I (...)
#define GPU_J (...)
```

so the code generator stays cleaner.

# Memory layout advice

This matters enormously for AMD GPUs.

Basilisk fields are already fairly GPU-friendly because they are SoA-ish.

But AMD likes:

* coalesced accesses
* aligned vector loads
* predictable strides

You should:

## preserve contiguous scalar storage

Current SSBO layout is already reasonable.

Do NOT switch to pointer-heavy structures.

# Watch out for dynamic branching

Basilisk stencils often contain:

```c
if (condition)
```

NVIDIA tolerates divergence fairly well.

AMD wavefronts are usually more sensitive.

You may want:

* stencil specialization
* branch flattening
* compile-time constants

later.

But do not optimize too early.

# Reductions

Current GLSL reductions are fairly custom.

For HIP:

* use block reductions
* shared memory
* wave intrinsics

AMD wavefront size is 64, not 32.

That matters a lot.

Use:

```cpp
__shfl_down
```

carefully through HIP abstractions.

# Important ROCm-specific issue

The current backend relies heavily on:

```c
runtime string compilation
```

ROCm runtime compilation exists, but:

* compile latency can be large
* caching matters a LOT

You absolutely want:

## persistent kernel cache

Hash:

* generated source
* compilation flags
* GPU architecture

Then cache binaries on disk.

Otherwise startup overhead will be painful.

# Suggested incremental roadmap

This is the safest route.

## Phase 1 — minimal HIP backend

Goal:

* execute simple `foreach`
* scalar fields only
* no reductions
* no adaptivity

Just prove:

```c
foreach()
  a[] = b[] + c[];
```

works.

## Phase 2 — synchronization model

Implement:

```c
gpu_cpu_sync()
```

using HIP memory copies.

## Phase 3 — reductions

Then:

* `reduction(+:x)`
* norms
* statistics

## Phase 4 — multigrid support

Only after basic kernels are stable.

## Phase 5 — optimize

Only then:

* streams
* async overlap
* wavefront tuning
* LDS/shared-memory tiling

# One thing I would NOT do

I would NOT attempt:

* HIPIFY of generated GLSL
* SPIR-V translation
* OpenCL interop

Those routes become extremely messy because Basilisk generates kernels dynamically with custom semantics.

Generate HIP source directly instead.

# A very practical suggestion

Start by extracting the backend-independent pieces from `grid/gpu/grid.h` into:

```text
gpu/common/
```

Then create:

```text
gpu/glsl/
gpu/hip/
```

Otherwise the file will become impossible to maintain.

# Overall assessment

This is a feasible project because Basilisk already has:

* explicit loop IR
* centralized kernel generation
* explicit field metadata
* clean synchronization semantics

Those are the hard parts.

The main work is:

1. building a HIP code generator,
2. reproducing runtime compilation,
3. reproducing the memory model,
4. preserving fallback behavior.

The current GLSL backend is actually an excellent prototype architecture for a HIP backend.

[1]: https://basilisk.fr/src/grid/gpu/grid.h?utm_source=chatgpt.com "Basilisk - src/grid/gpu/grid.h"
