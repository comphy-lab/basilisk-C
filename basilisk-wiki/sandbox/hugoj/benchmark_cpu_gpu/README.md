# Benchmark of Basilisk

This folder aims at benchmarking Basilisk on several hardware : my laptop, the
sandbox and Jean Zay HPC. Both cpu and gpu. The goal is to use the multilayer
solver, but other solver are tested.

## Method

The results should show : 
- the speedup of GPU vs CPU (if any)
- the evolution of the speed for the simulation
- the profiling trace 

Step 1 : gpu tests that comes with Basilisk. In a first time 'turbulence.c'
only. I use a modified version of the original example that only output a
snapshot at the beginning and the end of the simulation (instead of a movie).
Step 2 : "real" case using breaking.c example then bigger case with increased horizontal resolution (not vertical as the influence of the number of layer should be tested separately). The code is not exactly the one from the example: it has been striped down of the diagnostics to avoid any slow down by transfer from gpu to cpu.

The procedure of a benchmark is to compile the code first (with tracing enable,
without generating figure or visual output and in single precision), then run
the code a few times at different resolution.

Get the cpu name with
```bash
cat /proc/cpuinfo | grep -i 'model name'
```


Get the gpu name with
```bash
glxinfo -B
```
or if it's a nvidia gpu
```bash
nvidia-smi
```

if you have a laptop, select the dedicated gpu with
```bash
__NV_PRIME_RENDER_OFFLOAD=1 __GLX_VENDOR_LIBRARY_NAME=nvidia glxinfo -B
```
`


## Results

Here plot the figure

___

### Local


First, check that if your laptop has a way to reduce energy consumption with a "battery saver"
option, make sure you disable this to use all the power available. Also plug it
to the wall to not run on battery.

You also want to make sure that you have a dedicated graphic card. Basilisk is
known to run on both AMD and NVIDIA (gaming) gpus but nvidia cards are more combat
proven. Make sure you have the latest drivers too.

My hardware is [this](results_bench_hardware), see also [this](https://www.notebookcheck.net/Nvidia-RTX-PRO-1000-Blackwell-Generation-Laptop-Benchmarks-and-Specs.997331.0.html)


#### Compile

See this [link](https://basilisk.fr/src/INSTALL)

```bash
darcs clone https://basilisk.fr/basilisk
cd basilisk/src
```
Then
```bash
export BASILISK=path/to/basilisk/src
```

Compile the GPU (OpenGL) lib on Debian system
```bash
apt install libglfw3-dev
cd $BASILISK/grid/gpu
make
```

Compile the GPU (CUDA) lib on Debian system
```bash
apt install nvidia-cuda-toolkit # can be reduced in size ...
cd $BASILISK/grid/cuda
make
```
The compilation (you might need to modify the Makefile) command is
```bash
cc -I/home/jacqhugo/basilisk/src -I.. -DSINGLE_PRECISION -O2 -g -Wall -pipe -D_FORTIFY_SOURCE=2   -c -o cuda.o cuda.c
ar cr libbuda.a cuda.o
```

Compile the GPU (HIP) lib on Debian system. If you have a nvidia GPU, then you
will need the nvidia-cuda-toolkit package.
```bash
apt install hipcc
apt install cuda-driver-toolkit # if not done already
cd $BASILISK/grid/cuda
make
```
The compilation (you might need to modify the Makefile) command is
```bash
cc -I/home/jacqhugo/basilisk/src -I.. -DSINGLE_PRECISION -O2 -D__HIP_PLATFORM_NVIDIA__ -g -D_FORTIFY_SOURCE=2 -c -o hip.o hip.c
```
Note: you might need to `-I/path/to/hipcc/headers`. I used `hipcc` instead of
`cc` but it should not be necessary.


***Notes***

- if your GPU is recent, you may get this error:
```
(fragment shader):399: CUDA: nvrtc: error: invalid value for --gpu-architecture (-arch)
```
Then you need to modify cuda.c and hip.c to force an older
architecture

```c
sprintf(arch, "--gpu-architecture=compute_89");

```
just after the architecture detection, so after the line
```c
sprintf (arch, "--gpu-architecture=compute_%d%d", major, minor);
```

- make sure that the nvidia driver is updated
```bash
apt install nvidia-drivers
```

#### Run the tests

Run the simple script `benchmark_local.sh` (I'm using a NVIDIA card, you might
want to modify it if you use AMD). 
### Sandbox

The cpu is (16 cores)
```bash
model name	:
```

The gpu is [this](https://www.techpowerup.com/gpu-specs/geforce-rtx-4090.c3889)
```bash
OpenGL renderer string: NVIDIA GeForce RTX 4090 D/PCIe/SSE2 (stokes.lmm.jussieu.fr)
Dedicated video memory: 24564 MB
```

**Step 1**

See the results [here](https://basilisk.fr/src/grid/gpu/Benchmarks/turbulence)
for both cpu and gpu.

**Step 2**
- CPU
- GPU

___


### Jean Zay

The plan is to compile Basilisk on Jean-Zay (qcc), build the library GLFW,
compile the code and then run it.

The GPU is [this](https://www.techpowerup.com/gpu-specs/a100-pcie-40-gb.c3623)

- ***Step 1*** : compile Basilisk on the HPC

install darcs locally
```bash
sudo apt install darcs
```

clone locally the repo using darcs (the tar version is not updated)
```bash
darcs clone https://basilisk.fr/basilisk
```
updating the repo is as easy as
```bash
darcs pull
```

copy on JeanZay the `/src` folder with something like
```bash
cd
cd basilisk
tar -czvf src.tar.gz src/
scp -r src id@jean-zay.idris.fr:jean_zay_home/basilisk/src
```

then connect to Jean Zay and untar the `/src` folder
```bash
ssh -X id@jean-zay.idris.fr
cd basilisk/
tar -xvf src.tar.gz
```

you can now compile the code using the file
[compile_bas_JZ_cuda_opengl.slurm](compile_bas_JZ_cuda_opengl.slurm). You can do
it using a job
```bash
sbatch compile_bas_JZ_cuda_opengl.slurm
```
or just as a script (not ideal but working)
```bash
chmod +x compile_bas_JZ_cuda_opengl.slurm
./compile_bas_JZ_cuda_opengl.slurm > logs.compile_bas
```

export Basilisk's path by adding something like this to your $HOME/.bashrc
```bash
echo 'export BASILISK=$HOME/basilisk/src' >> ~/.bashrc
```

___

- ***Step 2*** : Compile GLFW (WIP)

Clone the repo (e.g in $HOME/lib)
```bash
cd
mkdir lib
git clone https://github.com/glfw/glfw.git
```

We need the shared libraries (libglfw.so) so you to set `-DBUILD_SHARED_LIBS=ON` in the submitted job. See the example file [here]{./compileJZ_glfw.slurm}. To compile GLFW, submit the job with
```bash
sbatch compileJZ_glfw.slurm
```

Copy this slurm file [compileJZ_glfw.slurm](compileJZ_glfw.slurm) in `~/lib/` and submit the job
```bash
sbatch compileJZ_glfw.slurm
```

GLFW libs are now available in `~/lib/glfw/install`

___

- ***Step 3*** : compile the code

We want to compile the code for each backend. Examples will use the CUDA
backend.

Copy to the HPC the simple test case
[`turbulence_benchmark.c`](turbulence_benchmark.c) (a modified version of this
[turbulence.c](https://basilisk.fr/src/examples/turbulence.c) that just outputs
snapshot instead of a movie)

Copy the slurm file [compile_JZ_cuda.slurm](compile_JZ_cuda.slurm) in the same
folder as the `.c` source file and
submit the job
```bash
sbatch compile_JZ_cuda.slurm
```

___

- ***Step 4***: run the code

Copy this slurm file [run_JZ_cuda.slurm](run_JZ_cuda.slurm)
Submit the job
```bash
sbatch run_JZ_cuda.slurm
```


***JEAN ZAY TIPS***
- see loading of machine : `module load python` then `slurmtop`

### Adastra

noeud de visu gpu nvidia ?

### Gricad 

I use BIGFOOT.
The plan is to use guix to get the desired compilation environnement, build
basilisk, build GLFW, compile the turbulence.c code, run the code.

The GPU is [this](https://www.techpowerup.com/gpu-specs/tesla-v100-pcie-16-gb.c2957)(16Go version) or [this](https://www.techpowerup.com/gpu-specs/tesla-v100-pcie-32-gb.c3184)(32Go version)

As of 12/06/26: the code runs with CUDA, not OpenGL

- ***Step 1***: install the right environnement

user doc : [link](https://gricad-doc.univ-grenoble-alpes.fr/en/hpc/softenv/guix/)

first time using guix, run (quite long)
```bash
source /applis/site/guix-start.sh
```
update packages (quite long)
```bash
guix pull
guix package -u
```
We now need to install packages in dedicated environnement. For gcc, you can run
```bash
guix install -p $GUIX_USER_PROFILE_DIR/basilisk_profile gcc-toolchain@14.3.0
```
This will install gcc-toolchain version 14.3.0 in the profile basilisk_profile.

You can also install imagemagick
```bash
guix install -p $GUIX_USER_PROFILE_DIR/basilisk_profile imagemagick
```

Using the nvidia-cuda-toolkit can be done with 
```bash
source /applis/environments/cuda_env.sh 11.7
```
and the available version can be listed with 
```bash
source /applis/environments/cuda_env.sh -l
```
<!--
You can install darcs, the version control system that Basilisk uses
```bash
guix install -p $GUIX_USER_PROFILE_DIR/basilisk_profile darcs
```
-->

Running Basilisk simulation requires:
on cpu
```bash
guix install -p $GUIX_USER_PROFILE_DIR/basilisk_profile gcc-toolchain@14.3.0
guix install -p $GUIX_USER_PROFILE_DIR/basilisk_profile openmpi@4.1.6
```
on gpu (OpenGL), add
```bash
guix install -p $GUIX_USER_PROFILE_DIR/basilisk_profile glfw@3.4
```
on CUDA, you need to source the cuda profile (see above), and run on a gpu node.
The code needs libcuda.so that is loaded **only** when a GPU ressource is
dedicated (this lib comes with the driver, not with the cuda tool kit). This
means that you won't be able to compile on the frontale node, you need to ask
for a GPU.

Then, to load the profile, do
```bash
refresh_guix basilisk_profile
```

Now `which gcc` will point to the version 14.3.0 that we installed previously.
Note that this loading procedure needs to be added to oar files ! see
FICHIER_COMPILE_BIGFOOT

___

- ***Step 2***: compile Basilisk (qcc)
Install steps : [link](https://basilisk.fr/src/INSTALL)
Question : could be done using darcs ? -> install darcs package guix ?

on your local machine (do not use the tar version, its outdated)
```bash
darcs clone https://basilisk.fr/basilisk
```
then copy it to bigfoot through ciment bastion (I copied on $HOME)
```bash
scp -r basilisk user@bigfoot.ciment:path/on/remote
```
copy also the config file (WIP, for now its the same as config.gcc)
```bash
scp config.bigfoot user@bigfoot:path/on/remote/src/
```
connect to Bigfoot

Copy the file `config.bigfoot` inside `basilisk/src`
`config file here`

then use the following [compileBAS_bigfoot.sh](compileBAS_bigfoot.sh
) script
```bash
cd basilisk/src/
```

Note that we load the profile with gcc

submit compile job
```bash
ln -sf config.bigfoot config
chmod +x compile_bas_BIGFOOT_cuda_opengl.sh
oarsub -S ./compile_bas_BIGFOOT_cuda_opengl.sh
```

At this point you can add the path to your bashrc
```bash
echo "export BASILISK=/home/user/basilisk/src >> ~/.bashrc
```
___

- ***Step 3***: compile GLFW lib MAYBE NOT NECESSARY -> GUIX PACKAGE

<!--

Doc : [link](https://www.glfw.org/docs/latest/compile_guide.html)
There might be a better way to do this using guix but I havent found it...



On your $HOME, you can create a `lib` folder:
```bash
cd 
mkdir lib
git clone https://github.com/glfw/glfw.git
```
Copy this `compile_glfw.sh` script in the `lib` folder
`script here`
Note that we load the profile with the same gcc version as the compilation for
qcc.

Submit the job
```bash
oarsub -S ./compile_glfw.sh
```
-->
___

- ***Step 4***

Copy the file [`turbulence_benchmark.c`](turbulence_benchmark.c) (a modified version of this [turbulence.c](https://basilisk.fr/src/examples/turbulence.c)) to the workspace (Bettik/Silenus), for me its
in `/bettik/PROJECTS/pr-data-ocean/jacqhugo/`.
```bash
scp turbulence_benchmark.c jacqhugo@bigfoot.ciment:/home/jacqhugo/work/bench_cleanr/turbulence.c
```


You need to compile the program. I use the development queue for this, see the
`-t devel` option in the script. If you want to compile the openGL version (Work
in progress), use the `_gpu` version of the compile and run script. If you want
to use the CUDA version, use `_cuda` and if you want to use the HIP version, use
`_hip`. Examples are shown using the CUDA version.
```bash
chmod +x cuda_BIGFOOT_compile.sh
./cuda_BIGFOOT_compile.sh
```

Then use the OAR script [cuda_BIGFOOT_run.sh](cuda_BIGFOOT_run.sh) 
```bash
chmod +x cuda_BIGFOOT_run.sh
oarsub -S ./cuda_BIGFOOT_run.sh`
```

This script will run the program at different resolutions

working logs:
- OpenGL version :
8/6/26: quand je run le program j'ai 
error: XDG_RUNTIME_DIR is invalid or not set in the environment.
11/6/26: je déclare la variable, donc plus le problème. Mais le code crash sans
erreur.

***GRICAD TIPS***
- user `oarstat -u user` to see pending jobs
- use `chandler` to see the current usage of the machine
- user `recap.py` (as a command) to see a recap of the machine hardware
- use `guix search packagename` to search for a package
- when doing development stuff, use the inscruction `#OAR -t devel` in your
script to launch jobs faster

### KRAKEN 

maybe easier than bigfoot, can run on cpu partition only

## Todo

- amd rocm sur gricad ??
- version of nvidia driver on the sandbox ?
- use darcs everywhere it runs. Less trouble for keeping up with version ...
- use containers for: portability, ease of use
- use profiling tools from nvidia, amd ?
- OpenGL code compiles but when run it crashes immediately without explicit line
  (both JeanZay and Gricad/Bigfoot)
