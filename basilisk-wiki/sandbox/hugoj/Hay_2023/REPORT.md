# Reproducing Hay et al. (2023) — Europa's ocean rotating convection in Basilisk

**Status:** development plan, awaiting approval before coding.
**Target paper:** Hay, Fenty, Pappalardo & Nakayama (2023), *JGR Planets* 128, e2022JE007648 — "Turbulent Drag at the Ice-Ocean Interface of Europa in Simulations of Rotating Convection" (`~/papiers/Hay_2023_Europa.pdf`).
**Codebase:** Basilisk (`$BASILISK`), multilayer solver; structure modeled on `$BASILISK/examples/held-suarez.c`, drag pattern from `examples/global-tides.c` / `examples/ocean.h`.

---

## 1. Goal

Reproduce, with **general (regime-level) similarity**, the flow of Hay et al. (2023): a global spherical-shell ocean on Europa, heated from below and cooled at the top by a fixed temperature contrast ΔT, rotating at Europa's spin rate, with **quadratic turbulent drag (cD = 0.002) at both the seafloor and the ice–ocean interface**.

Success criteria for the first milestone (go/no-go gate):

1. Rotating convection spins up **time-mean zonal jets** from uniform ΔT forcing;
2. Low-Ra runs show **equatorial superrotation**, preferentially outside the tangent cylinder latitude θ_tc ≈ 20.6°;
3. The diagnosed **ice–ocean torque** (eq. 5 of the paper) is of plausible magnitude, with the correct regime behavior vs Ra (sign reversal near Ra ≈ 1.5 × 10⁷ as a stretch goal).

The **non-hydrostatic extension** (`layered/nh.h`) is deliberately excluded from the first milestone and will be added later if the hydrostatic regime reproduction proves the framework sound.

## 2. What the paper does

MITgcm, global cubed-sphere, **non-hydrostatic**, Boussinesq, full 3D Coriolis, linear EOS in temperature only. 100 km ocean on a 1561 km body, 58 vertical levels, ~3.4 km horizontal. Free surface, rigid seafloor. Key numbers:

| Parameter | Symbol | Value | Basilisk mapping |
|---|---|---|---|
| Radius | R | 1561 km | `Radius = 1561000` (length unit = m) |
| Rotation period | P | 306 000 s | `Omega = 2*pi/306000` ≈ 2.053e-5 s⁻¹ |
| Gravity | g | 1.3 m/s² | `G = 1.3` |
| Ocean thickness | D | 100 km | flat `zb[] = -100000` |
| Density | ρ₀ | 1000 kg/m³ | reference (`rho0`) |
| Thermal expansion | α | 2 × 10⁻⁴ K⁻¹ | EOS `drho(T)` |
| Viscosity = diffusivity | ν = κ | 61.6 m²/s | vertical + horizontal diffusion |
| Drag coefficient | cD | 0.002 (top & bottom) | `K0()` macro |
| Ekman / Prandtl | Ek, Pr | 3 × 10⁻⁴, 1 | exact from the values above |
| Rayleigh | Ra | 6.67e5 → 6.67e7 | set via ΔT (Table 2 of paper) |
| ΔT ↔ Ra | | 9.73e-3 / 9.73e-2 / 0.973 K | Ra = 6.67e5 / 6.67e6 / 6.67e7 |

Check: Ra = α·g·ΔT·D³/(νκ) = 6.85e7·ΔT[K] ✓ reproduces the paper's Table 2 with these dimensional values.

**Forcing:** top and bottom grid volumes relaxed toward fixed temperatures over 2 timesteps — "in effect constant temperatures" (Dirichlet ΔT). Initial random T perturbations ~ 10⁻³ΔT. Convergence identified via global-mean kinetic energy.

**Drag:** τ_h = ρ₀·cD·|u_h|·u_h at top and bottom; viscosity deliberately *absent* from the drag law so the boundary stress stays ~100× smaller than a no-slip stress (paper §2.4).

**Diagnostics to match:** time/zonal-mean zonal jets (Figs 2–3), surface zonal stress τ_φ = ρ₀cD|u|u_φ, axial torque from the surface stress,
T_z = 2πR³ ∫ τ̄_φ cos²θ dθ (eq. 5), and the mean/turbulent KE decomposition (eqs. 7–9).

## 3. Solving stack (Held-Suarez-style custom build)

One self-contained `europa-ocean.c`, structured like `held-suarez.c`:

```
grid/multigrid.h      spherical.h
layered/hydro.h       (free surface, layered momentum/continuity)
layered/implicit.h    (implicit barotropic mode)
layered/dr.h          (buoyancy from a temperature tracer + user EOS)
layered/coriolis.h    (F0() + custom K0() — used for drag, see §5)
layered/diffusion.h   (vertical viscosity/diffusion, implicit)
layered/remap.h       (sigma remapping each step)
layered/perfs.h, profiling.h
```

- **Grid:** global lon–lat on the multigrid: `dimensions(2,1); size(360); origin(-180,-90); periodic(right);` poles closed as dry boundary rows (as in `ocean.h`/`global-tides.c`). N = 256 first (Δ ≈ 1.4° ≈ 24 km at equator).
- **Layers:** nl = 10–20 σ-layers kept equidistant by `remap.h`; flat bottom at −100 km; deformable free surface with implicit barotropic mode (EK… not needed: gravity-wave speed √(gD) ≈ 360 m/s handled implicitly).
- **Time:** DT ≈ 300–600 s; advective CFL is trivial (jet speeds ≪ grid speed); one Europa day = 510 steps at DT = 600 s.
- **Coriolis:** traditional approximation, `F0() = 2*Omega*sin(y)` with y the latitude.

## 4. Thermal forcing (constant ΔT, bottom hot / top cold)

Closest mapping of the paper's Dirichlet boundary temperatures, mirroring its "relax over 2 timesteps" implementation:

1. **Tracer + EOS:** `dr.h` carries a per-layer temperature tracer T; linear EOS
   `#define drho(T) (-2e-4*((T) - T0))` (dimensionless Δρ/ρ₀; warm ⇒ buoyant).
2. **Boundary relaxation event** (`i++`), top and bottom layers only:
   `T[] = (T[] + dt/tau*T_bc)/(1 + dt/tau)`, τ = 2·DT, T_bc = T0 ± ΔT/2 — verbatim MITgcm semantics.
3. **Vertical thermal diffusion** κ = 61.6 m²/s through `vertical_diffusion()` with zero-flux BCs at both ends — heat enters the interior only via the pinned boundary layers. Implicit ⇒ unconditionally stable.
4. **Convective adjustment event:** hydrostatic layers cannot overturn; after diffusion, mix T (thickness-weighted) across statically unstable adjacent layers. This is the convection surrogate — the key difference from the paper's resolved plumes.
5. **Init:** T = T0 + random noise ±10⁻³·ΔT (paper's seed); statistically steady state detected from global KE, restarts via `dump()`.
6. **ΔT per run** from the paper's Table 2; start with **Ra = 6.67e6 (ΔT = 0.0973 K)**.

## 5. Quadratic drag at bottom **and** top

`layered/coriolis.h` applies a per-layer linear coefficient K0; making it velocity-dependent produces quadratic drag (pattern proven in `global-tides.c` / `ocean.h`, which even use the paper's Cb = 2×10⁻³ — but bottom-only):

```c
const double Cd = 2e-3;
#define K0() ((point.l == 0 || point.l == nl - 1) ? \
        (h[] > dry ? Cd*norm(u)/h[] : HUGE) : 0.)
#include "layered/coriolis.h"
```

This is exactly τ = ρ₀cD|u|u recast as a per-unit-mass deceleration of the boundary layer (cD|u|u/h), applied semi-implicitly (α_H = 1). The "near-boundary parcel velocity" of the paper is the top/bottom **layer-mean** velocity — the layer-model analog of MITgcm's boundary cell.

**Companion setting (essential, per paper §2.4):** vertical viscosity must be **free-slip at both boundaries** (`lambda_b` ≈ HUGE at the bottom; the top is already free-slip by default) so the drag law is the *only* boundary stress. Otherwise the model slides back toward a no-slip regime the paper explicitly avoids.

## 6. Diagnostics & validation

- Top-layer (zonal-, time-mean) zonal velocity ū(φ): compare to Figs 2–3; check superrotation and its confinement outside θ_tc = arccos(1 − D/R) ≈ 20.6° at low Ra, reversal at high Ra.
- Surface zonal stress τ_φ = ρ₀cD|u_φ||u| on the top layer; axial torque
  `Tz = Σ_cells τ_φ · (R cosφ) · dv()` (paper eq. 5), as a time series; mean/turbulent KE split (ū vs u′) via Welford running averages (held-suarez pattern), feeding eqs. 7–9.
- Movies (`output_ppm`): top-layer u.x, T, vorticity; KE time series for the steady-state criterion; dumps for restart.
- Ra sweep **after** the gate: 6.67e5 / 6.67e6 / 6.67e7; torque-vs-Ra against Fig 4, including the sign reversal near Ra ≈ 1.5e7.

## 7. Development stages

| # | Stage | Content | Gate |
|---|---|---|---|
| 1 | Scaffold | write `europa-ocean.c`; N=64 smoke test (CPU), conservation checks | compiles, runs, conserves |
| 2 | Control | short runs: rotation off vs on; drag on/off | sensible qualitative response |
| 3 | Production 1 | **GPU**, N=256, nl = 10–20, DT ≈ 600 s, Ra = 6.67e6, run to KE plateau with dumps | **jets emerge? torque sign?** ← go/no-go |
| 4 | Sweep | Ra = 6.67e5, 6.67e7 (same config, restarted) | regime transition reproduced |
| 5 | Refine | N=512 check; longer averaging | robustness of stage 3–4 results |
| 6 | Later | `layered/nh.h` non-hydrostatic extension; plume-scale sensitivity | comparison against stage 3 |

## 8. Honest limitations (fixed expectations)

- Hydrostatic + traditional Coriolis vs the paper's non-hydrostatic full-Coriolis: plume-scale convection and Taylor columns are **parameterized (convective adjustment), not resolved**; with D = 100 km and Δx ≈ 24 km the hydrostatic assumption is formally weak at grid scale — hence stage 6.
- The tangent-cylinder discontinuity is a full-sphere columnar-geometry effect; a thin-shell layered model may reproduce it only partially.
- **Main scientific risk:** whether zonal jets spin up at all from uniform ΔT forcing in this framework — addressed head-on by the stage-3 gate.
- Target is **regime similarity** (jet structure, torque sign/regime with Ra), not quantitative torque magnitudes.
- MITgcm's constant boundary T reproduces Dirichlet exactly; our τ = 2·dt relaxation is the same in spirit (paper used 2 timesteps itself).

## 9. Deliverables (on approval)

1. `europa-ocean.c` — the model, in this directory (`sandbox/hugoj/Hay_2023/`).
2. Build/run notes for CPU smoke tests and GPU production (`qcc`/make, `CFLAGS` flags).
3. Post-processing (zonal means, stress, torque time series) and a comparison note against paper Figs 2–4.

**Next action on approval:** start stage 1 (write `europa-ocean.c`, compile, smoke test).
