# Debye-Hückel e-e Screening Implementation

## STATUS: COMPLETE

The production Debye e-e path has been exercised by a real JAC atomic calculation
(`Atomic.Computation` -> `Basics.perform`, He / Be / Ne, verified by a caller census AND a sampling profile of
the running calculation), and the unscreened regression is clean (**83 / 83**, one test added, no approved
reference re-approved). Both conditions set for COMPLETE are met.

**THE PERMANENT TECHNICAL RECORD IS `docs/debye-electron-electron-screening.md`.** That report -- 21 sections,
~2180 lines -- carries the theory, the physics-to-code mapping, the forensic record of every modification, the
full validation tables and the limitations, with every claim labelled SOURCE-VERIFIED / ANALYTICALLY VERIFIED /
NUMERICALLY VERIFIED / INFERENCE / NOT VERIFIED. **THIS FILE STAYS THE CHRONOLOGICAL DIARY** and is deliberately
not rewritten to match it: the two serve different readers, and the order in which things went wrong here is
information the cleaned report cannot carry.

Two genuine defects were found by review and fixed, neither of which had ever affected a published number:
an inverted `(-1)^(k+1)` in the closed-form branch of `RadialIntegrals.iScaledYukawa` (Debugging history,
item 7), and a rank-blind overflow floor that returned `Inf` at k = 12 and `NaN` at k >= 14 (item 8, added
05-Sep-2026 in the final review).

All changes are UNCOMMITTED, per CLAUDE.md.

## Baseline

```
Repository:        <JAC-root>   (OpenJAC/JenaAtomicCalculator.jl)
Branch:            master
Starting commit:   75f963acac367e7365eed27854199a572a465d2c
Julia version:     1.12.7   (juliaup)
Platform:          Linux 6.8.0-138-generic, x86_64
Working tree:      clean except untracked CLAUDE.local.md (a user file -- PRESERVED, not touched)
Baseline test:     cd test && julia --project=.. runtests.jl        [see "Tests" below]
```

Note on Rule 11 (CLAUDE.md): the suite MUST be started from `test/`, never from the repository root.

## Repository architecture (Phase 1 reconnaissance)

All names below were read out of `src/`; nothing here is invented.

| Physics object | JAC representation |
|---|---|
| relativistic orbital | `Radial.Orbital` with tabulated `P::Vector`, `Q::Vector`, `subshell::Subshell` |
| radial grid | `Radial.Grid` -- `grid.r` (points), `grid.wr` (quadrature weights), `grid.NoPoints`, plus B-spline data `nsL`, `nsS` |
| radial integration | `sum_i f[i]*grid.wr[i]`; `RadialIntegrals.V0`; per-cell 8-pt Gauss-Legendre `RadialIntegrals.cellIntegral` |
| radial pair density | `rho_ac(r) = P_a P_c + Q_a Q_c` (formed inline everywhere) |
| electron density | `rho_t(r) = sum_a occ_a (P_a^2 + Q_a^2)` (formed inline in `Basics.computePotential`) |
| Coulomb kernel `r_<^k / r_>^(k+1)` | inline in `RadialIntegrals.SlaterRk`, `SlaterRkComponent`, `Yk_ab`, `buildScreenedPotential`, `buildScreenedPotentialPair`, and in the `max(r_i,r_j)` loops of `Basics.computePotential` |
| screened potential `Y_k`/`V_k` | `RadialIntegrals.buildScreenedPotential(k, b, d, grid)` -- **O(N) forward/backward sweep** |
| Slater integral `R^k(abcd)` | `RadialIntegrals.SlaterRkKinkAware` (production), `SlaterRk` (naive O(N^2)), `SlaterRkReference` (oracle) |
| direct **and** exchange `F^k`/`G^k` | the SAME `R^k(abcd)`; direct vs exchange is decided purely by the *slot assignment* `(a,b,c,d)` in the angular coefficient `SpinAngular.Coefficient2p`. **One kernel covers both.** |
| angular coefficients | `SpinAngular.computeCoefficients(TwoParticleOperator(0,plus), csf_r, csf_s, subshells)` -> `coeff.V`, `coeff.nu`, `coeff.a..d` |
| two-electron matrix element | `InteractionStrength.XL_Coulomb(L,a,b,c,d,grid[,cache])` = angular `C^L` prefactor x `SlaterRkKinkAware` |
| CI Hamiltonian | `Hamiltonian.setupMatrix` / `setupMatrixKinkAware` -> `Hamiltonian.performCI` |
| SCF driver | `SelfConsistent.performSCF` -> `Basics.scfProcedure(scField)` -> one of `:meanFieldIteration`, `:averageLevel`, `:optimizedLevel`, `:hydrogenicStartOnly` |
| Hartree/direct SCF potential | `Basics.computePotential(scField, grid, level)` for `DFSField` (JAC's default), `CHField`, `KSField`, `HSField` |
| caching | `InteractionStrength.XLCache` (per CI matrix); `RadialIntegrals.ScreenedPotentialCache` (per SCF run, AL/EOL only) |
| configuration | `ManyElectron.AsfSettings` (`eeInteraction`, `eeInteractionCI`, `scField`, ...); `Basics.AbstractPlasmaModel`; module globals in `Defaults` |

### Coulomb e-e call graph (the actual one, traced)

```
AsfSettings.eeInteractionCI  ::CoulombInteraction
        |
SelfConsistent.performSCF(configs, nm, grid, settings)
        |-- (SCF)  scfProcedure(DFSField) == :meanFieldIteration
        |             -> SelfConsistent.solveMeanFieldBasis[Anderson]
        |                   -> Basics.computePotential(DFSField, grid, wLevel)      <-- k=0 DIRECT Hartree,
        |                        rho_t(r) built from the orbitals;  V(r)=int rho_t(r')/r_> dr'    O(N^2)
        |                   -> Basics.add(nuclearPotential, wp) -> Bsplines.setupLocalMatrix -> new orbitals
        |
        `-- (CI)   Hamiltonian.performCI -> Hamiltonian.setupMatrix
                      -> SpinAngular.computeCoefficients (angular)
                      -> InteractionStrength.XL_Coulomb(nu,a,b,c,d,grid,cache)      <-- DIRECT *and* EXCHANGE
                            -> RadialIntegrals.SlaterRkKinkAware(k,a,b,c,d,grid)
                                  -> RadialIntegrals.buildScreenedPotential(k,b,d,grid)   O(N) sweep
                                  -> contract V_k(r) against rho_ac(r) with grid.wr
```

**Critical architectural answer (section 7 of the task): JAC is Path A + Path C.**
* **Path A** -- one generic radial kernel (`buildScreenedPotential` / `SlaterRkKinkAware`) serves both direct
  and exchange Slater integrals; they differ only in the orbital slots the angular coefficient supplies.
* **Path C** -- a *separate*, independent k=0 direct-potential calculation
  (`Basics.computePotential`) drives the SCF. It does **not** call the Slater machinery at all; it re-derives
  `int rho(r')/r_>` with its own `max(r_i,r_j)` double loop. Screening the Slater path alone therefore leaves
  the orbitals unscreened -- exactly the failure mode the task warns about.
* Third, partly independent path: the B-spline `ScreenedPotentialCache` line
  (`buildScreenedPotentialPair` / `buildScreenedPotentialCache` -> `XL_CoulombTensor` / `XL_CoulombKinkAware`)
  used *only* by `scField = ALField()` and `EOLField()`.

### What already existed before this work

JAC already carried a partial Debye-Hückel vocabulary, and it is important not to mistake it for a finished path:

* `Basics.DebyeHueckelModel(debyeLength)` -- the parameter object (lambda_D in a_0). **Reused as-is.**
* `Nuclear.nuclearPotentialDH` -- Debye screening of the **nucleus**. Out of scope (task section 20); untouched.
* `InteractionStrength.XL_Coulomb_DH` -> `RadialIntegrals.SlaterRkDebyeHueckel` -- a screened Slater integral
  reached from `Basics.compute(JP, basis, nm, grid, settings, plasmaModel)` (the `Plasma.LineShiftScheme` CI
  step) and from `AutoIonization.computeLinesPlasma`.
  **It was mathematically the right expansion but truncated at p <= 2 and O(N_r^2)** -- see "Debugging history".
* `Plasma.perform(::LineShiftScheme, ...)` computes **field-free** (unscreened) orbitals and then re-does only
  the CI matrix under the plasma model. So the SCF path was never screened.

## Mathematical derivation

The e-e interaction is replaced by

    V_ee(r_12) = exp(-mu r_12) / r_12 ,        mu = 1 / lambda_D   [a_0^-1],  r in a_0, energy in Hartree.

Its multipole (Gegenbauer) expansion is

    exp(-mu r_12)/r_12  =  mu sum_k (2k+1) i_k(mu r_<) k_k(mu r_>) P_k(cos omega),

with the modified spherical Bessel functions normalised as

    i_0(x) = sinh(x)/x ,           k_0(x) = exp(-x)/x .

so the whole implementation is the single substitution, at fixed rank k,

    r_<^k / r_>^(k+1)   -->   mu (2k+1) i_k(mu r_<) k_k(mu r_>).

**Coulomb limit.** i_k(x) -> x^k/(2k+1)!! and k_k(x) -> (2k-1)!!/x^(k+1) as x -> 0, and (2k+1)!! = (2k+1)(2k-1)!!, so

    mu (2k+1) i_k(mu r_<) k_k(mu r_>)  ->  (2k+1) (2k-1)!!/(2k+1)!! * r_<^k/r_>^(k+1)  =  r_<^k/r_>^(k+1).   [exact]

**The screened potential.** The Coulomb sweep JAC already uses,

    V_k(r) = r^-(k+1) int_0^r s^k rho(s) ds  +  r^k int_r^inf s^-(k+1) rho(s) ds,

becomes, with the SAME separable two-region structure (this is what keeps the algorithm O(N)):

    V_k^Y(r) = mu(2k+1) [ k_k(mu r) int_0^r i_k(mu s) rho(s) ds  +  i_k(mu r) int_r^inf k_k(mu s) rho(s) ds ].

**Weak-screening check.** exp(-mu r12)/r12 = 1/r12 - mu + mu^2 r12/2 - ... so to first order every electron
PAIR contributes -mu, and for N electrons  Delta E_ee ~ -N(N-1)/2 * mu.

## Numerical stability

`i_k(x)` grows like exp(x)/(2x) and `k_k(x)` falls like exp(-x)/x, so neither may be formed on its own over a
grid that reaches mu*r ~ 10^2..10^5. The implementation works throughout in the SCALED functions

    ihat_k(x) = exp(-x) i_k(x)          khat_k(x) = exp(+x) k_k(x)

both of which are O(1/x) for large x, and it carries the residual exponential as a cell-local factor
`exp(-mu (r_> - r_<))` which is always <= 1. See `RadialIntegrals.iScaledYukawa` / `kScaledYukawa` /
`yukawaKernel` / `buildScreenedPotentialYukawa` for the details and the branch thresholds.

## The user-facing API

```julia
using JenaAtomicCalculator

# ON.  The argument is the DEBYE LENGTH lambda_D, in a_o (Bohr radii).  1/r_12 becomes exp(-r_12/lambda_D)/r_12.
Defaults.setDefaults("e-e screening", Basics.DebyeHueckelModel(5.0))

wa = Atomic.Computation(Atomic.Computation(), name="Be", grid=Radial.Grid(true),
                        nuclearModel=Nuclear.Model(4.), configs=[Configuration("1s^2 2s^2")])
wb = perform(wa; output=true)

# OFF -- the default, and JAC's behaviour before this switch existed.
Defaults.setDefaults("e-e screening", Basics.NoPlasmaModel())

Defaults.getDefaults("e-e screening")   # the model currently in force
Defaults.eeScreeningMu()                # mu = 1/lambda_D [a_o^-1], EXACTLY 0.0 when unscreened
```

**Units, stated once and unambiguously.** `lambda_D` is a LENGTH in `a_o`; `mu = 1/lambda_D` is an INVERSE
length in `a_o^-1`; `mu*r` is dimensionless; energies are Hartree. The only place the two are easy to confuse
is the pre-existing `RadialIntegrals.SlaterRkDebyeHueckel(..., lambda)`, whose argument named `lambda` has
always in fact been `mu` -- its callers pass `1/plasmaModel.debyeLength`. That misnomer was NOT renamed (it is
a public signature and renaming it is outside this task) but is now stated in the first lines of its docstring.

**What the switch does and does not reach.** It reaches the two-electron Slater integrals of the Hamiltonian,
direct and exchange alike, AND the k = 0 direct potential that drives the mean-field SCF. It does NOT touch the
electron-nucleus potential; nuclear Debye screening remains the separate `Nuclear.nuclearPotentialDH`.

**Invalid input is refused where it is set**, not carried into the kernels: `lambda_D <= 0`, `NaN` and `Inf`
all raise from `setDefaults`, and any plasma model other than `NoPlasmaModel`/`DebyeHueckelModel` raises too.

## Design decision

**The smallest change that reaches both Coulomb paths with ONE authoritative mu.**

Four candidate architectures were listed in the task; the code answered the question itself (see the call graph
above). JAC is **Path A + Path C**, so a minimal implementation needs exactly TWO insertion points and no more:

| | insertion point | reaches |
|---|---|---|
| A | `RadialIntegrals.buildScreenedPotential` | every kink-aware Slater integral -- DIRECT and EXCHANGE alike, since the two differ only in the orbital slots the angular coefficient supplies |
| C | `Basics.computePotential` (5 methods) | the k = 0 direct potential the mean-field SCF iterates the orbitals against |

Three decisions were taken and are worth stating because each had a plausible alternative.

1. **The screening parameter is a `Defaults` global, not a field of `AsfSettings`.** Neither
   `buildScreenedPotential(k, b, d, grid)` nor `computePotential(scField, grid, level)` is handed a settings
   object, and neither is called from anywhere holding one. Threading a parameter to both would touch eight
   modules and change the positional `AsfSettings` constructor that ~200 call sites use. Against that: a run
   whose SCF and whose Hamiltonian disagreed about mu would be silently wrong. One authoritative value, read
   only through `Defaults.eeScreeningMu()`, is the smaller risk -- and is how JAC already treats
   `GBL_FRAMEWORK` and `GBL_CONT_POTENTIAL`. The reasoning is recorded in `src/module-Defaults.jl` beside the
   global, not only here.
2. **The unscreened case takes the ORIGINAL code path, not an equivalent one.** Every branch tests
   `mu == 0.` or `mu > 0.` for identity, so with screening off no Yukawa code is entered at all and the former
   numbers are reproduced bit for bit rather than to within a tolerance. This is what makes the regression
   statement below a statement about code and not about quadrature.
3. **Where the screening cannot reach, the code REFUSES rather than half-screens.** A partly screened
   calculation converges, returns plausible energies and says nothing about the half that was left bare. The
   two such places are named in `SelfConsistent.checkEeScreeningIsConsistent`.

## Implementation changes

Seven files, 849 insertions and 44 deletions at the close (`git diff --numstat`). Per the format asked for in
the task:

### `src/module-Defaults.jl`  (+74 / -1)

| item | role | change | why | side effects | tests |
|---|---|---|---|---|---|
| `GBL_EE_SCREENING` | new global | added, `= Basics.NoPlasmaModel()` | the single authoritative screening state | none while it holds the default | (5) |
| `setDefaults(::String, ::Basics.AbstractPlasmaModel)` | new method | validates `isfinite(lambda) && lambda > 0`, refuses any model but `NoPlasmaModel`/`DebyeHueckelModel`, prints a notice | invalid parameters must be rejected where they are set, not carried into the kernels as Inf/NaN | a new `setDefaults` method; the key `"e-e screening"` was unused | (8) |
| `eeScreeningMu()` | new function | returns `1/lambda_D`, or EXACTLY `0.0` | one accessor, so no caller re-derives mu | newly exported name | census |
| `getDefaults("e-e screening")` | new branch | returns the model | symmetry with every other Defaults key | none | -- |

### `src/module-RadialIntegrals.jl`  (+351 / -35)

| item | role | change | why | side effects | tests |
|---|---|---|---|---|---|
| `GBL_YUKAWA_X_SWITCH = 35.0` | new const | branch between the ascending series and the closed form of `ihat_k` | neither form is accurate over the whole range | see the warning now written beside it | (8c) |
| `GBL_YUKAWA_X_FLOOR = 1.0e-25` | new const | below it the Coulomb kernel is returned unchanged | `ihat_k ~ x^k` and `khat_k ~ x^-(k+1)` are formed separately and only their PRODUCT is well scaled; NOT a physical clamp -- at x = 1e-25 the two kernels differ by <= 1e-25 relative, 1e-9 times the Float64 epsilon | none | (3) |
| `iScaledYukawa(k, x)` | new | `exp(-x) i_k(x)` | `i_k` alone overflows past x ~ 709 | -- | (1), (8a), (8c) |
| `kScaledYukawa(k, x)` | new | `exp(x) k_k(x)`, an exact finite all-positive sum | `k_k` alone underflows | -- | (1), (8a) |
| `yukawaKernel(k, r_<, r_>, mu)` | new | `mu (2k+1) i_k(mu r_<) k_k(mu r_>)` | THE substitution | asserts `r_< <= r_>`; the kernel is not symmetric | (2), (3), (4), (8) |
| `buildScreenedPotentialYukawa` | new | the O(N) forward/backward sweep, in scaled form | keeps the separable structure; an EXPONENTIAL tail beyond the source extent, not the Coulomb power law | -- | (7), 9(c), 10 |
| `SlaterRkYukawa` | new | contracts that potential against rho_ac | -- | -- | (5), (6) |
| `buildScreenedPotential` | **the one branch point** | 4 lines added at the top | every kink-aware consumer funnels through it, so direct and exchange are both reached | none at mu = 0: the added test is an identity test | (5), census |
| `SlaterRkDebyeHueckel` | pre-existing | body replaced by a call to `SlaterRkYukawa` | the old body truncated the `i_k` series at p <= 2 (90 % error at mu r_< = 10) and was O(N_r^2) | callers of `Plasma.LineShiftScheme` and `AutoIonization.computeLinesPlasma` get DIFFERENT, more accurate numbers -- see Known limitations | -- |
| `Yk_ab(k, r, rho, mtp, grid, mu)` | new method | screened counterpart, falls through to the unscreened method at mu <= 0 | `HSField` builds its mean field from `Yk_ab` | none | (6) |

### `src/module-BasicsAZ-inc-compute.jl`  (+65 / -6)

`Basics.computePotential` for **`CHField`, `HSField`, `KSField`, `DFSField(level)`, `DFSField(basis)`** -- the
five e-e mean fields -- now read `mu = Defaults.eeScreeningMu()` and replace `1/r_>` by
`yukawaKernel(0, r_<, r_>, mu)` (`HSField` instead passes `mu` to `Yk_ab`). `mu == 0.` keeps the original line.

**`AaDFSField` and `AaHSField` are deliberately untouched** -- they are the average-atom plasma line, out of
scope, and their own `mu` argument is the CHEMICAL POTENTIAL. **`ThomasFermiField` is untouched** -- it is a
statistical starting potential built from Z and the electron number, not a Hartree term from the density, and
Thomas-Fermi screening is explicitly out of scope.

### `src/module-InteractionStrength.jl`  (+12 / -2)

Both `XLCache` keys now carry `Defaults.eeScreeningMu()` in the `Float64` slot that previously held a constant
`0.`. Without it a value computed under one screening model could be handed back under another; a plasma-shift
run computes field-free first and screened second in one session. At mu = 0 the key is byte-for-byte the former
one.

### `src/module-SelfConsistent.jl`  (+58)

`checkEeScreeningIsConsistent(settings)` -- returns immediately at mu = 0, otherwise raises a naming error for
`ALField`/`EOLField` (their B-spline tensor line is NOT reached by the switch, so their SCF would run bare
while the Hamiltonian was screened) and for `BreitInteraction`/`CoulombBreit`/`CoulombGaunt` (the transverse
part has its own kernels; a screened Coulomb term beside a bare Breit term is not a defined approximation).
Called at both `performSCF` entry points.

### `src/module-TestFrames-inc-plasma.jl` and `test/runtests.jl`  (added by this audit)

`testModule_PlasmaScreening()` and its entry in the previously empty `@testset "JAC plasma"`. See Tests.

## Tests

Every row below was EXECUTED; the "Executed?" column exists so that nothing can be read as a plan. Diagnostics
live in `tools/diag-debye-screening-{kernel,atomic,performance,production,stability}.jl` and are re-runnable.

### The empty plasma testset -- what it actually was

Not an oversight and not a placeholder. `src/module-TestFrames-inc-plasma.jl` carried a note: the former
`testModule_PlasmaShift` was **deleted on 09-Aug-2026** because PlasmaShift stopped being an `Atomic` property
and became a `Plasma.Computation` scheme, "so a test for it belongs with the Plasma module and has to be
written afresh". `test/approved/test-PlasmaShift-approved.sum` is the orphan of that deletion and is referenced
by nothing (SOURCE-VERIFIED by grep over `src/` and `test/`). So: a production-level Debye regression test was
genuinely MISSING, and this is the right place for it. `testModule_PlasmaScreening` was added there and wired
into the previously empty `@testset "JAC plasma"`.

It deliberately uses **no approved data and cannot**: every check is an exact closed form, an identity between
two independent constructions, or the mu -> 0 limit, so it cannot pass on a stale stored reference. The orphaned
`.sum` was left untouched -- deleting it is an editorial act for the maintainer, not for this task.

### Layered results

| # | Level | What it checks | Executed? | Result |
|---|---|---|---|---|
| 1 | 1 Bessel | `ihat_k`, `khat_k` against `SpecialFunctions` (`besselix`/`besselkx`, with the (pi/2) conversion) | YES | 2.56e-14 (ihat), 8.93e-16 (khat) |
| 2 | 1 Bessel | against an independent **256-bit BigFloat** reference (series below x = 60, upward recurrence above), 150 (k, x_<, x_>) combinations | YES | worst 2.90e-16 excluding Float64 underflow |
| 3 | 1 Bessel | continuity/accuracy across the internal branch switch at x = 35, each branch against the 256-bit reference at its own x | YES | series <= 8.81e-16, closed form <= 3.28e-16, k = 0..8 |
| 4 | 2 kernel | **Test C** -- kernel against an angular Legendre projection `(2k+1)/2 int P_k(w) exp(-mu r12)/r12 dw` (QuadGK; uses NO Bessel function) | YES | ~1e-15 |
| 5 | 2 kernel | k = 0 against the closed form `sinh(mu r_<) exp(-mu r_>)/(mu r_< r_>)` | YES | <= 2.06e-16 |
| 6 | 2 kernel | **Test B** -- Coulomb limit: exact at mu = 0, and the OBSERVED order of approach | YES | p = 0.996-0.999 (k = 0), 1.995-2.000 (k >= 1) |
| 7 | 2 kernel | **Test E** -- sign and the Coulomb bound, 7830 combinations, k = 0..8, mu*r = 1e-30..1e6 | YES | **0 violations** |
| 8 | 3/4 radial | O(N) sweep vs an explicit O(N^2) double sum over the same kernel | YES | 5.7e-6 .. 1.6e-4, the same magnitude at mu = 0 (quadrature, not physics) |
| 9 | 5 two-electron | **Test F** -- He 1s2s singlet-triplet splitting = 2 G^0(1s,2s) exactly, i.e. pure EXCHANGE | YES | moves 3.4558e-2 -> 7.3395e-2 Ha |
| 10 | 5 two-electron | direct F^0, exchange G^0 and exchange G^1 each respond to the switch | YES | all three move by > 1e-6 relative |
| 11 | 6 SCF | **Test G** -- screened SCF+CI vs screened CI on Coulomb-converged orbitals | YES | see "SCF consistency" below |
| 12 | 7 atomic | **Test D** -- `Delta E = -mu N(N-1)/2` at mu = 1e-5 | YES | He 0.999982, Be 0.999987, Ne 0.999994; all negative |
| 13 | 7 atomic | **Test H** -- grid convergence, compact and diffuse | YES | Be ~1e-9 Ha (rbox 25 -> 40); diffuse Li 1s^2 3s ~1e-7 Ha (rbox 70 -> 110) |
| 14 | 7 atomic | **Test I** -- progressively stronger screening | YES | finite, no NaN/Inf to mu = 200 (lambda_D = 0.005 a_o); saturates at -19.8897 Ha |
| 15 | -- | **Test A** -- full suite with screening OFF | YES | see "Regression" below |
| 16 | -- | guards: 5 refusals + 4 invalid lambda_D + a mis-ordered kernel call | YES | all refuse; the same settings are still ACCEPTED with screening off |

### Production atomic validation (section 2 of the audit) -- NUMERICALLY VERIFIED

Run through the user-facing route `Atomic.Computation` -> `Basics.perform`, grid 1582 points, rbox 25 a_o.

```
                                    N    E(Coulomb) [Ha]  E(screened) [Ha]   Delta E [Ha]   -mu N(N-1)/2
  mu = 0.20 a_o^-1  (lambda_D = 5 a_o)
  He  1s^2                          2       -2.857967131      -3.034819576      -0.176852      -0.200000
  Be  1s^2 2s^2                     4      -14.570977349     -15.533391190      -0.962414      -1.200000
  Ne  1s^2 2s^2 2p^6               10     -128.672259270    -136.690308921      -8.018050      -9.000000
  mu = 1.00 a_o^-1  (lambda_D = 1 a_o)
  He  1s^2                          2       -2.857967131      -3.444526208      -0.586559      -1.000000
  Be  1s^2 2s^2                     4      -14.570977349     -17.390701542      -2.819724      -6.000000
  Ne  1s^2 2s^2 2p^6               10     -128.672259270    -156.961366968     -28.289108     -45.000000
```

The compact single-pair case (He) and two cases with many pairs (Be, 6; Ne, 45) were required by the audit. The
last column is the FIRST-ORDER law and is printed as a scale, not a tolerance: at mu = 0.2 and 1.0 the typical
r_12 is of order 1 a_o, so mu r_12 is not small and the higher orders are meant to be visible. That the shift
is a monotonically decreasing FRACTION of -mu N(N-1)/2 as mu grows (88 %, 80 %, 89 % at mu = 0.2 against 59 %,
47 %, 63 % at mu = 1.0) is the expected behaviour of a series in mu r_12, not a defect.

The mu -> 0 limit ON THE SAME PRODUCTION ROUTE (Be, N = 4, first-order shift -6 mu):

```
    mu [1/a_o]               E [Ha]       Delta E [Ha]   DeltaE/(-6mu)     residual
   0 (Coulomb)     -14.570977349023                  -               -            -
       1.0e-02     -14.630234506802  -5.9257157779e-02     0.987619296    7.428e-04
       1.0e-03     -14.576969843895  -5.9924948720e-03     0.998749145    7.505e-06
       1.0e-04     -14.571577273895  -5.9992487239e-04     0.999874787    7.513e-08
       1.0e-05     -14.571037348248  -5.9999224945e-05     0.999987082    7.751e-10
       1.0e-06     -14.570983349030  -6.0000073940e-06     1.000001232   -7.394e-12
```

The residual falls by exactly two decades per decade of mu -- second order, as the theory requires -- until it
reaches the SCF convergence floor at mu = 1e-6, where it changes sign. **The Coulomb limit is recovered.**

### The tested path IS the production path (section 1) -- NUMERICALLY VERIFIED

Two independent instruments, neither of which changes a number.

**(a) A caller census.** `Defaults.eeScreeningMu()` was replaced by a version returning the identical value and
additionally recording its caller from `stacktrace()`. Over the six production runs above:

```
  XL_Coulomb                        148 calls      <- the CI Hamiltonian
  computePotential                  108 calls      <- the SCF direct potential
  #buildScreenedPotential#1          88 calls      <- the Slater kernel
  checkEeScreeningIsConsistent       12 calls      <- the consistency guard
```

**(b) A sampling profile of one real `Basics.perform` run** (Be, two configurations, mu = 0.2):

```
  computePotential               1264 samples        SlaterRkKinkAware             195
  yukawaKernel                   1202 samples        buildScreenedPotential        195
  iScaledYukawa                   977 samples        buildScreenedPotentialYukawa  189
  kScaledYukawa                    63 samples
```

The Yukawa kernels are executed inside a real atomic calculation. Together with the SOURCE-VERIFIED chain
`Basics.perform` (module-BasicsAZ-inc-perform.jl:21) -> `SelfConsistent.performSCF` ->
{`Basics.computePotential`; `Hamiltonian.performCI` -> `InteractionStrength.XL_Coulomb` ->
`RadialIntegrals.SlaterRkKinkAware` -> `buildScreenedPotential`}, the chain is closed at both ends.

### Direct AND exchange (section 5) -- SOURCE-VERIFIED and NUMERICALLY VERIFIED

The architecture SHARES one kernel, and this was established by reading the code, not assumed:
`InteractionStrength.XL_Coulomb` computes an angular `C^L` prefactor and then calls
`RadialIntegrals.SlaterRkKinkAware(L, a, b, c, d, grid)` -- **the same call for both**. Direct and exchange
differ only in which orbitals the angular coefficient `SpinAngular.Coefficient2p` puts in the four slots:
`(a,b,a,b)` is `F^k`, `(a,b,b,a)` is `G^k`. There is no second radial routine to modify.

Confirmed numerically three ways: (i) `F^0(1s,2s)`, `G^0(1s,2s)` and `G^1(1s,2p_3/2)` all move under the switch
(new suite test); (ii) **Test F**, the He 1s2s singlet-triplet splitting, which is `2 G^0(1s,2s)` and contains
NO direct part, moves from 3.4558e-2 to 7.3395e-2 Ha -- an implementation that screened only `F^k` would leave
it exactly unchanged; (iii) the Ne production run, whose 45 pairs include `G^1(1s,2p)` and `G^2(2p,2p)`.

### SCF consistency (section 6) -- NUMERICALLY VERIFIED

JAC's default `scField = DFSField()` runs `scfProcedure == :meanFieldIteration`, i.e. it genuinely re-optimises
the orbitals; this is not a frozen-orbital implementation. The measurement below screens EVERYTHING in column 2
and screens only the CI matrix on Coulomb-converged orbitals in column 3 -- the latter is what a frozen-orbital
implementation would return.

```
          mu        full (SCF+CI)     CI only (frozen)       difference   % of Delta E
       0.010     -14.630234506802     -14.630219463692       -1.504e-05         0.0254
       0.100     -15.104234227457     -15.102385763470       -1.848e-03         0.3466
       0.500     -16.458003379469     -16.376148172566       -8.186e-02         4.3378
       1.000     -17.390701542441     -17.135405589875       -2.553e-01         9.0539
```

The orbital-relaxation part is second order in mu, as it must be, and reaches 9 % of the whole shift at
mu = 1 a_o^-1. Screening the Slater integrals alone would therefore have been wrong by that much -- which is
exactly the Path C failure the task warned about, and it is real here rather than hypothetical.

### Numerical stability (section 8) -- NUMERICALLY VERIFIED

* **Accuracy**: worst relative error 2.90e-16 against a 256-bit reference over 150 combinations spanning
  mu*r_> = 1e-12 .. 1e5 and r_</r_> = 1e-6, 0.5, 1.0, for k = 0, 1, 2, 4, 8.
* **Finiteness and sign**: 7830 combinations (k = 0..8, mu*r_> = 1e-30..1e6, six ratios). **Zero** NaN, zero
  Inf, zero negative values, zero values exceeding the Coulomb kernel.
* **Underflow**: 459 of those 7830 return an exact `0.0` where the Coulomb kernel is finite. These are
  `exp(-(x_> - x_<))` underflowing at a separation of more than 745 Debye lengths, where the true value is
  below 1e-323. Zero is the correct Float64 answer there, not a failure -- and it is the ONLY place the two
  disagree.
* **Cancellation**: `khat_k` is a finite all-positive sum, so none is possible. `ihat_k` below x = 35 is an
  all-positive ascending series; above 35 it is an alternating closed form whose terms are dominated by the
  first. Measured at the crossover: 8.81e-16 (series) and 3.28e-16 (closed form) against 256-bit.
* **No tolerance was widened and no clamp added to make anything pass.** The single threshold,
  `GBL_YUKAWA_X_FLOOR = 1e-25`, exists because `ihat_k ~ x^k` and `khat_k ~ x^-(k+1)` are formed separately and
  only their product is well scaled; at x = 1e-25 the two kernels differ by <= 1e-25 relative, nine orders
  below the Float64 epsilon.

### Caching (section 9) -- NUMERICALLY VERIFIED

| cache | lifetime | carries mu? | verdict |
|---|---|---|---|
| `InteractionStrength.XLCache` | per CI matrix, but a plasma run builds two in one session | **now YES** (added) | one cache object across a mu change: 4.208057097553249 -> 3.346223116076047 -> 4.208057097553249, `===` identical on return, 2 entries. Without mu in the key the third value would have repeated the first. |
| `vkCache` (V_k vectors, keyed `(k, subshell_b, subshell_d, mtp)`) | **per call by contract** -- keyed on subshell LABELS, so it already may not outlive the orbitals | not needed | exact within one mu at mu = 0 and mu = 0.7; a mu change forces re-convergence, which by that same contract forces a new cache |
| `ScreenedPotentialCache` / `tensorCaches` (B-spline tensor line) | per SCF | n/a | AL/EOL only, and both are REFUSED under screening |
| `rkCache`, `radial1pCache`, `radial2pCache` | per call | n/a | EOL only; refused |
| `Dierckx.Spline1D` | per build | n/a | interpolates the DENSITY rho_bd, carries no kernel; the Yukawa sweep rebuilds it identically |

**The asymptotic tail was the real risk and is handled.** The Coulomb branch continues `V_k` beyond the source
extent as `fullInner/r^(k+1)` -- a POWER LAW, which under screening would restore precisely the long-range
interaction the screening removes. The Yukawa branch continues it as `khat_k(mu r) exp(-mu r)` instead.
Measured on a deliberately truncated source, r = 7.243 .. 20.127 a_o beyond an extent of 7.236 a_o, mu = 0.8:

```
  apparent power-law exponent  d ln V / d ln r =  -11.085     (Coulomb k = 0 would give -1.000)
  chord decay constant         d ln V / d r    =   -0.879
  predicted for exp(-mu r)/(mu r)              =   -0.879
```

### Performance (section 10) -- NUMERICALLY VERIFIED

Minimum of 15 repetitions per point (a single sweep is milliseconds and one run in four catches a garbage
collection; the mean gave a spurious ratio of 0.74 at one size before this was fixed).

```
    N points   Coulomb [ms]    Yukawa [ms]    ratio
         812          2.888          3.254     1.13
        1330          6.812          7.630     1.12
        2338         18.619         20.777     1.12
        4354         61.059         66.216     1.08
```

**The ratio is flat over a 5.4x range of grid sizes.** That is the audit's question: the screened sweep has the
SAME complexity as the Coulomb one, so no O(N_r^2) has been introduced -- the separable two-region structure is
preserved and the algorithm is still one forward and one backward pass over the cells. (The per-point cost of
BOTH grows with N on this machine, identically; that belongs to the pre-existing unscreened algorithm and is
outside this task.)

The retired truncated `SlaterRkDebyeHueckel`, for comparison, was O(N_r^2) and cost **13.2x** the screened
sweep at 1330 points (13.983 ms against 1.062 ms).

End to end, on the production path (`SelfConsistent.performSCF`, Be 1s^2 2s^2, 1582 points):
**0.92 s unscreened against 1.77 s screened, i.e. 1.9x.** The cost sits in `Basics.computePotential`, whose
`i, j` double loop over the grid was ALREADY O(N^2) before this work; screening makes each of those N^2 kernel
evaluations more expensive but does not change the exponent.

### Regression with screening OFF (section 3 / Test A) -- NUMERICALLY VERIFIED

| run | when | invocation | result |
|---|---|---|---|
| baseline, before any edit | commit 75f963ac | `cd test && julia --project=.. runtests.jl` | **82 / 82**, 9m00.2s |
| after the implementation | -- | same | **82 / 82**, 7m42.6s |
| after the new plasma test | final | same | **83 / 83**, 7m42.8s |

The count rose from 82 to 83 because exactly one test was ADDED (`testModule_PlasmaScreening`, `[OK]`). Zero
Fail, zero Error, zero Broken in every run. **No `test/approved/*.sum` reference was modified** --
`git status --porcelain test/` reports only `test/runtests.jl`, which is this task's deliberate one-line
addition. Nothing was re-approved to make anything pass.

**On the invocation.** The audit asked for `julia --project=. -e 'using Pkg; Pkg.test()'`. That was NOT used,
and the deviation is deliberate: **CLAUDE.md Rule 11** requires the suite to be started from `test/`, because
every `TestFrames` include opens its summary with a relative path and a run from the repository root scatters
`test-*.sum` and `zzz-*.sum` across it. The Rule 11 form runs the identical `test/runtests.jl` with the
identical project.

## Debugging history

Kept because each entry cost real time and each would otherwise be rediscovered.

1. **`AmosException with id 2: overflow`** at `besseli(k+0.5, 1000.0)` in the prototype reference. Cause: an
   unscaled library call. Fixed by moving to `besselix`/`besselkx` and dropping the manual `exp` factors. This
   is the same lesson the implementation itself is built on.
2. **The Coulomb limit is NOT second order for all k.** The first expectation was that it would be. It is not:
   `mu i_0(mu a) k_0(mu b) = sinh(mu a) exp(-mu b)/(mu a b) ~ 1/b - mu`, so **k = 0 carries a constant -mu**
   and its RELATIVE error is O(mu r_>), i.e. order 1; only k >= 1 is order 2. That constant is precisely the
   physics behind `Delta E = -mu N(N-1)/2`. The test was then written to check the OBSERVED order rather than
   agreement at a fixed mu.
3. **A patch that would have been silently wrong, caught before it ran.** A `str.replace(old, new, 1)` on the
   `'rhot'` marker patched `Basics.computePotential(::AaDFSField, grid, orbitals, mu, temp)` FIRST -- and that
   method's `mu` is the CHEMICAL POTENTIAL. The inserted `mu = Defaults.eeScreeningMu()` shadowed it, and the
   method that should have been patched, `DFSField(basis)`, was left alone. Reverted and redone, then verified
   by listing the enclosing function of every `Defaults.eeScreeningMu()` occurrence.
4. **An apparent Test C failure that was the REFERENCE.** k = 0, r1 = 1, r2 = 1.0000001, mu = 0.7, relative
   1.54e-9. The QuadGK integrand is near-singular when r12 -> 1e-7. Resolved NOT by loosening the tolerance but
   by adding an independent exact closed form for k = 0, which agreed to 2.06e-16, and removing that point from
   the quadrature-reference list.
5. **The same trap again, in the new suite test.** `(1 - exp(-2x))/(2x)` as the reference for `ihat_0` at
   x = 1e-6 subtracts two numbers agreeing in eleven digits; the REFERENCE was wrong by 5.3e-12 and failed a
   correct implementation. Fixed with `-expm1(-2x)/(2x)`, again not by widening a tolerance.
6. **And a third time, in the 256-bit reference.** Its ascending series was capped at 20000 terms while it
   needs about x/2 of them, so at x = 1e5 it truncated and reported the implementation as infinitely wrong.
   Replaced above x = 60 by the upward recurrence
   `ihat_(k+1) = ihat_(k-1) - (2k+1)/x ihat_k` seeded by two exact closed forms.
7. **A real sign error in `iScaledYukawa`, found on 05-Sep-2026 by this audit.** The closed-form branch wrote
   `(-1)^(k+1)` as `isodd(k+1) ? 1. : -1.`, which is the inverse. Lifted out and evaluated BELOW its switch it
   was wrong by **3.7 % at x = 2**, 9.1e-5 at x = 5, 4.1e-9 at x = 10. **No result was ever affected**: the
   branch is only reached above x = 35, where the mis-signed term is `exp(-70) = 4e-31`, fifteen orders below
   the Float64 epsilon -- which is also why no test through the public API could have caught it, and why the
   finding is recorded as a warning beside `GBL_YUKAWA_X_SWITCH` rather than as a new test.
8. Smaller ones: `Radial.Grid(true; rnt=...)` takes only `printout` (use
   `Radial.Grid(Radial.Grid(false); ...)`); `Radial.Grid` has no `rbox` property (use `grid.r[end]`);
   `Basics.EOLField` and `CoulombBreit` do not have zero-argument constructors; Julia soft-scope warnings in
   the diagnostics, fixed with `Ref`; a `let a, b = f()` binds only `b`.

## Final validation

See "Tests" above for the numbers. Status labels:

| claim | label |
|---|---|
| one shared kernel serves direct and exchange | SOURCE-VERIFIED + NUMERICALLY VERIFIED |
| the SCF direct potential is screened with the same mu | SOURCE-VERIFIED + NUMERICALLY VERIFIED |
| the tested path is the path a real `perform()` takes | NUMERICALLY VERIFIED (census + profile) |
| the Coulomb limit is recovered, at kernel and at atomic level | ANALYTICALLY + NUMERICALLY VERIFIED |
| the electron-nucleus potential is untouched | SOURCE-VERIFIED (no edit outside the five e-e `computePotential` methods) + NUMERICALLY (at mu = 200 the Be energy saturates at -19.8897 Ha, the bare hydrogenic 4-electron sum) |
| no NaN/Inf/negative/over-Coulomb value over 7830 combinations | NUMERICALLY VERIFIED |
| no O(N_r^2) introduced | NUMERICALLY VERIFIED (flat ratio over 5.4x in N) |
| no Coulomb cache or tail reused under screening | SOURCE-VERIFIED + NUMERICALLY VERIFIED |
| unscreened behaviour unchanged | NUMERICALLY VERIFIED (see Regression) |

## Known limitations

1. **`ALField` and `EOLField` cannot be screened** and are refused. They optimise their orbitals through the
   B-spline tensor line (`buildScreenedPotentialPair` / `buildScreenedPotentialCache` -> `XL_CoulombTensor`),
   a separate construction the switch does not reach. Screening them means writing a Yukawa counterpart of
   that line; that is the natural next piece of work and was deliberately not attempted here.
2. **Breit and Gaunt are not screened** and are refused in combination with screening. The transverse
   photon-exchange kernels are a different operator; a Debye-screened Coulomb term beside a bare Breit term is
   not a defined approximation. (The task also placed Breit screening out of scope.)
3. **The average-atom line (`AaDFSField`, `AaHSField`) is not screened.** It is a separate plasma model with
   its own chemical potential, and out of the scope this task set.
4. **`Basics.AtomicFeatures`** computes `F^k`/`G^k` with the naive `RadialIntegrals.SlaterRk`, which has the
   Coulomb kernel written inline and is NOT reached by the switch. It builds feature vectors for the
   `DeepLearning` module and takes no part in any Hamiltonian, so this is a documentation point rather than an
   inconsistency -- but a feature vector computed during a screened session is a Coulomb one.
5. **`Plasma.LineShiftScheme` and `AutoIonization.computeLinesPlasma` NOW RETURN DIFFERENT NUMBERS.** They
   reach `SlaterRkDebyeHueckel`, whose body was replaced. The old body summed the correct expansion but
   TRUNCATED the `i_k` series at p <= 2: measured against the exact kernel it is off by 1.7e-4 at
   mu r_< = 1, 7.4e-3 at 2, 0.15 at 4 and **0.908 at mu r_< = 10** -- i.e. it lost 90 % of the interaction in
   exactly the strong-screening regime it exists for. The new numbers are the correct ones. No approved
   reference is affected (`test/approved/test-PlasmaShift-approved.sum` is referenced by nothing), and those
   two schemes remain what they were in every other respect -- in particular they still screen only the CI
   matrix and not the SCF, which is their own documented design and was not changed here.
6. **The screening parameter is a session global.** Two concurrent calculations in one Julia session cannot use
   different Debye lengths, and a threaded CI matrix reads one value. This follows from the design decision
   above and is the price of guaranteeing that the SCF and the Hamiltonian agree.
7. **`mu` is constant in r.** This is the Debye-Hueckel model as specified. Ion-sphere, Stewart-Pyatt,
   Ecker-Kroell, dynamic and Thomas-Fermi screening were all explicitly out of scope and none was attempted.
8. **Nothing here screens the electron-NUCLEUS interaction**, which remains the separate business of
   `Nuclear.nuclearPotentialDH`. A physically complete Debye plasma model would screen both; combining them is
   a decision for the maintainer, not a defect in this work.

### Future work (recorded, deliberately NOT done)

* a Yukawa counterpart of the B-spline tensor line, which would lift limitation 1;
* screening the Breit kernels (limitation 2);
* a `Plasma.Computation` scheme that screens the SCF as well as the CI matrix, now that the machinery exists;
* deleting the orphaned `test/approved/test-PlasmaShift-approved.sum` -- an editorial act for the maintainer.


---

# 05-Sep-2026 — FINAL ADVERSARIAL REVIEW AND RELEASE-READINESS PASS

Chronological, as everything above. The cleaned, reorganised version of all of this is
`docs/debye-electron-electron-screening.md`; what follows is the diary of the review itself.

## What was reviewed

`git status`, `git diff --stat`, `git diff`, `git diff --check`, then every modified file in full, plus
`progress.md`, `CLAUDE.md`, `CLAUDE.local.md` and the three earlier `tools/diag-debye-screening-*.jl`. The
eighteen review categories requested were worked through; the classification of every proposed change is below.

## Findings, classified

| # | Finding | Class | Action |
|---|---|---|---|
| 1 | `GBL_YUKAWA_X_FLOOR = 1e-25` is rank-blind. `kScaledYukawa` overflows from k = 11; `yukawaKernel` returned **Inf at k = 12 and NaN at k >= 14** where the Coulomb value is 1.0. The stability sweep had stopped at k = 8 and never saw it. | **REQUIRED FOR RELEASE SAFETY** | **FIXED** — rank-aware `RadialIntegrals.yukawaXFloor(k)`; identical (1e-25) for every k <= 10, so nothing reachable changed. A matching guard added to `buildScreenedPotentialYukawa`, which calls `kScaledYukawa` directly and so does not inherit the kernel's short-circuit; it RAISES rather than returning Inf. |
| 2 | A screened run left **no trace in its own output**. The switch is a session global and defaults to off, so a screened `.sum` or console log was indistinguishable from an unscreened one. | **REQUIRED FOR RELEASE SAFETY** | **FIXED** — `checkEeScreeningIsConsistent` gained a `printout` keyword and prints one line naming lambda_D and mu and stating what is and is not affected; both `performSCF` sites pass their own `printout` through. |
| 3 | Test gap: only `DFSField` was ever exercised screened. Four of the five modified `computePotential` methods were covered by inspection alone — including `HSField`, which reaches the screening by a **different** route (`Yk_ab`). | **REQUIRED FOR RELEASE SAFETY** | **FIXED** — suite checks (10) and (11) added: all five mean fields must respond, and in the right direction (Z(r) larger everywhere, since removing e-e repulsion makes the electronic part smaller). |
| 4 | Test gap: rank coverage stopped at k = 4 in the suite and k = 8 in the sweep — which is exactly why finding 1 survived. | **REQUIRED FOR RELEASE SAFETY** | **FIXED** — suite check (4b) covers k up to 20 at vanishing mu; the sweep now runs k = 0..16. |
| 5 | `git diff --check` reports trailing whitespace on 2 docstring lines. | **DO NOT CHANGE** | House style — 107 such lines already in `module-Defaults.jl`, 113 in `module-RadialIntegrals.jl`, and the two trailing spaces are a Markdown hard line break. Removing them would change rendering and make the new code inconsistent with its own file. |
| 6 | `rtol` is accepted but unused by `buildScreenedPotentialYukawa`. | **DO NOT CHANGE** | It is unused by the **pre-existing** Coulomb `buildScreenedPotential` too; the Yukawa version mirrors it for signature symmetry. Removing it is unrelated cleanup. |
| 7 | `Defaults.eeScreeningMu` was added to the module's `export` list. | **OPTIONAL CLEANUP** — not done | A legitimate query, documented, no name collision. Un-exporting functioning code for tidiness is outside the remit. |
| 8 | Accidental unrelated changes / dead code / duplicated logic / leftover debug output. | **none found** | Every diff hunk sits in an expected function (`git diff -U0` hunk-context listing); the only two `println`s added are the intentional user notice and a test summary line. |

## Diagnostics that were themselves wrong — corrected, not tolerated

Three test *criteria* condemned correct code during this pass. Each was repaired by making the comparison
honest, never by widening it. All three are written up in full in §12 of the report.

* **The exponential tail** was compared against `-mu` when the tail is `exp(-mu r)/(mu r)`; the 1/r prefactor
  contributes, and the correct chord prediction is `-mu - ln(r3/r1)/(r3-r1) = -0.879`. Measured: -0.879.
* **`Inf`/`NaN` at k = 10, x ~ 1e-30** turned out to be the *unscreened* Coulomb kernel's own Float64 limit:
  `(1e-30)^11` underflows to 0, so `r_<^k/r_>^(k+1)` is `Inf` screened or not, and JAC's pre-existing code
  behaves identically. Verified by evaluating both. Now excluded and **counted** in the sweep's output.
* **23 "screened exceeds Coulomb" violations** turned out to be Float64 **subnormals**: `ihat_10(5.6e-32)` is
  2.5e-323, which carries about two significant bits. Verified against `floatmin`. Now excluded and counted.

After all three corrections the sweep reports **13 544 combinations, 0 violations**, with 1 098 + 148
exclusions each stated with its reason in the diagnostic's own output.

## Tests re-run after the changes

| What | When | Result |
|---|---|---|
| `testModule_PlasmaScreening` standalone | after findings 1–4 | **[OK]**, 34.3 s |
| Full suite, `cd test && julia --project=.. runtests.jl` | after all source changes | **83 / 83, 7m42.1 s**, zero Fail/Error |
| The same, re-run after the documentation was written | final confirmation | **83 / 83, 7m42.4 s**, zero Fail/Error |
| `tools/diag-debye-screening-stability.jl` | after the diagnostic corrections | 0 violations over 13 544 combinations; worst accuracy 7.58e-15 vs 256-bit |
| Guard check (5 refusals + acceptance with screening off) | after finding 2 | all five refuse; `ALField` still accepted unscreened |

**The suite's evidence window was checked**, per the CLAUDE.md discipline: the latest `src/` or `test/`
modification time was 02:27:04 and the suite started at 02:35:29, so the tree stood still for the whole run and
the 83/83 applies to the final source state.

## Final numerical results (unchanged by the review except where noted)

* Regression: 82/82 baseline → 82/82 implemented → 83/83 with the test → **83/83 after the review**. No
  approved reference re-approved at any stage.
* Kernel accuracy: worst **7.58e-15** against a 256-bit reference over 210 combinations, k <= 16,
  mu*r = 1e-12..1e5. (Improved coverage; the previous figure, 2.90e-16, was over k <= 8.)
* Stability: **0 NaN, 0 Inf, 0 negative, 0 exceeding Coulomb** over 13 544 combinations.
* Rank robustness: finite and equal to the Coulomb kernel for **k = 0..30** at vanishing mu (was Inf/NaN
  from k = 12 before finding 1).
* Production: He/Be/Ne through `perform` at mu = 0.2 and 1.0; Coulomb limit `DeltaE/(-6mu)` → 1.000001 with a
  residual falling exactly as mu^2; weak-screening ratios 0.999982 / 0.999987 / 0.999994.
* SCF consistency: full vs frozen-orbital differs by 0.025 % at mu = 0.01, **9.05 % at mu = 1**.
* Performance: screened sweep 1.08–1.13x the Coulomb sweep, **flat** over a 5.4x range in N; end-to-end SCF
  ~1.7x (single runs, ±10 %).

## Known limitations at close

Unchanged from the list above, plus one added by finding 1: above rank 16 the small-mu Coulomb substitution
becomes approximate rather than exact (a graceful degradation replacing a NaN, in a regime no atomic Slater
integral reaches). The full, classified list — implementation / physical-model / numerical / untested /
unsupported — is §17 of `docs/debye-electron-electron-screening.md`.

## Repository state at close

Seven modified files, plus `docs/debye-electron-electron-screening.md`, `progress.md` and five
`tools/diag-debye-screening-*.jl`. **Nothing committed**, per CLAUDE.md. `CLAUDE.local.md`,
`OVERNIGHT_TASK.md`, `progress_before_resume.md`, `launch_claude.py` and `automation_logs/` are pre-existing
and are NOT part of this work.
