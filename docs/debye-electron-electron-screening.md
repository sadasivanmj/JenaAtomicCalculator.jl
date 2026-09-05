# Debye–Hückel Screening of the Electron–Electron Interaction in JAC

**A forensic implementation report**

| | |
|---|---|
| Repository | JenaAtomicCalculator.jl |
| Base commit | `75f963acac367e7365eed27854199a572a465d2c` (branch `master`) |
| Julia | 1.12.7 |
| Platform | Linux 6.8.0-138-generic, x86_64 |
| Report date | 05-Sep-2026 |
| Chronological development log | `progress.md` (this report is the cleaned permanent record; `progress.md` remains the diary) |

## How to read the evidence labels

Every substantive claim below carries one of the following. They are not decoration: a claim without a label is
a statement of intent or of context, and a claim with `[INFERENCE]` has **not** been executed or proved.

| Label | Meaning |
|---|---|
| `[PHYSICS]` | Standard physics or mathematics, stated as background. Not specific to this repository. |
| `[CODE EVIDENCE]` | Established by reading the source. The file and function are named so the reader can check. |
| `[ANALYTICAL VERIFICATION]` | Derived on paper, with the derivation shown or reproducible from what is shown. |
| `[NUMERICAL VERIFICATION]` | Measured by running code. The number quoted is an actual output, and the command that produced it is given in §18. |
| `[INFERENCE]` | A reasoned conclusion that was not directly executed. Treated as provisional. |
| `[NOT VERIFIED]` | Explicitly unproven. Listed so that it cannot be mistaken for the others. |

---

# 0. Executive Summary

## What Debye electron–electron screening means

`[PHYSICS]` An atom embedded in a plasma is not in vacuum. The free charges of the plasma rearrange around
every charge in the system and reduce its influence at long range. In the Debye–Hückel model this replaces the
bare Coulomb repulsion between two bound electrons,

$$V_{ee}(r_{12}) = \frac{1}{r_{12}},$$

by a Yukawa (screened) repulsion,

$$V_{ee}(r_{12}) = \frac{e^{-r_{12}/\lambda_D}}{r_{12}} = \frac{e^{-\mu r_{12}}}{r_{12}},
\qquad \mu \equiv \frac{1}{\lambda_D},$$

with $\lambda_D$ the Debye length. The interaction is unchanged at $r_{12} \ll \lambda_D$ and is cut off
exponentially beyond it.

## What JAC did before this work

`[CODE EVIDENCE]` JAC already carried a partial Debye–Hückel vocabulary, and it is important not to mistake it
for a finished capability:

* `Basics.DebyeHueckelModel(debyeLength)` — the parameter object existed and is reused unchanged.
* `Nuclear.nuclearPotentialDH` — Debye screening of the **nucleus** (a one-body effect). Complete, and
  deliberately untouched by this work.
* `InteractionStrength.XL_Coulomb_DH` → `RadialIntegrals.SlaterRkDebyeHueckel` — a screened *two-electron*
  Slater integral, reachable only from `Plasma.LineShiftScheme` and `AutoIonization.computeLinesPlasma`.

## What was missing

Four distinct deficiencies, of three different kinds:

1. **Missing functionality.** Nothing screened the e–e interaction on the ordinary calculation path. A
   `perform(Atomic.Computation)` had no way to request it; the switch simply did not exist.
2. **Architecture limitation.** `Plasma.LineShiftScheme` screened only the CI matrix, on orbitals converged in
   the *bare* Coulomb field. The self-consistent field was never screened, so the calculation was internally
   inconsistent by construction. `[CODE EVIDENCE]` — the pre-existing approved output
   `test/approved/test-PlasmaShift-approved.sum` says so in its own header: *"Plasma screening included
   perturbatively in CI matrix but not in SCF field."*
3. **Implementation limitation.** `SlaterRkDebyeHueckel` summed the correct multipole expansion but
   **truncated the $i_k$ power series after three terms** ($p \le 2$), which is accurate only while
   $\mu r_< \ll 1$ — i.e. it failed in exactly the strong-screening regime the function exists for.
4. **Numerical limitation.** The same routine performed an explicit $O(N_r^2)$ double integration over the
   radial grid, while JAC's unscreened Coulomb machinery is $O(N_r)$.

## What was changed

One substitution, applied at two places, controlled by one authoritative parameter:

$$\frac{r_<^k}{r_>^{k+1}} \;\longrightarrow\; \mu\,(2k+1)\, i_k(\mu r_<)\, k_k(\mu r_>).$$

| # | File | What |
|---|---|---|
| 1 | `src/module-Defaults.jl` | the screening state, its setter, validator and accessor |
| 2 | `src/module-RadialIntegrals.jl` | the Bessel primitives, the kernel, the screened $O(N)$ potential sweep, and the **single branch point** that every Slater integral passes through |
| 3 | `src/module-BasicsAZ-inc-compute.jl` | the five e–e mean-field potentials that drive the SCF |
| 4 | `src/module-InteractionStrength.jl` | the two two-electron caches, whose keys now carry $\mu$ |
| 5 | `src/module-SelfConsistent.jl` | the consistency guard, and the notice that makes a screened run self-identifying |
| 6 | `src/module-TestFrames-inc-plasma.jl`, `test/runtests.jl` | the regression test and its entry in the previously empty plasma testset |

## What was validated

`[NUMERICAL VERIFICATION]` A seven-layer hierarchy, from the Bessel functions to a full atomic energy, plus
two mathematically independent references (a 256-bit arbitrary-precision implementation, and an angular
Legendre projection that uses no Bessel function at all). The production path was proved to be the tested path
by a caller census **and** a sampling profile of a real running calculation. Details in §13 and §14.

## Final test status

`[NUMERICAL VERIFICATION]` See §14.1 for the full table and §20 for the command.

```
Baseline (before any edit)      82 / 82   9m00.2 s
After the implementation        82 / 82   7m42.6 s
After the regression test       83 / 83   7m42.8 s
After the final review          83 / 83   7m42.1 s
```

No `test/approved/*.sum` reference was modified or re-approved at any stage.

## Major limitations

Stated in full in §17. In brief: `ALField` and `EOLField` cannot be screened and are **refused**, not
half-screened; Breit and Gaunt are likewise refused in combination with screening; the average-atom fields and
the machine-learning feature extractor are unscreened by design; and the parameter is a session global.

```
Status: COMPLETE
```

This label is used because both conditions set for it are met and were measured, not asserted: the production
Debye e–e path has been exercised by a real JAC atomic calculation (§14.5, §3.3), and the unscreened
regression is clean (§14.1). It is **not** a claim that every JAC calculation mode supports screening — §17
lists the modes that do not, and those refuse rather than answer.

---

# 1. Motivation

## 1.1 Why this is needed

`[PHYSICS]` Atomic structure calculations for plasma spectroscopy — for laser-produced plasmas, tokamak edge
and divertor regions, stellar interiors, inertial-confinement targets — cannot use vacuum atomic data without
qualification. The plasma environment shifts level energies, changes transition rates, moves ionization
thresholds and can push weakly bound states into the continuum. The Debye–Hückel model is the simplest
description of that environment that is still derived rather than fitted, and it is the standard first
approximation in the field.

## 1.2 The physical origin of the screening

`[PHYSICS]` In a plasma of temperature $T$ and free-electron density $n_e$, a test charge polarises the
surrounding mobile charges. Linearising the Poisson–Boltzmann equation for the induced density gives a
potential falling as $e^{-r/\lambda_D}/r$ rather than $1/r$, with

$$\lambda_D = \sqrt{\frac{k_B T}{4\pi n_e e^2}} \quad \text{(Gaussian units)}.$$

The Debye length is the distance over which the plasma neutralises a disturbance. Two important consequences
follow, and both matter for the code: the model is a **linear-response** result, so it is valid for weak
coupling; and $\lambda_D$ is a property of the *environment*, entering the atomic problem as a single external
parameter, not as something the calculation determines.

## 1.3 One-body nuclear screening is a different thing from two-body e–e screening

This distinction is the single most important one in this report, because the two are easy to conflate and JAC
contains both.

| | nuclear (one-body) | electron–electron (two-body) |
|---|---|---|
| what is screened | the electron–nucleus attraction | the electron–electron repulsion |
| operator | $-\dfrac{Z e^{-\mu r}}{r}$, a function of **one** coordinate | $\dfrac{e^{-\mu r_{12}}}{r_{12}}$, a function of the **separation** of two |
| where it acts | the one-particle Dirac Hamiltonian, i.e. a local potential | the two-particle interaction, i.e. the Slater integrals and the mean field |
| in JAC | `Nuclear.nuclearPotentialDH` — pre-existing, **untouched** | this work |
| sign of the energy shift | binds *less*, raises the energy | repels *less*, **lowers** the energy |

`[CODE EVIDENCE]` The two are implemented in entirely separate places in JAC and share no code. This report
concerns only the second. `[PHYSICS]` A physically complete Debye plasma model would apply both; combining
them is a modelling decision, not a defect in either.

## 1.4 Why $r_{12}$ cannot simply be treated as a radial grid coordinate

This is the heart of why the change belongs where it does, and it is where a naive implementation goes wrong.

`[PHYSICS]` A relativistic atomic-structure code does not represent the two-electron wavefunction on a
two-dimensional $(r_1, r_2)$ grid, let alone on an $r_{12}$ grid. It represents each orbital as a radial
function on a one-dimensional grid and handles the angular structure algebraically. The quantity $r_{12}$ —
the distance between two electrons — is **not a coordinate of that representation at all**. It is

$$r_{12} = \sqrt{r_1^2 + r_2^2 - 2 r_1 r_2 \cos\omega},$$

where $\omega$ is the angle between the two position vectors, and it therefore depends on the *angular*
variables that the algebra has already integrated away.

The only way a two-electron interaction enters such a code is through its **multipole expansion**: an
expansion in $P_k(\cos\omega)$ whose coefficients are functions of $r_1$ and $r_2$ separately. The angular
factor is absorbed into angular-momentum algebra (Wigner coefficients), and what is left for the radial code is
one function of two radii per multipole rank $k$. Screening the interaction therefore means **replacing that
radial function**, and nothing else. Any implementation that tries to attach $e^{-\mu r}$ to a single orbital,
or to a single radial coordinate, is screening a different operator — see §5.2, where this is worked out
explicitly.

`[INFERENCE]` This is why the modification belongs at the radial-kernel level and why the angular machinery
needs no change at all: the substitution replaces the coefficient of $P_k(\cos\omega)$ and leaves the
$P_k(\cos\omega)$ itself, hence every Wigner coefficient, exactly as it was. §9.5 turns this inference into
`[CODE EVIDENCE]` by showing where in the source the two factors are formed.

---

# 2. Theoretical Background

## 2.1 The Coulomb electron–electron interaction and its multipole expansion

`[PHYSICS]` The generating function of the Legendre polynomials,

$$\frac{1}{\sqrt{1 - 2 t \cos\omega + t^2}} = \sum_{k=0}^{\infty} t^k P_k(\cos\omega), \qquad |t| < 1,$$

applied with $t = r_</r_>$ where

$$r_< = \min(r_1, r_2), \qquad r_> = \max(r_1, r_2),$$

gives the Laplace expansion

$$\boxed{\;\frac{1}{r_{12}} = \sum_{k=0}^{\infty} \frac{r_<^k}{r_>^{k+1}} P_k(\cos\omega).\;}$$

The choice $t = r_</r_>$ rather than $r_1/r_2$ is what makes the series converge everywhere, and it is the
origin of the $r_<$ / $r_>$ structure that pervades every atomic-structure code. The radial object the code
must supply is therefore the **rank-$k$ Coulomb kernel**

$$U_k^{\text{Coul}}(r_1, r_2) = \frac{r_<^k}{r_>^{k+1}},$$

a function that is continuous but has a *kink* — a discontinuous derivative — on the line $r_1 = r_2$. That
kink is not incidental; JAC's production quadrature is built around it (`SlaterRkKinkAware`), and the screened
kernel inherits it unchanged.

## 2.2 The Debye/Yukawa interaction

`[PHYSICS]` The screened interaction is

$$V_{ee}(r_{12}) = \frac{e^{-\mu r_{12}}}{r_{12}}, \qquad \mu = \frac{1}{\lambda_D},$$

which is the Green's function of the modified Helmholtz operator $(\nabla^2 - \mu^2)$, exactly as $1/r_{12}$
is the Green's function of $\nabla^2$. This is why the expansion below exists and has the form it does: the
$r^k$ and $r^{-(k+1)}$ solutions of the radial Laplace equation are replaced by the two solutions of the
radial modified Helmholtz equation, which are precisely the modified spherical Bessel functions.

## 2.3 The Yukawa multipole (Gegenbauer) expansion

`[PHYSICS]` The corresponding expansion is

$$\boxed{\;\frac{e^{-\mu r_{12}}}{r_{12}} = \mu \sum_{k=0}^{\infty} (2k+1)\, i_k(\mu r_<)\, k_k(\mu r_>)\,
P_k(\cos\omega),\;}$$

so the entire implementation is the single substitution, at fixed rank $k$,

$$\frac{r_<^k}{r_>^{k+1}} \;\longrightarrow\; U_k^{\text{Yuk}}(r_1,r_2) = \mu\,(2k+1)\, i_k(\mu r_<)\,
k_k(\mu r_>).$$

### The normalisation, and why it is the most dangerous line in this report

`[PHYSICS]` The modified spherical Bessel functions used here are fixed by

$$i_0(x) = \frac{\sinh x}{x}, \qquad k_0(x) = \frac{e^{-x}}{x}.$$

`[CODE EVIDENCE]` **A widely used competing convention differs by a factor $\pi/2$ in the second kind.** Both
GSL's `gsl_sf_bessel_kl_scaled` and the combination `sqrt(pi/2x) * besselkx(k+1/2, x)` available through
Julia's `SpecialFunctions` return

$$k_k^{\text{(library)}}(x) = \frac{\pi}{2}\, k_k(x),$$

i.e. $\pi/2$ times the function required above. A code that adopts the library convention without the
conversion is wrong by a constant factor $2/\pi \approx 0.6366$ in **every** screened two-electron matrix
element. Such an error is unusually hard to catch: it is a smooth constant, so the calculation still converges,
the energies still look plausible, the Coulomb limit is still approached in the right *direction*, and only an
absolute check against a closed form reveals it.

This is treated as a first-class validation problem in this work, and the defence is deliberately redundant:
the convention is stated in the module banner
(`src/module-RadialIntegrals.jl`), restated in the docstring of `kScaledYukawa`, and — most importantly —
pinned by a test against a closed form that contains **no Bessel function at all** (§14.2, check 2). The same
convention issue is documented in the literature of other atomic-structure codes; FAC (Gu, *Can. J. Phys.*
**86**, 675 (2008)) uses the same $k_0 = e^{-x}/x$ normalisation for the identical expansion.

## 2.4 The Coulomb limit

`[ANALYTICAL VERIFICATION]` The small-argument behaviour of the two functions is

$$i_k(x) \xrightarrow{x \to 0} \frac{x^k}{(2k+1)!!}, \qquad
k_k(x) \xrightarrow{x \to 0} \frac{(2k-1)!!}{x^{k+1}}.$$

Substituting,

$$\mu (2k+1)\, i_k(\mu r_<)\, k_k(\mu r_>)
\;\longrightarrow\;
\mu (2k+1)\, \frac{(\mu r_<)^k}{(2k+1)!!} \cdot \frac{(2k-1)!!}{(\mu r_>)^{k+1}}
= (2k+1)\frac{(2k-1)!!}{(2k+1)!!}\cdot \frac{r_<^k}{r_>^{k+1}}.$$

The powers of $\mu$ cancel exactly: $\mu \cdot \mu^k / \mu^{k+1} = 1$. Since $(2k+1)!! = (2k+1)(2k-1)!!$, the
prefactor is $(2k+1)\cdot\frac{1}{2k+1} = 1$, and

$$\boxed{\;\mu(2k+1) i_k(\mu r_<) k_k(\mu r_>) \;\xrightarrow{\mu \to 0}\; \frac{r_<^k}{r_>^{k+1}}\;}$$

**exactly**, with no residual constant. `[NUMERICAL VERIFICATION]` §14.3.

### The order of the approach, which is not what one first expects

`[ANALYTICAL VERIFICATION]` It is natural to assume the correction is $O(\mu^2)$ at every rank. **It is not.**
For $k = 0$ the kernel is available in closed form,

$$\mu\, i_0(\mu a)\, k_0(\mu b) = \mu \cdot \frac{\sinh(\mu a)}{\mu a}\cdot\frac{e^{-\mu b}}{\mu b}
= \frac{\sinh(\mu a)\, e^{-\mu b}}{\mu\, a\, b},$$

and expanding for small $\mu$ with $\sinh(\mu a) = \mu a + O(\mu^3)$ and $e^{-\mu b} = 1 - \mu b + O(\mu^2)$,

$$\frac{\mu a (1 - \mu b)}{\mu a b} + O(\mu^2) = \frac{1}{b} - \mu + O(\mu^2).$$

**The $k=0$ correction is the constant $-\mu$**, so its *relative* error is $O(\mu r_>)$ — first order. Only
for $k \ge 1$ does the linear term cancel and the approach become second order. This is not a curiosity: that
constant $-\mu$ per electron pair *is* the physics of §2.5, and a test written to expect order 2 at $k=0$
would fail a correct implementation. `[NUMERICAL VERIFICATION]` The observed orders are $0.996$–$0.999$ for
$k=0$ and $1.995$–$2.000$ for $k \ge 1$ (§14.3), and the suite test asserts the *observed order*, not
agreement at a fixed $\mu$.

## 2.5 The weak-screening limit, and why it is the best single test in this report

`[ANALYTICAL VERIFICATION]` Expanding the interaction itself rather than its multipoles,

$$\frac{e^{-\mu r_{12}}}{r_{12}} = \frac{1}{r_{12}}\left(1 - \mu r_{12} + \frac{\mu^2 r_{12}^2}{2} - \cdots\right)
= \frac{1}{r_{12}} - \mu + \frac{\mu^2 r_{12}}{2} - \cdots$$

The first-order term is $-\mu$, **a constant** — it does not depend on $r_{12}$, and therefore not on the
wavefunction. Summing over all distinct pairs of $N$ electrons,

$$\boxed{\;\Delta E_{ee} = \Big\langle \sum_{i<j}\big(V^{\text{scr}}_{ij} - V^{\text{Coul}}_{ij}\big)\Big\rangle
\;\xrightarrow{\mu \to 0}\; -\mu\,\frac{N(N-1)}{2}.\;}$$

Three properties make this the most valuable validation target available:

1. **It is parameter-free.** No orbital, basis, grid or convergence setting enters. The right-hand side is
   arithmetic.
2. **It is independent of the wavefunction**, because a constant operator has the same expectation value in
   every normalised state. Orbital relaxation is therefore *second* order and cannot contaminate the leading
   term — which is also why it is a valid test of a self-consistent calculation, not only a frozen one.
3. **It fixes the absolute normalisation and the sign at once.** A kernel wrong by $2/\pi$ (§2.3) fails it; a
   kernel with the wrong sign fails it; a kernel that screens only the direct term fails it, because the
   counting $N(N-1)/2$ includes exchange pairs.

`[NUMERICAL VERIFICATION]` §14.4 and §14.5.

---

# 3. JAC Architecture Before Modification

Everything in this section was established by reading `src/`. Nothing is generic atomic-physics description.

## 3.1 The objects

`[CODE EVIDENCE]`

| Physics object | JAC representation |
|---|---|
| relativistic orbital | `Radial.Orbital`, with tabulated `P`, `Q` and a `subshell::Subshell` |
| radial grid | `Radial.Grid` — `grid.r`, quadrature weights `grid.wr`, `grid.NoPoints` |
| radial integration | `sum_i f[i]*grid.wr[i]`; `RadialIntegrals.V0`; per-cell 8-point Gauss–Legendre `RadialIntegrals.cellIntegral` |
| radial pair density | $\rho_{ac}(r) = P_aP_c + Q_aQ_c$, formed inline at each use |
| total electron density | $\rho_t(r) = \sum_a \mathrm{occ}_a (P_a^2 + Q_a^2)$, formed inline in `Basics.computePotential` |
| screened potential $Y_k$ | `RadialIntegrals.buildScreenedPotential(k, b, d, grid)` — an $O(N)$ forward/backward sweep |
| Slater integral $R^k(abcd)$ | `RadialIntegrals.SlaterRkKinkAware` (production), `SlaterRk` (naive $O(N^2)$), `SlaterRkReference` (oracle) |
| angular coefficients | `SpinAngular.computeCoefficients(TwoParticleOperator(0,plus), ...)` → `coeff.V`, `coeff.nu`, `coeff.a..d` |
| two-electron matrix element | `InteractionStrength.XL_Coulomb(L,a,b,c,d,grid[,cache])` |
| CI Hamiltonian | `Hamiltonian.setupMatrix` → `Hamiltonian.performCI` |
| SCF driver | `SelfConsistent.performSCF` → `Basics.scfProcedure(scField)` |
| mean-field potential | `Basics.computePotential(scField, grid, level)` |

**A naming warning that matters.** `RadialIntegrals.buildScreenedPotential` has nothing to do with plasma
screening. In this code "screened potential" is Hartree's $Y_k$ (GRASP's `YZK`) — the potential produced by one
orbital-pair density, "screened" in the sense that the charge inside radius $r$ shields the charge outside.
It is the *Coulomb* object. This report writes $V_k$ for it and reserves "screening" for the plasma effect.

## 3.2 The two Coulomb paths, traced

`[CODE EVIDENCE]` This is the traced call graph, not a schematic. The key discovery is that there are **two
independent constructions of the same physics**, and only one of them goes through the Slater machinery.

```
                          AsfSettings (eeInteraction, eeInteractionCI, scField)
                                              |
                          SelfConsistent.performSCF(configs, nm, grid, settings)
                                              |
              +-------------------------------+--------------------------------+
              |                                                                |
        PATH C: the SCF                                              PATH A: the Hamiltonian
              |                                                                |
   scfProcedure(DFSField) == :meanFieldIteration                    Hamiltonian.performCI
              |                                                                |
   SelfConsistent.solveMeanFieldBasis[Anderson]                     Hamiltonian.setupMatrix
              |                                                                |
   Basics.computePotential(DFSField, grid, level)                   SpinAngular.computeCoefficients
              |                                                       (angular; UNCHANGED)
   rho_t(r) from the orbitals;                                                  |
   V(r) = int rho_t(r')/r_>  --  its OWN i,j double loop            InteractionStrength.XL_Coulomb
              |                                                        (DIRECT *and* EXCHANGE)
   Basics.add(nuclearPotential, wp)                                             |
              |                                                     RadialIntegrals.SlaterRkKinkAware
   Bsplines.setupLocalMatrix -> NEW ORBITALS                                    |
              |                                                     RadialIntegrals.buildScreenedPotential
              +---------------- feeds ------------------------------->    (the O(N) V_k sweep)
                                                                                |
                                                                        contract with rho_ac
                                                                                |
                                                                        energies
```

`[CODE EVIDENCE]` The mean-field loop, in isolation:

```
SCF  ->  orbitals  ->  rho_t(r)  ->  Basics.computePotential  ->  local potential
                            ^                                            |
                            +--- Bsplines.setupLocalMatrix <-------------+
                                   (diagonalise -> new orbitals; iterate)
```

**The architectural finding.** JAC is neither purely "one shared kernel" nor purely "separate direct and
exchange". It is **Path A *and* Path C simultaneously**:

* **Path A** — one generic radial kernel serves the direct **and** the exchange Slater integrals. They differ
  only in which orbitals the angular coefficient places in the four slots. One insertion point covers both.
* **Path C** — a **separate, independent** $k=0$ direct-potential calculation drives the SCF. It never calls
  the Slater machinery; it re-derives $\int \rho(r')/r_>$ with its own `max(r_i, r_j)` double loop. Screening
  Path A alone would leave the *orbitals* unscreened while the Hamiltonian built from them was screened.
* A **third** line, the B-spline tensor construction
  (`buildScreenedPotentialPair` / `buildScreenedPotentialCache` → `XL_CoulombTensor`), is used **only** by
  `scField = ALField()` and `EOLField()`. It is not reached by this work — see §17.

`[NUMERICAL VERIFICATION]` The Path C omission is not hypothetical. Measured on Be $1s^2 2s^2$, screening the
Slater integrals alone (frozen Coulomb orbitals) is wrong by 0.025 % of the shift at $\mu = 0.01$, rising to
**9.05 %** at $\mu = 1\,a_0^{-1}$ — §14.7.

## 3.3 The user-facing entry point

`[CODE EVIDENCE]` `src/module-BasicsAZ-inc-perform.jl:15` defines
`Basics.perform(computation::Atomic.Computation; output)`, and line 21 of the same file calls
`SelfConsistent.performSCF(computation.configs, nModel, computation.grid, computation.asfSettings)`. The
ordinary user route and the route the tests exercise are therefore the same route, and both guards installed by
this work sit at `performSCF`.

---

# 4. Physics-to-Code Mapping

`[CODE EVIDENCE]` Every row corresponds to real source. Line numbers are as of the final state of this work.

| Physics quantity | Mathematical expression | JAC object | File | Function / type | Role |
|---|---|---|---|---|---|
| e–e Coulomb interaction | $1/r_{12}$ | implicit; never formed | — | — | expanded in multipoles before it reaches the radial code |
| Coulomb multipole kernel | $r_<^k / r_>^{k+1}$ | written inline | `module-RadialIntegrals.jl` | `SlaterRk`, `SlaterRkComponent`, `Yk_ab`, `buildScreenedPotential` | the object replaced by this work |
| the same kernel, again | $1/r_>$ at $k=0$ | `max(r_i,r_j)` loop | `module-BasicsAZ-inc-compute.jl` | 5 × `Basics.computePotential` | the SCF's independent copy (Path C) |
| **Yukawa kernel** | $\mu(2k+1) i_k(\mu r_<) k_k(\mu r_>)$ | **new** | `module-RadialIntegrals.jl` | `yukawaKernel` | the replacement |
| scaled Bessel, 1st kind | $\hat\imath_k(x) = e^{-x} i_k(x)$ | **new** | `module-RadialIntegrals.jl` | `iScaledYukawa` | overflow-free $i_k$ |
| scaled Bessel, 2nd kind | $\hat k_k(x) = e^{x} k_k(x)$ | **new** | `module-RadialIntegrals.jl` | `kScaledYukawa` | underflow-free $k_k$ |
| Hartree $Y_k$ / potential | $V_k(r)$ | $O(N)$ sweep | `module-RadialIntegrals.jl` | `buildScreenedPotential` (Coulomb), `buildScreenedPotentialYukawa` (**new**) | the separable two-region integral |
| Slater integral | $R^k(abcd)$ | | `module-RadialIntegrals.jl` | `SlaterRkKinkAware`, `SlaterRkYukawa` (**new**) | contracts $V_k$ against $\rho_{ac}$ |
| **Direct** Slater integral | $F^k(ab) = R^k(abab)$ | slot pattern $(a,b,a,b)$ | `module-SpinAngular.jl` | `Coefficient2p` slots | **same** radial routine |
| **Exchange** integral | $G^k(ab) = R^k(abba)$ | slot pattern $(a,b,b,a)$ | `module-SpinAngular.jl` | `Coefficient2p` slots | **same** radial routine |
| two-electron ME | $X^L(abcd)$ | angular $C^L$ × radial | `module-InteractionStrength.jl` | `XL_Coulomb`, `XL_CoulombKinkAware` | where angular meets radial |
| Hartree/direct potential | $\int \rho_t(r')\,U_0(r,r')\,dr'$ | | `module-BasicsAZ-inc-compute.jl` | `computePotential(::DFSField/::HSField/::KSField/::CHField, …)` | drives the SCF |
| electron density | $\rho_t(r)=\sum_a \mathrm{occ}_a(P_a^2+Q_a^2)$ | `rhot` | `module-BasicsAZ-inc-compute.jl` | inline in the above | the source of the mean field |
| pair density | $\rho_{bd}=P_bP_d+Q_bQ_d$ | `rhoBd` | `module-RadialIntegrals.jl` | inline in the sweeps | the source of $V_k$ |
| **Screening parameter** | $\mu = 1/\lambda_D$ | `Float64` | `module-Defaults.jl` | `eeScreeningMu()` | the single authority |
| Screening state | the model | `Basics.AbstractPlasmaModel` | `module-Defaults.jl` | `GBL_EE_SCREENING` | default `NoPlasmaModel()` |
| Debye length | $\lambda_D$ [$a_0$] | `debyeLength` | `module-Basics-inc-abstract.jl` | `Basics.DebyeHueckelModel` | **pre-existing**, reused |
| ME cache | memoised $X^L$ | `Dict` | `module-InteractionStrength.jl` | `XLCache` | key now carries $\mu$ |
| $V_k$ cache | memoised potential | `Dict` | `module-RadialIntegrals.jl` | `vkCache` argument | per-call by contract |
| Consistency guard | — | — | `module-SelfConsistent.jl` | `checkEeScreeningIsConsistent` | refuses what cannot be screened |
| Nuclear screening | $-Ze^{-\mu r}/r$ | pre-existing | `module-Nuclear.jl` | `nuclearPotentialDH` | **not touched** (§16) |

---

# 5. Why the Existing Implementation Was Insufficient

## 5.1 The four deficiencies, classified

**(a) Missing functionality — no e–e screening on the ordinary path.**
`[CODE EVIDENCE]` Before this work, `RadialIntegrals.SlaterRkDebyeHueckel` was reachable from exactly two
places: `Basics.compute(JP, basis, nm, grid, settings, plasmaModel)` (the `Plasma.LineShiftScheme` CI step) and
`AutoIonization.computeLinesPlasma`. An ordinary `Atomic.Computation` had no route to it and no setting that
would select it. The capability existed as a component, not as a feature.

**(b) Architecture limitation — the SCF was never screened.**
`[CODE EVIDENCE]` `Plasma.perform(::LineShiftScheme, ...)` computes field-free orbitals and then re-does only
the CI matrix under the plasma model. The pre-existing approved output states this in its own header:

```
  Plasma shifts for Debye-Hueckel model:
     + lambda = 0.25
     + Plasma screening included perturbatively in CI matrix but not in SCF field.
```

`[NUMERICAL VERIFICATION]` The size of that omission was measured in this work: 9.05 % of the total shift at
$\mu = 1\,a_0^{-1}$ for Be (§14.7). It is a deliberate and documented approximation in that scheme, not a bug —
but it is not a self-consistent screened calculation, which is what this work provides.

**(c) Implementation limitation — a truncated series.**
`[CODE EVIDENCE]` The retired body of `SlaterRkDebyeHueckel` read (abbreviated):

```julia
function ul_DH(L::Int64, s::Float64, r::Float64)
    sum = 0.;   suma = 0.
    for  p = 0:2                       # <-- the power series of i_k, cut off after three terms
        for q = 0:L
            sum = sum + (2^(L-q)) * (lambda^(L+p+p-q)) * factorial(L+q) * factorial(L+p) /
                        ( factorial(L+L+p+p+1) * factorial(L-q) * factorial(p) * factorial(q)) *
                        (s^(L+p+p)) * exp(-lambda*r) / (r^(q+1))
        end
        if (p == 2) suma = sum   end
    end
    return( (L+L+1) * sum )
end
```

The expansion is the correct one — it is $\mu(2k+1) i_k k_k$ written out — but the $p$ loop is the ascending
power series of $i_k$, and stopping at $p = 2$ is accurate only while $\mu r_< \ll 1$.

`[NUMERICAL VERIFICATION]` Measured against the exact kernel, the retired form is in error by

| $\mu r_<$ | 0.02 | 1 | 2 | 4 | 10 |
|---|---|---|---|---|---|
| relative error, $k=0$ | 1.3e-14 | 1.7e-4 | 7.4e-3 | 0.15 | **0.908** |

i.e. it lost **91 % of the interaction** at $\mu r_< = 10$ — precisely the strong-screening regime the routine
exists to describe. (`tools/diag-debye-screening-performance.jl`, which reproduces the retired code verbatim
and compares it pointwise with `yukawaKernel`.)

**(d) Numerical limitation — $O(N_r^2)$ where JAC is $O(N_r)$.**
`[CODE EVIDENCE]` The retired routine ran an explicit double loop `for r = 2:mtp_ac, for s = 2:mtp_bd` over the
grid. `[NUMERICAL VERIFICATION]` It cost 13.983 ms against 1.062 ms for the screened $O(N)$ sweep at 1330 grid
points — a factor **13.2** — and the gap grows with $N$.

## 5.2 Why multiplying the Coulomb interaction by an exponential is not the screened interaction

This is the error the whole design is arranged to make impossible, and it is worth stating precisely because
the wrong expressions look superficially reasonable.

`[PHYSICS]` The required operator is $e^{-\mu r_{12}}/r_{12}$, where $r_{12} = |\mathbf r_1 - \mathbf r_2|$.
Consider the three plausible-looking alternatives:

**(i) $e^{-\mu r_1}/r_1$ or $e^{-\mu r_2}/r_2$ — a one-body operator.** These depend on a *single* electron's
distance from the *nucleus*. They are not a two-electron interaction at all; they are the nuclear screening of
§1.3 applied to the wrong term. Adding such a factor screens each electron's attraction to the origin, which is
`Nuclear.nuclearPotentialDH`'s business and is explicitly out of scope here. The physical content is entirely
different: one binds less, the other repels less, and the two shift the energy in *opposite* directions.

**(ii) $e^{-\mu r_>}\,r_<^k/r_>^{k+1}$ — the Coulomb kernel times an exponential of the larger radius.** This is
the most tempting error, because it has the right shape and the right limit at $\mu \to 0$. It is nevertheless
wrong: an exponential in $r_>$ is not an exponential in $r_{12}$. The correct rank-$k$ coefficient is the
projection

$$U_k(r_1,r_2) = \frac{2k+1}{2}\int_{-1}^{1} P_k(w)\,\frac{e^{-\mu r_{12}(w)}}{r_{12}(w)}\,dw,
\qquad r_{12}(w) = \sqrt{r_1^2 + r_2^2 - 2r_1r_2 w},$$

and performing that integral gives $\mu(2k+1) i_k(\mu r_<) k_k(\mu r_>)$ — a product of two *different*
functions, not a single exponential factor. `[NUMERICAL VERIFICATION]` This projection is computed
independently by quadrature in `tools/diag-debye-screening-kernel.jl`, using **no Bessel function of any
kind**, and agrees with the implemented kernel to ~1e-15 (§14.2, check 3). Form (ii) does not.

**(iii) $e^{-\mu(r_1+r_2)}/r_{12}$ or similar separable guesses.** Separable in $r_1$ and $r_2$, hence
expressible as a product of one-body factors, hence — by the same argument as (i) — not the two-body operator.

`[INFERENCE]` The general statement: $e^{-\mu r_{12}}$ is *not* separable in $r_1$ and $r_2$, and the only way
to represent a non-separable two-body operator in a one-dimensional radial code is through the multipole
expansion. This is why the modification must replace the radial kernel and can do nothing else.

---

# 6. Implementation Design

## 6.1 The screening state

`[CODE EVIDENCE]` `src/module-Defaults.jl` holds

```julia
GBL_EE_SCREENING             = Basics.NoPlasmaModel()
```

beside JAC's other `GBL_` settings, with a long `##` block recording the reasoning below. The only accepted
values are `Basics.NoPlasmaModel()` (the default) and `Basics.DebyeHueckelModel(lambda_D)`.

### Why a module global rather than a field of `AsfSettings`

This was the one genuinely contested design decision, and it was taken on the following evidence.

`[CODE EVIDENCE]` The two functions that must agree about $\mu$ are
`RadialIntegrals.buildScreenedPotential(k, b, d, grid)` and
`Basics.computePotential(scField, grid, level)`. **Neither receives a settings object, and neither is called
from anywhere that holds one** — the first takes a rank, two orbitals and a grid; the second takes a field
type, a grid and a level. Threading a parameter to both would touch eight modules and would change the
positional `AsfSettings` constructor used by roughly 200 call sites.

`[INFERENCE]` Against that cost stands a specific hazard: a screened run whose SCF and whose Hamiltonian
disagreed about $\mu$ would produce a plausible, converged, *wrong* number with nothing in the output to say
so. A single authoritative value removes that failure mode by construction. The pattern is not foreign to the
codebase — `GBL_FRAMEWORK` and `GBL_CONT_POTENTIAL` are handled the same way.

`[CODE EVIDENCE]` The same design choice was made independently in FAC, whose Debye e–e screening is likewise
carried by a file-scope flag and read through an accessor rather than passed down the call chain.

The cost is stated honestly in §17: two concurrent calculations in one Julia session cannot use different
Debye lengths.

## 6.2 The user-facing API

```julia
using JenaAtomicCalculator

# ON.  The argument is the DEBYE LENGTH lambda_D, in a_o (Bohr radii).
Defaults.setDefaults("e-e screening", Basics.DebyeHueckelModel(5.0))

wa = Atomic.Computation(Atomic.Computation(), name="Be", grid=Radial.Grid(true),
                        nuclearModel=Nuclear.Model(4.), configs=[Configuration("1s^2 2s^2")])
wb = perform(wa; output=true)

# OFF -- the default, and JAC's behaviour before this switch existed.
Defaults.setDefaults("e-e screening", Basics.NoPlasmaModel())

Defaults.getDefaults("e-e screening")   # the model currently in force
Defaults.eeScreeningMu()                # mu = 1/lambda_D [a_o^-1], EXACTLY 0.0 when unscreened
```

`[CODE EVIDENCE]` Setting the screening prints a notice naming both $\lambda_D$ and $\mu$ and stating what is
and is not affected. In addition, **every screened `performSCF` announces itself**:

```
> Debye-Hueckel e-e screening is IN FORCE for this calculation:  lambda_D = 5.0 a_o, mu = 0.2 a_o^-1.
>   1/r_12 -> exp(-mu r_12)/r_12 in the Slater integrals AND in the direct SCF potential; the
>   electron-nucleus potential is unchanged.
```

That second line was added during the final review (§7.6). Its purpose is narrow and important: because the
switch is a session global and defaults to off, without it a screened log or summary file would be
indistinguishable from an unscreened one, and a number recovered months later would carry no record of which it
was.

## 6.3 Variables and conventions

| Code variable | Physical meaning | Units | Where defined | Where used |
|---|---|---|---|---|
| `debyeLength` | $\lambda_D$ | $a_0$ | `Basics.DebyeHueckelModel` (pre-existing) | read once, by `eeScreeningMu` |
| `GBL_EE_SCREENING` | the screening model | — | `module-Defaults.jl` | read only via `eeScreeningMu()` |
| `mu` | $\mu = 1/\lambda_D$ | $a_0^{-1}$ | returned by `Defaults.eeScreeningMu()` | `yukawaKernel`, `buildScreenedPotentialYukawa`, `computePotential`, `Yk_ab`, cache keys |
| `x`, `xSmall`, `xLarge` | $\mu r$, $\mu r_<$, $\mu r_>$ | dimensionless | local | the Bessel arguments |
| `rSmall`, `rLarge` | $r_<$, $r_>$ | $a_0$ | arguments of `yukawaKernel` | — |
| `rhoBd` | $\rho_{bd} = P_bP_d + Q_bQ_d$ | $a_0^{-1}$ | local in the sweeps | the source of $V_k$ |
| `rhot` | $\rho_t = \sum_a \mathrm{occ}_a(P_a^2+Q_a^2)$ | $a_0^{-1}$ | local in `computePotential` | the source of the mean field |
| `Vk` | $V_k(r)$ | — | returned by the sweeps | contracted against $\rho_{ac}$ |
| `lambda` (argument of `SlaterRkDebyeHueckel`) | **$\mu$, not $\lambda_D$** | $a_0^{-1}$ | pre-existing signature | see the warning below |

### The four ambiguities this section exists to remove

1. **$\lambda_D$ versus $\mu$.** $\lambda_D$ is a *length* in $a_0$; $\mu = 1/\lambda_D$ is an *inverse length*
   in $a_0^{-1}$. The user API takes $\lambda_D$; every internal function takes $\mu$; the conversion happens in
   exactly one place, `Defaults.eeScreeningMu()`.
   **The one exception is a pre-existing misnomer.** `RadialIntegrals.SlaterRkDebyeHueckel(..., lambda)` has an
   argument *named* `lambda` that has always in fact been $\mu$ — `[CODE EVIDENCE]` its callers
   (`InteractionStrength.XL_Coulomb_DH`, and through it `Basics.compute(..., plasmaModel)` and
   `AutoIonization.computeLinesPlasma`) all pass `1/plasmaModel.debyeLength`. It was **not renamed**: it is a
   public signature and renaming it is outside this task. It is now stated in the first lines of that
   function's docstring so no reader can be misled.
2. **$r$ versus $r_{12}$.** $r$ is a radial grid coordinate, an electron's distance from the nucleus.
   $r_{12}$ is the separation of two electrons and is *not* a coordinate of this representation at all
   (§1.4). No function in this work takes $r_{12}$; the multipole expansion has already removed it.
3. **Radial density versus orbital product.** $\rho_t$ (a sum over occupied orbitals, weighted by occupation)
   is what the *mean field* is built from; $\rho_{bd}$ (one orbital *pair*) is what a *Slater integral* is built
   from. They enter different code paths and must not be confused: the first is `computePotential`'s `rhot`,
   the second is the sweeps' `rhoBd`.
4. **Coulomb versus Yukawa kernel.** `buildScreenedPotential` is the *Coulomb* $Y_k$ construction despite its
   name (§3.1); `buildScreenedPotentialYukawa` is the plasma-screened one.

## 6.4 Default-off behaviour, and why it is an identity test

`[CODE EVIDENCE]` `Defaults.eeScreeningMu()` returns **exactly `0.0`** when no screening is in force, and every
branch tests that value for identity (`mu == 0.` or `mu > 0.`). The consequence is stronger than "the
unscreened answer is recovered": with screening off, **no Yukawa code is entered at all** and the original
lines execute unchanged, so the former numbers are reproduced bit for bit rather than to within a tolerance.
`[NUMERICAL VERIFICATION]` The suite test asserts this with `!=` rather than a tolerance (§14, check 5), and
the regression is 83/83 with no approved reference moved (§14.1).

## 6.5 Invalid parameters

`[CODE EVIDENCE]` Rejected at the point where they are set, not carried into the kernels:

```julia
if  !isfinite(lambda)  ||  lambda <= 0.
    error("Defaults.setDefaults(\"e-e screening\"): the Debye length must be finite and > 0; got " *
          "lambda_D = $lambda a_o.")
end
```

and any plasma model other than `NoPlasmaModel`/`DebyeHueckelModel` is refused with a message naming the type
it received. `yukawaKernel` additionally asserts $r_< \le r_>$ — the kernel is **not symmetric** in its two
arguments, and silently repairing a mis-ordered call would hide a caller's bug behind a wrong number.
`[NUMERICAL VERIFICATION]` $\lambda_D \in \{0, -1, \mathrm{NaN}, \mathrm{Inf}\}$ and a mis-ordered kernel call
all raise (§14, check 8).

---

# 7. Modification-by-Modification Forensic Record

## 7.1 The screening state and its accessor

```
File:            src/module-Defaults.jl
Module:          Defaults
Function/type:   GBL_EE_SCREENING (new global)
                 setDefaults(sa::String, model::Basics.AbstractPlasmaModel) (new method)
                 eeScreeningMu() (new function)
                 getDefaults("e-e screening") (new branch)
Original:        no e-e screening state existed anywhere in JAC
Math object:     mu = 1/lambda_D
Problem:         two independent code paths must agree on one value of mu, and neither is reachable
                 from a settings object (§6.1)
Required:        one authoritative parameter; off by default; invalid input refused at the setter
Modification:    a validated global, its setter, its getter and a single accessor
Why here:        Defaults is where JAC already keeps session-wide settings, and it is the only module
                 both RadialIntegrals and BasicsAZ already depend on
Consumers:       RadialIntegrals.buildScreenedPotential; 5 x Basics.computePotential;
                 2 x InteractionStrength cache key; SelfConsistent.checkEeScreeningIsConsistent
```

```julia
function eeScreeningMu()
    global GBL_EE_SCREENING
    if  typeof(GBL_EE_SCREENING) == Basics.DebyeHueckelModel    return( 1.0 / GBL_EE_SCREENING.debyeLength )
    else                                                        return( 0.0 )
    end
end
```

The `else` branch returns the literal `0.0`, not a computed zero. That is deliberate and is what makes the
identity tests of §6.4 meaningful.

## 7.2 The branch point — the whole of Path A

```
File:            src/module-RadialIntegrals.jl
Function:        buildScreenedPotential(k, b, d, grid; rtol, mtpOut)
Original:        builds the Coulomb V_k(r) by an O(N) forward/backward cell sweep
Math object:     V_k(r) = r^-(k+1) int_0^r s^k rho_bd ds + r^k int_r^inf s^-(k+1) rho_bd ds
Problem:         this is the ONLY radial construction the production Slater path uses, and it was
                 unconditionally Coulomb
Required:        dispatch to the Yukawa sweep when screening is in force, and be untouched otherwise
Modification:    four lines at the top of the function
Why here:        it is the single point every kink-aware consumer funnels through -- SlaterRkKinkAware
                 with and without its V_k cache, and XL_CoulombKinkAwareKernel -- so ONE insertion
                 reaches the direct and the exchange integral alike (§9)
Consumers:       SlaterRkKinkAware -> XL_Coulomb / XL_CoulombKinkAware -> Hamiltonian.setupMatrix
```

```julia
mu = Defaults.eeScreeningMu()
if  mu > 0.
    return( RadialIntegrals.buildScreenedPotentialYukawa(k, b, d, grid, mu; rtol=rtol, mtpOut=mtpOut) )
end
```

Line by line: `eeScreeningMu()` is read **once per potential build**, not once per grid point, so the cost is
negligible; the test is `> 0.` against a value that is exactly `0.0` when unscreened, so the unscreened path is
byte-identical to what it was; and the `return` is unconditional, so there is no possibility of the Coulomb and
Yukawa contributions being added twice.

## 7.3 The screened potential sweep

```
File:            src/module-RadialIntegrals.jl
Function:        buildScreenedPotentialYukawa (new)
Math object:     V_k^Y(r) = mu(2k+1)[ k_k(mu r) int_0^r i_k(mu s) rho ds + i_k(mu r) int_r^inf k_k(mu s) rho ds ]
Problem:         the unscaled moments overflow (i_k) and underflow (k_k) over a real grid
Required:        the SAME O(N) two-pass structure, with no O(N_r^2) double integration
Modification:    both running moments accumulated in exponentially scaled form
Why here:        the Yukawa kernel is still SEPARABLE -- i_k(mu r_<) k_k(mu r_>) is a product of a
                 function of r_< and a function of r_> -- which is exactly the property the sweep needs
Consumers:       SlaterRkYukawa, and buildScreenedPotential's branch
```

Derivation of the recurrences is in §9.2. The point that is not obvious from the code is the accumulator
update:

```julia
acc = acc * exp(-mu * (ri - grid.r[i-1])) +
      RadialIntegrals.cellIntegral(s -> exp(-mu*(ri - s)) * iScaledYukawa(k, mu*s) * splBd(s),
                                   grid.r[i-1], ri)
```

The scaled moment is $\widetilde I(r) = e^{-\mu r}\int_0^r i_k(\mu s)\rho(s)\,ds$. Advancing from $r_{i-1}$ to
$r_i$ multiplies the accumulated part by $e^{-\mu(r_i - r_{i-1})}$ — a factor $\le 1$ — and adds one cell whose
integrand carries $e^{-\mu(r_i - s)} \le 1$ for $s \in [r_{i-1}, r_i]$. **Every factor formed is bounded by 1
times something of order the answer**, which is what makes the sweep overflow-free without any clamp.

**The tail beyond the source extent.** `[CODE EVIDENCE]` The Coulomb branch continues $V_k$ past the last
tabulated point of the density as `fullInner / r^(k+1)` — a **power law**. The Yukawa branch must not reuse it:

```julia
if  r >= rmaxBd    Vk[i] = pre * kScaledYukawa(k, x) * exp(-mu*(r - rmaxBd)) * fullInnerScaled
```

`[INFERENCE]` A power-law tail attached to screened data would restore precisely the long-range interaction the
screening removes — the multipole field of the enclosed charge, undamped, out to the box edge.
`[NUMERICAL VERIFICATION]` Measured on a deliberately truncated source (§14.8): the decay is exponential, with
a chord slope of $-0.879$ against $-0.879$ predicted for $e^{-\mu r}/(\mu r)$, and an apparent power-law
exponent of $-11.085$ where a Coulomb $k=0$ tail would give $-1.000$.

`[CODE EVIDENCE]` The same hazard was identified independently in FAC, whose screened branch bypasses its
cached $Y_k$ with the comment that *"the cached tail model yk ~ (r0/r)^k is invalid for a screened kernel"*.

## 7.4 The five mean-field potentials — the whole of Path C

```
File:            src/module-BasicsAZ-inc-compute.jl
Functions:       Basics.computePotential for CHField, HSField, KSField, DFSField(level), DFSField(basis)
Original:        each forms  int rho_t(r')/r_>  with its own i,j double loop over the grid
Math object:     the k = 0 direct (Hartree) potential
Problem:         this is a SECOND, independent implementation of the same Coulomb physics; screening
                 the Slater integrals alone leaves the ORBITALS unscreened (§3.2)
Required:        the same k = 0 kernel, with the same mu
Modification:    1/r_> replaced by yukawaKernel(0, r_<, r_>, mu); mu == 0. keeps the original line
Why here:        it is what SelfConsistent.solveMeanFieldBasis iterates the orbitals against
Consumers:       Basics.add(nuclearPotential, wp) -> Bsplines.setupLocalMatrix -> new orbitals
```

```julia
mu = Defaults.eeScreeningMu()
for i = 1:npoints
    if  mu == 0.
        for j = 1:npoints    rg = max( grid.r[i],  grid.r[j]);    wx[j] = rhot[j]/rg   end
    else
        for j = 1:npoints
            rl = min( grid.r[i],  grid.r[j]);    rg = max( grid.r[i],  grid.r[j])
            wx[j] = rhot[j] * RadialIntegrals.yukawaKernel(0, rl, rg, mu)
        end
    end
    wb[i] = RadialIntegrals.V0(wx, npoints, grid::Radial.Grid)
```

The unscreened branch is the original line, character for character. `HSField` differs: it builds its field
from `RadialIntegrals.Yk_ab` rather than from an inline loop, so it receives $\mu$ through a new method
`Yk_ab(k, r, rho_ab, mtp, grid, mu)` which falls through to the unscreened method when `mu <= 0.`.

## 7.5 The caches

```
File:            src/module-InteractionStrength.jl
Functions:       XL_Coulomb(..., cache::XLCache), XL_CoulombKinkAware(..., cache::XLCache)
Original:        key = (:Coulomb, L, a.subshell, b.subshell, c.subshell, d.subshell, 0., :direct)
Problem:         the Float64 slot held a CONSTANT 0.; a value computed under one screening model
                 could be returned under another
Required:        the key must identify the interaction that produced the value
Modification:    the constant 0. replaced by Defaults.eeScreeningMu()
Why here:        this is the only memo in the two-electron path that can outlive a change of mu
Consumers:       Hamiltonian.setupMatrix and setupMatrixKinkAware
```

`[INFERENCE]` The scenario is concrete rather than theoretical: a plasma-shift calculation computes the
field-free multiplet first and the screened one second, in one session. At $\mu = 0$ the key is byte-for-byte
the former one, so no unscreened behaviour changes. `[NUMERICAL VERIFICATION]` §14 and §11.

## 7.6 The consistency guard and the run-time notice

```
File:            src/module-SelfConsistent.jl
Function:        checkEeScreeningIsConsistent(settings; printout) (new)
Problem:         (a) some SCF fields and some interaction operators cannot be screened by this work;
                 (b) a screened run left no trace in its own output
Required:        refuse rather than half-screen; and say plainly when screening is in force
Modification:    a guard raising a naming error, plus a one-line notice
Why here:        performSCF is the single entry point of both the ordinary and the tested route (§3.3)
Consumers:       both performSCF methods
```

`[INFERENCE]` The reasoning for refusing rather than half-screening: a partly screened calculation runs,
converges, and returns plausible energies, with nothing in the output to say that half of the interaction was
left bare. Refusal converts a silent wrong answer into a loud one. `[NUMERICAL VERIFICATION]` All five
refusals fire, and the same settings are still accepted with screening off (§14, check 8).

---

# 8. Yukawa/Bessel Implementation

## 8.1 Why the functions are computed in scaled form

`[PHYSICS]` The asymptotics are $i_k(x) \sim e^{x}/(2x)$ and $k_k(x) \sim e^{-x}/x$. `[INFERENCE]` Over a real
JAC grid, $\mu r$ can reach $10^2$–$10^5$: $i_k$ overflows Float64 beyond $x \approx 709$ and $k_k$ underflows
below $x \approx -745$, so **neither may be formed on its own**. Their *product* is perfectly well scaled,
because the exponentials cancel except for the residual $e^{-\mu(r_> - r_<)} \le 1$.

The implementation therefore works throughout in

$$\hat\imath_k(x) = e^{-x} i_k(x), \qquad \hat k_k(x) = e^{+x} k_k(x),$$

both of which are $O(1/x)$ at large $x$, and carries the residual explicitly:

$$U_k = \mu(2k+1)\,\hat\imath_k(\mu r_<)\,\hat k_k(\mu r_>)\,e^{-\mu(r_> - r_<)}.$$

`[INFERENCE]` This removes the need for any clamp on the argument. FAC's implementation of the same expansion
clamps $\mu r$ at 400 and justifies the clamp physically ("bound orbital densities vanish long before
$\mu r$ reaches _EE_XMAX"); the scaled formulation used here needs no such assumption, which is a genuine
advantage of the present design and is stated here because it is the kind of difference that is easy to lose.

## 8.2 $\hat k_k$ — exact, finite, no branch

$$\hat k_k(x) = \frac{1}{x}\sum_{q=0}^{k} \frac{(k+q)!}{q!\,(k-q)!\,(2x)^q}.$$

`[ANALYTICAL VERIFICATION]` The sum is **finite** (it terminates at $q=k$) and every term is **positive**.
There is therefore no truncation error, no branch, and no cancellation whatever: the result is exact to
rounding at every $x$. The coefficients are built by the recurrence
$c_{q+1} = c_q (k+q+1)(k-q)/(q+1)$, which avoids forming factorials.

The normalisation check is immediate: at $k=0$ the sum has one term, $c_0 = 1$, so $\hat k_0(x) = 1/x$, i.e.
$k_0(x) = e^{-x}/x$ — the required convention, and $2/\pi$ times the GSL/`SpecialFunctions` one (§2.3).

## 8.3 $\hat\imath_k$ — two branches, and why the crossover sits at 35

**Below $x = 35$: the ascending series.**

$$i_k(x) = \sum_{p=0}^{\infty}\frac{x^{k+2p}}{2^p\,p!\,(2k+2p+1)!!}$$

`[ANALYTICAL VERIFICATION]` Every term is positive, so there is **no cancellation**, and the sum converges once
$p \gtrsim x/2$. The implementation advances by the term ratio $x^2/\big(2p(2k+2p+1)\big)$ and stops when a
term falls below $10^{-18}$ of the running total, with a hard cap of 400 iterations (reached only above
$x \approx 800$, which this branch never sees). At the top of its range $i_k(35) \approx 2\times10^{13}$ — far
from overflow.

**Above $x = 35$: the exact closed form.**

$$i_k(x) = \frac{1}{2x}\left[e^{x}\sum_{q=0}^{k}\frac{(-1)^q c_q}{(2x)^q}
+ (-1)^{k+1} e^{-x}\sum_{q=0}^{k}\frac{c_q}{(2x)^q}\right],
\qquad c_q = \frac{(k+q)!}{q!\,(k-q)!},$$

which in scaled form is

$$\hat\imath_k(x) = \frac{1}{2x}\left[\sum_q \frac{(-1)^q c_q}{(2x)^q}
+ (-1)^{k+1} e^{-2x}\sum_q \frac{c_q}{(2x)^q}\right].$$

`[ANALYTICAL VERIFICATION]` The first sum **alternates**, so it is the only place in this implementation where
cancellation is possible at all. It loses digits only while its largest term exceeds the result, i.e. while
$(2k)!/k! > (2x)^k$, which for $k \le 12$ means $x < 8$.

### Why the threshold is 35 and why that is safe

`[INFERENCE]` The two branches are individually accurate on overlapping ranges: the series is exact up to at
least $x \approx 800$ (limited only by the iteration cap), and the closed form is accurate above $x \approx 8$.
The crossover is placed at 35 — well inside both — **so that neither branch is ever used near its own limit**.
That is the whole justification, and it is a deliberate margin rather than a tuned value.

`[NUMERICAL VERIFICATION]` Measured at the crossover against a 256-bit reference, each branch at its own
argument (§14.8): the series gives $\le 8.81\times10^{-16}$ and the closed form $\le 3.28\times10^{-16}$
relative error, for $k = 0 \ldots 8$. There is no step.

## 8.4 Every numerical threshold, and why each is safe

| Threshold | Value | What it does | Why it is safe |
|---|---|---|---|
| `GBL_YUKAWA_X_SWITCH` | 35.0 | series ↔ closed form for $\hat\imath_k$ | both branches are accurate well beyond it in their own direction; measured error $\le 8.8\times10^{-16}$ on both sides (§14.8). **Carries a warning in the source** — see §12, problem 7 |
| `GBL_YUKAWA_X_FLOOR` | 1.0e-25 | below this $\mu r_>$, return the Coulomb kernel | the two kernels differ by $O(\mu r_>) \le 10^{-25}$ relative there, which is $10^{-9}$ times the Float64 epsilon — the substitution is exact as far as Float64 can tell |
| `yukawaXFloor(k)` | 1e-25 for $k \le 10$, rising above | the same floor, made rank-aware | **added in the final review**; the constant floor is safe only to rank 10 — see §12, problem 8 |
| series relative cut | 1e-18 | stops the ascending series | two orders below Float64 epsilon; all terms positive so the truncated tail is bounded by the last term |
| series iteration cap | 400 | hard stop | reached only above $x \approx 800$; this branch is used only below 35 |

`[INFERENCE]` **No threshold in this table was introduced to make a test pass.** Each was derived from the
Float64 exponent range or from the convergence of a series, and each is stated with the range over which it was
verified. Where a threshold's justification turned out to be incomplete — `GBL_YUKAWA_X_FLOOR` — it was
corrected during the review rather than re-argued (§12, problem 8).

## 8.5 Underflow, and what the correct answer is

`[NUMERICAL VERIFICATION]` Of 7830+ sampled kernel evaluations, a number return an exact `0.0` where the
Coulomb kernel is finite. These are **not failures**: they are $e^{-(x_> - x_<)}$ underflowing when the two
electrons are separated by more than about 745 Debye lengths, where the true screened interaction is below
$10^{-323}$. Zero is the correct Float64 answer there. This is the only regime in which the implemented kernel
and the arbitrary-precision reference disagree, and the disagreement is a property of Float64, not of the code.

## 8.6 The independent references used to validate all of this

| Reference | What it is | Shares with the implementation |
|---|---|---|
| `SpecialFunctions.besselix` / `besselkx` | library cylindrical Bessel, converted by $\sqrt{\pi/2x}$ and $\sqrt{2/\pi x}$ | nothing but the mathematics; **note the $\pi/2$** |
| 256-bit `BigFloat` re-implementation | ascending series below $x=60$, upward recurrence $\hat\imath_{k+1} = \hat\imath_{k-1} - \frac{2k+1}{x}\hat\imath_k$ above, seeded by two exact closed forms | nothing; written separately in `tools/diag-debye-screening-stability.jl` |
| Angular Legendre projection | $U_k = \frac{2k+1}{2}\int_{-1}^{1}P_k(w)\,\frac{e^{-\mu r_{12}}}{r_{12}}\,dw$ by QuadGK | **no Bessel function of any kind** — the strongest of the three, because it tests the expansion itself |
| $k=0$ closed form | $\sinh(\mu r_<)e^{-\mu r_>}/(\mu r_< r_>)$ | no Bessel function; pins the normalisation absolutely |

---

# 9. Slater Integral Modification

## 9.1 What was replaced

`[CODE EVIDENCE]` The production Slater integral is

$$R^k(abcd) = \int_0^\infty\!\!\int_0^\infty \rho_{ac}(r)\;U_k(r,s)\;\rho_{bd}(s)\,dr\,ds,
\qquad \rho_{xy} = P_xP_y + Q_xQ_y,$$

evaluated by JAC as a contraction of $\rho_{ac}$ against the potential $V_k$ built from $\rho_{bd}$. The
substitution is confined to $U_k$:

$$U_k(r,s) = \frac{r_<^k}{r_>^{k+1}} \;\longrightarrow\; \mu(2k+1)\,i_k(\mu r_<)\,k_k(\mu r_>).$$

## 9.2 The new radial expression, derived

`[ANALYTICAL VERIFICATION]` Splitting the inner integral at $s = r$ and using $r_< = \min$, $r_> = \max$:

$$V_k(r) = \int_0^{r} U_k(r,s)\rho(s)\,ds + \int_r^{\infty} U_k(r,s)\rho(s)\,ds.$$

In the first integral $s < r$, so $r_< = s$ and $r_> = r$; in the second the roles swap. For Coulomb this gives
JAC's existing two-region form

$$V_k(r) = \frac{1}{r^{k+1}}\int_0^{r} s^k \rho(s)\,ds + r^k\int_r^{\infty} \frac{\rho(s)}{s^{k+1}}\,ds,$$

and for Yukawa, identically,

$$\boxed{\;V_k^{Y}(r) = \mu(2k+1)\left[k_k(\mu r)\int_0^{r} i_k(\mu s)\rho(s)\,ds
+ i_k(\mu r)\int_r^{\infty} k_k(\mu s)\rho(s)\,ds\right].\;}$$

**The structural point is that this is still separable.** In each region the kernel factorises into a function
of $r$ times a function of $s$, so the two inner integrals are *running moments*: the value at $r_i$ differs
from the value at $r_{i-1}$ by exactly one grid cell. That is what preserves the $O(N)$ forward/backward sweep
and is why no $O(N_r^2)$ double integration appears. `[NUMERICAL VERIFICATION]` §14.10 confirms the cost ratio
is flat over a 5.4× range of grid sizes.

The scaled recurrences actually implemented are

$$\widetilde I(r_i) = e^{-\mu\Delta r}\,\widetilde I(r_{i-1}) + \int_{\text{cell}} e^{-\mu(r_i - s)}\hat\imath_k(\mu s)\rho(s)\,ds,$$
$$\widetilde O(r_i) = e^{-\mu\Delta r}\,\widetilde O(r_{i+1}) + \int_{\text{cell}} e^{-\mu(s - r_i)}\hat k_k(\mu s)\rho(s)\,ds,$$

with $V_k^Y(r_i) = \mu(2k+1)\big[\hat k_k(\mu r_i)\widetilde I(r_i) + \hat\imath_k(\mu r_i)\widetilde O(r_i)\big]$.

## 9.3 All multipole orders

`[CODE EVIDENCE]` The rank $k$ is a plain argument threaded from the angular coefficient (`coeff.nu`) through
`XL_Coulomb` to `SlaterRkKinkAware` to `buildScreenedPotential`. Nothing in the substitution is special to
$k = 0$; `iScaledYukawa` and `kScaledYukawa` take $k$ as their first argument and `yukawaKernel` carries the
$(2k+1)$ factor. `[NUMERICAL VERIFICATION]` Ranks 0–16 are exercised in the stability sweep and ranks 0–20 in
the suite test (§14.8).

## 9.4 Direct and exchange — the proof from source

This is the question the audit asks to be *proved*, not asserted.

`[CODE EVIDENCE]` `InteractionStrength.XL_Coulomb` (in `src/module-InteractionStrength.jl`) computes the
angular prefactor and then makes exactly one radial call:

```julia
xc = AngularMomentum.CL_reduced_me(a.subshell, L, c.subshell) *
     AngularMomentum.CL_reduced_me(b.subshell, L, d.subshell)
if   rem(L,2) == 1    xc = - xc    end
XL_Coulomb = xc * RadialIntegrals.SlaterRkKinkAware(L, a, b, c, d, grid)
```

There is **no branch on direct versus exchange anywhere in this function**, and no second radial routine. The
distinction lives entirely in *which orbitals the angular coefficient places in the four slots*:

| | slot pattern | name |
|---|---|---|
| direct | $(a,b,a,b)$ | $F^k(ab) = R^k(abab)$ |
| exchange | $(a,b,b,a)$ | $G^k(ab) = R^k(abba)$ |

`[CODE EVIDENCE]` Those slots are supplied by `SpinAngular.computeCoefficients(TwoParticleOperator(0,plus),
…)` as `coeff.a`, `coeff.b`, `coeff.c`, `coeff.d`, and `Hamiltonian.setupMatrix` passes them straight through.

`[INFERENCE]` Therefore a change to `SlaterRkKinkAware`'s kernel necessarily reaches both, and it is *not
possible* to screen one without the other through this path.

`[NUMERICAL VERIFICATION]` Confirmed three independent ways (§14.6): $F^0(1s,2s)$, $G^0(1s,2s)$ and
$G^1(1s,2p_{3/2})$ all move under the switch; the He $1s2s$ singlet–triplet splitting — which is $2G^0(1s,2s)$
and contains **no direct contribution at all** — moves from 3.4558e-2 to 7.3395e-2 Ha; and the Ne production
run carries 45 pairs including $G^1(1s,2p)$ and $G^2(2p,2p)$.

## 9.5 Why the angular coefficients need not change

`[PHYSICS]` The substitution replaces the *coefficient of* $P_k(\cos\omega)$ in the multipole expansion and
leaves $P_k(\cos\omega)$ itself untouched. Every angular-momentum factor — the reduced matrix elements
$\langle a\|C^{(L)}\|c\rangle$, the Wigner coefficients, the $(-1)^L$ phase, the triangular-delta selection
rules — depends only on the rank $k$ and on the subshell quantum numbers, never on the radial form of the
kernel.

`[CODE EVIDENCE]` This is visible directly in the code fragment above: `xc` (angular) and
`SlaterRkKinkAware(...)` (radial) are computed by separate calls and multiplied. The screening enters the
second factor only. `[NUMERICAL VERIFICATION]` The suite's pre-existing selection-rule checks — that
$X^L$ vanishes *exactly* where parity forbids it — continue to pass, screened and unscreened alike.

---

# 10. SCF / Hartree Modification

## 10.1 Was the SCF modified? Yes, and it had to be

`[CODE EVIDENCE]` JAC's default `scField = DFSField()` selects
`Basics.scfProcedure(DFSField) == :meanFieldIteration`, which runs
`SelfConsistent.solveMeanFieldBasis[Anderson]`: it builds a local potential from the current density,
diagonalises, and iterates. **The orbitals genuinely move.** This is not a frozen-orbital calculation, so a
frozen-orbital implementation of the screening would have been wrong for the default settings.

## 10.2 What density, and what is integrated

`[CODE EVIDENCE]` The mean field is built from the *total* electron density

$$\rho_t(r) = \sum_a \mathrm{occ}_a\,\big(P_a(r)^2 + Q_a(r)^2\big),$$

accumulated inline as `rhot` in each `computePotential` method (`occ` from
`Basics.computeMeanSubshellOccupation`). The quantity integrated is

$$V^{\text{dir}}(r) = \int_0^\infty \rho_t(r')\,U_0(r,r')\,dr',$$

evaluated as an explicit $i,j$ double loop over the grid with `RadialIntegrals.V0` performing the quadrature.

## 10.3 The $k=0$ Yukawa potential in JAC's own normalisation

`[ANALYTICAL VERIFICATION]` At $k = 0$ the substitution reads

$$\frac{1}{r_>} \;\longrightarrow\; \mu\, i_0(\mu r_<)\, k_0(\mu r_>)
= \mu\,\frac{\sinh(\mu r_<)}{\mu r_<}\cdot\frac{e^{-\mu r_>}}{\mu r_>}
= \frac{\sinh(\mu r_<)\,e^{-\mu r_>}}{\mu\,r_<\,r_>},$$

so the screened direct potential is

$$V^{\text{dir}}_{\text{scr}}(r) = \int_0^\infty \rho_t(r')\,
\frac{\sinh\big(\mu \min(r,r')\big)\,e^{-\mu \max(r,r')}}{\mu\,r\,r'}\,dr'.$$

Two sanity checks on this expression. As $\mu \to 0$, $\sinh(\mu r_<) \to \mu r_<$ and $e^{-\mu r_>} \to 1$,
giving $1/r_>$ — the Coulomb form. And for $r$ far outside the charge distribution it behaves as
$e^{-\mu r}/(\mu r)$ times the enclosed charge, i.e. a screened monopole rather than a bare one.

`[CODE EVIDENCE]` In the code this is `RadialIntegrals.yukawaKernel(0, rl, rg, mu)` with `rl = min`,
`rg = max`; and in the `HSField` method it arrives instead through the new
`RadialIntegrals.Yk_ab(0, r, rho, mtp, grid, mu)`. `[NUMERICAL VERIFICATION]` The $k=0$ closed form is checked
against `yukawaKernel` directly in the suite test to $\le 10^{-11}$ relative, and to $2.06\times10^{-16}$ in
the kernel diagnostic — this is the check that pins the normalisation without any Bessel function (§2.3).

## 10.4 How the screened potential enters the orbital equations

`[CODE EVIDENCE]` The chain is

```
computePotential -> Basics.add(nuclearPotential, wp) -> Bsplines.setupLocalMatrix -> generalized eigenproblem
                 -> new orbitals -> new rho_t -> (iterate to self-consistency)
```

`[INFERENCE]` The screening therefore changes the *shape* of the converged orbitals, not merely the energy
evaluated on fixed ones. Physically the electrons see less mutual repulsion, so each is bound more tightly by
the (unscreened) nucleus and the orbitals contract.

## 10.5 Why this is necessary for self-consistency

`[NUMERICAL VERIFICATION]` The size of the effect was measured rather than assumed. Screening everything
against screening only the CI matrix on Coulomb-converged orbitals, Be $1s^2 2s^2$:

| $\mu$ [$a_0^{-1}$] | full (SCF+CI) [Ha] | CI only, frozen [Ha] | difference [Ha] | % of $\Delta E$ |
|---|---|---|---|---|
| 0.010 | −14.630234506802 | −14.630219463692 | −1.504e−05 | 0.0254 |
| 0.100 | −15.104234227457 | −15.102385763470 | −1.848e−03 | 0.3466 |
| 0.500 | −16.458003379469 | −16.376148172566 | −8.186e−02 | 4.3378 |
| 1.000 | −17.390701542441 | −17.135405589875 | −2.553e−01 | **9.0539** |

The relaxation is second order in $\mu$, as §2.5 requires, and reaches 9 % of the whole shift at
$\mu = 1\,a_0^{-1}$. A Slater-only implementation would have been wrong by that amount, silently.

## 10.6 A change that was investigated and rejected

`[CODE EVIDENCE]` `Basics.computePotential` has **eight** methods. Five were modified; three were deliberately
not:

| Method | Decision | Reason |
|---|---|---|
| `AaHSField`, `AaDFSField` | **not modified** | the average-atom plasma line, out of the scope set for this work — and their own `mu` argument is the **chemical potential**, a different quantity entirely. A patch that touched `AaDFSField` was written, detected and reverted; see §12, problem 3 |
| `ThomasFermiField` | **not modified** | it takes `(Z, noElectrons)` and no orbitals: it is a statistical *starting* potential, not a Hartree term built from the density. Thomas–Fermi screening was explicitly out of scope |

`[NUMERICAL VERIFICATION]` The five that were modified are all exercised screened, with the correct sign, by
the suite test (§14, check 10) — added during the final review precisely because a production run exercises
only whichever one its `AsfSettings` names.

---

# 11. Cache and Performance Changes

## 11.1 What was cached before, and why the old key was insufficient

`[CODE EVIDENCE]` `InteractionStrength.XLCache` memoises the two-electron matrix element $X^L(abcd)$ for the
duration of one CI matrix. Its key was

```julia
key = (:Coulomb, L, a.subshell, b.subshell, c.subshell, d.subshell, 0., :direct)
```

The `Float64` slot held a **literal constant `0.`**. `[INFERENCE]` A key must identify everything the value
depends on. Once the interaction itself can change within a session, a key that does not name the interaction
can return a value computed under a different Hamiltonian — and the failure is silent, because a cached number
is indistinguishable from a computed one. The concrete scenario is not hypothetical: a plasma-shift study
computes the field-free multiplet first and the screened one second, in one Julia session.

The change is one token in each of two functions:

```julia
key = (:Coulomb, L, a.subshell, b.subshell, c.subshell, d.subshell, Defaults.eeScreeningMu(), :direct)
```

`[NUMERICAL VERIFICATION]` One cache object carried across a change of $\mu$ (§14.8):

```
mu = 0.0  (fresh)        X^0 = 4.208057097553249
mu = 0.5  (same cache)   X^0 = 3.346223116076047      differs from the first: true
mu = 0.0  (same cache)   X^0 = 4.208057097553249      identical (===) to the first: true
cache holds 2 entries
```

Without $\mu$ in the key the second line would have repeated the first. At $\mu = 0$ the key is byte-for-byte
the former one, so no unscreened behaviour changes.

## 11.2 Can Coulomb and Yukawa results coexist safely?

`[NUMERICAL VERIFICATION]` Yes, and the test above is exactly that question: the two values coexist in one
`Dict` as two entries and each is returned for its own $\mu$.

## 11.3 The complete cache audit

`[CODE EVIDENCE]` Every memo in or near the two-electron path, and its verdict:

| Cache | Where constructed | Lifetime | Carries $\mu$? | Verdict |
|---|---|---|---|---|
| `InteractionStrength.XLCache` | `BasicsAZ-inc-compute.jl:38`, `Hamiltonian.jl:111,148,365` | one CI matrix, but a plasma run builds two per session | **now yes** | fixed; verified §11.1 |
| `vkCache` — $V_k$ vectors, keyed `(k, subshell_b, subshell_d, mtp)` | `SelfConsistent-inc-orbitals.jl:97,145` | **per call, by contract** | not needed | keyed on subshell *labels*, so its own docstring already forbids it outliving the orbitals; a change of $\mu$ forces re-convergence, hence a new cache. Verified exact within one $\mu$ (§14.8) |
| `ScreenedPotentialCache` / `tensorCaches` | `SelfConsistent-inc-averagelevel.jl:139`, `-optimizedlevel.jl:1600` | per SCF | n/a | the B-spline tensor line, used only by `ALField`/`EOLField` — **both refused** under screening (§17) |
| `rkCache`, `radial1pCache`, `radial2pCache` | `SelfConsistent-inc-orbitals.jl:87,144`, `-optimizedlevel.jl:1565` | per call | n/a | EOL only; refused |
| `Dierckx.Spline1D` | 4 sites in `module-RadialIntegrals.jl` | per build | n/a | interpolates the **density** $\rho_{bd}$, not the kernel; the Yukawa sweep rebuilds it identically |

## 11.4 The Coulomb-specific asymptotic extrapolation that had to be bypassed

`[CODE EVIDENCE]` This was the most dangerous item in the audit, because it is not a cache in the usual sense
and would not be found by grepping for "cache". Beyond the last tabulated point of the source density, the
Coulomb branch of `buildScreenedPotential` continues $V_k$ analytically as

$$V_k(r) = \frac{\text{fullInner}}{r^{k+1}}, \qquad r > r_{\max},$$

the exterior multipole field of the enclosed charge — a **power law**. `[INFERENCE]` Reusing it under screening
would reinstate exactly the long-range interaction the screening exists to remove: the potential would fall as
$r^{-(k+1)}$ out to the box edge instead of being cut off at the Debye length. The Yukawa branch therefore
carries its own tail,

$$V_k^Y(r) = \mu(2k+1)\,\hat k_k(\mu r)\,e^{-\mu(r - r_{\max})}\,\widetilde I_{\text{full}},$$

i.e. $\propto e^{-\mu r}/(\mu r)$ at $k=0$.

`[NUMERICAL VERIFICATION]` Measured on a deliberately truncated source (the density cut at half the grid, so a
tail region exists at all), $r = 7.243 \ldots 20.127\,a_0$ beyond an extent of $7.236\,a_0$, $\mu = 0.8$:

```
apparent power-law exponent  d ln V / d ln r =  -11.085     (a Coulomb k=0 tail would give -1.000)
chord decay constant         d ln V / d r    =   -0.879
predicted for exp(-mu r)/(mu r)              =   -0.879
```

The predicted chord slope is $-\mu - \ln(r_3/r_1)/(r_3-r_1)$; the $1/r$ prefactor contributes the second term,
and comparing against $-\mu$ alone would have reported a correct tail as wrong by 10 % (§12, problem 9).

`[CODE EVIDENCE]` The identical hazard was identified independently in FAC, whose screened path bypasses its
cached $Y_k$ with the comment that the cached tail model is invalid for a screened kernel. Two independent
implementations arriving at the same conclusion is the strongest evidence available that this is a real trap
rather than an over-cautious reading.

## 11.5 Performance

**Method.** `[NUMERICAL VERIFICATION]` Each kernel/sweep figure is the **minimum of 15 repetitions**, not the
mean. A single sweep is a few milliseconds and roughly one run in four catches a garbage collection; the mean
produced a spurious ratio of 0.74 at one grid size before this was corrected (§12, problem 10). The SCF figures
are single warm runs and are quoted to two significant figures only — they are **not** precision benchmarks and
should not be read as such.

| Calculation | Coulomb | Yukawa | Ratio |
|---|---|---|---|
| Kernel + potential sweep, 812 pts | 2.894 ms | 3.273 ms | 1.13 |
| Kernel + potential sweep, 1330 pts | 6.875 ms | 7.697 ms | 1.12 |
| Kernel + potential sweep, 2338 pts | 18.691 ms | 20.912 ms | 1.12 |
| Kernel + potential sweep, 4354 pts | 61.549 ms | 66.696 ms | 1.08 |
| Full SCF, Be $1s^2 2s^2$, 1582 pts (single run, ±10 %) | ≈ 1.1 s | ≈ 1.9 s | ≈ 1.7 |
| Retired truncated $O(N_r^2)$ form, 1330 pts | — | 13.983 ms | **13.2× the screened sweep** |

**Is the $O(N)$ structure preserved?** `[NUMERICAL VERIFICATION]` Yes. The ratio is 1.08–1.13 and **flat** over
a 5.4× range in $N$. A ratio growing with $N$ would be the signature of an accidental $O(N_r^2)$; it does not
grow.

`[INFERENCE]` A caveat stated because the table shows it: the cost *per 1000 grid points* rises for **both**
kernels together (3.56 → 14.14 ms/1000 pts for Coulomb, 4.03 → 15.32 for Yukawa). That super-linearity belongs
to the pre-existing unscreened algorithm — most likely the spline evaluation and cache behaviour — and is
identical in both. It is outside the scope of this work and is recorded here only so that the flat *ratio* is
not mistaken for a claim that the absolute cost is linear.

The end-to-end SCF factor of roughly 1.7–1.9 is larger than the sweep ratio because it is dominated by
`Basics.computePotential`, whose $i,j$ double loop over the grid was **already** $O(N^2)$ before this work;
screening makes each of those $N^2$ kernel evaluations more expensive (a Bessel series instead of a division)
but does not change the exponent.

---

# 12. Debugging History

Ten problems are recorded. Each cost real time, and each would otherwise be rediscovered by whoever implements
this next. Four of them (4, 5, 6, 9) share one lesson worth stating in advance: **more of these were failures
of the reference than of the implementation**, and each time the temptation was to widen a tolerance.

## Problem 1 — Bessel overflow in the prototype reference

```
1. Expected      a library value for i_k(x) at x = 1000, to check the implementation against.
2. Observed      AmosException with id 2: overflow, at besseli(k+0.5, 1000.0).
3. Suspicious    the implementation had been written specifically to avoid this overflow; the
                 REFERENCE had walked into it.
4. Hypothesis    the reference used the unscaled library call.
5. Diagnostic    read the call: besseli, not besselix.
6. Evidence      i_k(1000) ~ e^1000/2000 ~ 10^431, far outside Float64.
7. Change        moved the reference to besselix / besselkx and removed the manual exp factors.
8. Result        agreement to 2.56e-14 (ihat) and 8.93e-16 (khat) over k = 0,1,2,4,6,8 and
                 x = 1e-8 ... 1000.
9. Lesson        a reference for an overflow-avoiding implementation must itself avoid the overflow.
                 The scaled functions are not an optimisation; they are the only way to write this.
```

## Problem 2 — the Coulomb limit is not second order at every rank

```
1. Expected      the Coulomb limit to be approached at O(mu^2) for all k, and a test written to check it.
2. Observed      k = 0 approached at O(mu^1); the test would have failed a correct implementation.
3. Suspicious    every other rank behaved as expected, so it was not a coding error.
4. Hypothesis    the k = 0 term carries a constant correction.
5. Diagnostic    expand the closed form:
                 mu i_0(mu a) k_0(mu b) = sinh(mu a) e^(-mu b)/(mu a b) ~ 1/b - mu.
6. Evidence      the leading correction is the CONSTANT -mu, so the RELATIVE error is O(mu r_>).
                 Measured orders: 0.996-0.999 for k = 0, 1.995-2.000 for k >= 1.
7. Change        the test now checks the OBSERVED ORDER (1 for k = 0, 2 for k >= 1) rather than
                 agreement at a fixed mu.
8. Result        passes, and now asserts something true rather than something convenient.
9. Lesson        that constant -mu is not an inconvenience; it IS the physics of Delta E = -mu N(N-1)/2.
                 A test that had been "fixed" by loosening its tolerance would have hidden the one
                 feature of the k = 0 kernel that matters most.
```

## Problem 3 — the `mu` name collision in `AaDFSField` (caught before it ran)

```
1. Expected      a scripted edit to insert  mu = Defaults.eeScreeningMu()  into the FIVE e-e mean-field
                 methods of Basics.computePotential.
2. Observed      the patch applied to Basics.computePotential(::AaDFSField, grid, orbitals, mu, temp)
                 -- and left DFSField(basis) unpatched.
3. Suspicious    the count of patched methods was right but the LIST was not; a marker-based
                 replace(old, new, 1) had matched the first 'rhot' in the file.
4. Hypothesis    the methods are ordered differently in the file than in the edit list.
5. Diagnostic    list the ENCLOSING FUNCTION of every occurrence of Defaults.eeScreeningMu().
6. Evidence      AaDFSField's own second argument is NAMED mu -- and it is the CHEMICAL POTENTIAL of
                 the average-atom model, a completely different quantity. The inserted line would have
                 SHADOWED it silently, corrupting an unrelated plasma model, while the method that
                 should have been patched was untouched.
7. Change        reverted the AaDFSField segment to its original text; patched DFSField(basis).
8. Result        the census now lists exactly CHField, HSField, KSField, DFSField x 2.
9. Lesson        two quantities that share a symbol get conflated by whatever edits the file, human or
                 scripted. Verify a bulk edit by listing what it TOUCHED, never by counting how many
                 times it fired. This is the same failure mode CLAUDE.md Rule 18 records for the
                 sqrt(2j+1) factor, arriving by a different route.
```

## Problem 4 — an apparent kernel failure that was the reference

```
1. Expected      the angular Legendre projection to agree with the kernel everywhere.
2. Observed      k = 0, r1 = 1, r2 = 1.0000001, mu = 0.7: relative 1.54e-9, far outside the ~1e-15
                 seen at every other point.
3. Suspicious    the discrepancy appeared only where r1 and r2 nearly coincide, i.e. r12 -> 1e-7.
4. Hypothesis    the QuadGK integrand is near-singular there, so the REFERENCE is inaccurate.
5. Diagnostic    add a THIRD, exact route for k = 0: sinh(mu r_<) e^(-mu r_>)/(mu r_< r_>).
6. Evidence      the closed form agrees with the kernel to 2.06e-16 at that very point.
7. Change        removed the near-coincident point from the quadrature-reference list and kept the
                 closed-form check, which is exact there.
8. Result        both references agree with the implementation to rounding.
9. Lesson        THE FIX WAS NOT TO LOOSEN THE TOLERANCE. When a reference and an implementation
                 disagree, the reference is a suspect too, and the way to settle it is a third route.
```

## Problem 5 — cancellation in the test's own reference formula

```
1. Expected      ihat_0(x) = (1 - exp(-2x))/(2x) as the reference in the new suite test.
2. Observed      at x = 1e-6 the test FAILED: implementation 0.9999990000006668,
                 reference 0.9999989999953662, relative 5.3e-12.
3. Suspicious    every other x passed to 1e-15; only the smallest failed.
4. Hypothesis    catastrophic cancellation: exp(-2e-6) = 0.999998..., so 1 - exp(-2x) subtracts two
                 numbers agreeing in eleven digits.
5. Diagnostic    recompute with expm1: -expm1(-2x)/(2x).
6. Evidence      that gives 0.9999990000006668 -- exactly the implementation's value.
7. Change        the reference now uses expm1. The tolerance was NOT touched.
8. Result        the check passes at every x from 1e-6 to 400.
9. Lesson        problem 4 again, in a different disguise, and inside a test written by someone who had
                 just been bitten by it. A reference formula is code and needs the same numerical care
                 as the thing it tests.
```

## Problem 6 — the 256-bit reference truncated its own series

```
1. Expected      an arbitrary-precision reference to be right by construction.
2. Observed      at x = 1e5 it returned 0 while the implementation returned 5e-11, i.e. "infinitely
                 wrong".
3. Suspicious    5e-11 is exactly 1/(2x^2) = the correct value of mu i_0(x) k_0(x) at x_< = x_> = 1e5;
                 the IMPLEMENTATION looked right and the reference did not.
4. Hypothesis    the reference's ascending series was capped at 20000 terms while it needs ~x/2.
5. Diagnostic    at x = 1e5 the series needs ~50000 terms; 20000 truncates it far from convergence.
6. Evidence      raising the cap changed the reference and not the implementation.
7. Change        above x = 60 the reference now uses the upward recurrence
                 ihat_(k+1) = ihat_(k-1) - (2k+1)/x ihat_k, seeded by two exact closed forms and stable
                 in that regime because x >> k.
8. Result        worst error 7.58e-15 over 210 combinations spanning mu*r = 1e-12 to 1e5, k = 0..16.
9. Lesson        the third time the reference was the defect. Arbitrary precision protects against
                 ROUNDING, not against a wrong algorithm or a truncated loop.
```

## Problem 7 — a real sign error, invisible through the public API

```
1. Expected      the closed-form branch of iScaledYukawa to be correct; it agreed with every reference.
2. Observed      nothing. Every test passed. It was found only by lifting the branch out of the module
                 and evaluating it BELOW its own switch.
3. Suspicious    the final review asked whether each branch is correct ON ITS OWN, independently of
                 where it happens to be used.
4. Hypothesis    read the code against the formula:
                    i_k = [e^x sum(-1)^q c_q/(2x)^q + (-1)^(k+1) e^-x sum c_q/(2x)^q]/(2x).
                 The code wrote (-1)^(k+1) as  isodd(k+1) ? 1. : -1.  -- which is the INVERSE:
                 (-1)^odd = -1, not +1.
5. Diagnostic    evaluate the lifted branch at x = 2, 5, 10, 18 against the exact
                 ihat_0 = (1-exp(-2x))/(2x).
6. Evidence      relative error 3.7e-2 at x = 2, 9.1e-5 at x = 5, 4.1e-9 at x = 10, 3.7e-16 at x = 18.
7. Change        isodd(k+1) ? -1. : 1.
8. Result        exact to Float64 (0.0 or ~1e-16 relative) for k = 0 and k = 1 at every x from 2 to 40.
9. Lesson        THE MOST IMPORTANT ENTRY HERE. The branch is only ever entered above x = 35, where the
                 mis-signed term is exp(-70) = 4e-31 -- FIFTEEN orders below the Float64 epsilon. No
                 result JAC had ever produced was affected, and NO TEST THROUGH THE PUBLIC API COULD
                 HAVE CAUGHT IT, because through the API the term is not representable at all.
                 A test that stays inside the operating range cannot validate a branch whose error
                 term is negligible inside that range. The defence is therefore not a test but a
                 WARNING AT THE CONSTANT: GBL_YUKAWA_X_SWITCH now carries a note saying that lowering
                 it requires re-verifying the closed form on its own, and how.
```

## Problem 8 — the constant floor was rank-blind: Inf, then NaN

```
1. Expected      GBL_YUKAWA_X_FLOOR = 1e-25 to make the kernel safe at any rank, as its comment claimed
                 ("only their PRODUCT is well scaled, so for absurdly weak screening the second factor
                 would overflow while the first underflows").
2. Observed      in the final review, at mu*r_> = 1e-24 with r_< = r_>:
                    k = 12 -> Inf,  k = 14 -> NaN,  where the correct answer is 1.0.
3. Suspicious    the stability sweep had reported ZERO violations -- but it stopped at k = 8.
4. Hypothesis    khat_k(x) is dominated at small x by (2k)!/(k! (2x)^k x), which diverges faster the
                 higher the rank; a single constant floor can only be safe up to some k.
5. Diagnostic    tabulate kScaledYukawa(k, 1e-25) for k = 0 ... 16.
6. Evidence      finite to k = 10 (6.5e283); Inf from k = 11. So 1e-25 was safe only to rank 10, and
                 the comment claiming otherwise was wrong.
7. Change        a rank-aware floor, RadialIntegrals.yukawaXFloor(k), derived from the condition that
                 the dominant term reach 1e300:  x^(k+1) 2^k >= (2k)!/(k! 1e300).  It returns exactly
                 1e-25 for every k <= 10, so NOTHING changes for any rank an atomic Slater integral
                 carries (k <= 2 l_max). A matching guard was added to buildScreenedPotentialYukawa,
                 which calls kScaledYukawa directly and so does not inherit the kernel's short-circuit;
                 it RAISES rather than returning Inf.
8. Result        finite and equal to the Coulomb kernel for k = 0 ... 30. The suite test now covers
                 k up to 20 and the stability sweep k up to 16, so the gap cannot reopen silently.
9. Lesson        a test sweep defines the domain over which "no NaN" has been established, and nothing
                 more. The sweep stopped at k = 8 because that is where atomic ranks stop -- which is
                 exactly the reasoning that leaves an edge unexamined. Widen the sweep past the
                 physically reachable range, precisely because that is where nobody is looking.
```

## Problem 9 — three test criteria that were wrong, not the code

Grouped because they share one shape: a check that flagged correct code, and was corrected by making the
comparison honest rather than by loosening it.

```
(a) THE EXPONENTIAL TAIL.
1. Expected      d ln V / d r = -mu beyond the source extent.
2. Observed      -0.879 against mu = 0.800; the check reported "not exponential".
3. Hypothesis    the tail is exp(-mu r)/(mu r), not exp(-mu r): the 1/r prefactor contributes.
4. Evidence      the predicted CHORD slope is -mu - ln(r3/r1)/(r3-r1) = -0.879. Measured -0.879.
5. Change        the criterion now compares against the correct prediction.
6. Lesson        a criterion derived from a simplified model of the code will condemn the code.

(b) THE COULOMB KERNEL'S OWN FLOAT64 LIMIT.
1. Observed      after widening the sweep to k = 16, Inf and NaN at k = 10, x_> ~ 1e-30.
2. Suspicious    the rank-aware floor of problem 8 was supposed to have removed exactly this.
3. Diagnostic    evaluate the UNSCREENED Coulomb kernel r_<^k/r_>^(k+1) at the same points.
4. Evidence      it returns the SAME Inf and the SAME NaN: (1e-30)^11 underflows to 0, so the Coulomb
                 kernel is Inf there whether or not anything is screened, and JAC's pre-existing code
                 has the identical behaviour. r = 1e-30 a_o is twenty orders below the nuclear radius.
5. Change        the sweep now EXCLUDES and COUNTS combinations where the Coulomb kernel is itself
                 non-finite, and says so in its output.
6. Lesson        "the screened kernel returns Inf" is only a finding if the unscreened one does not.

(c) SUBNORMALS.
1. Observed      23 cases where the screened kernel slightly EXCEEDS the Coulomb kernel -- by 8e-8 to
                 7 % -- which is physically impossible.
2. Diagnostic    check whether ihat_k(x_<) is subnormal at those points.
3. Evidence      it is, in every one: ihat_10(5.6e-32) = 2.5e-323, which carries about TWO significant
                 BITS. A 7 % discrepancy is what two bits of mantissa produce.
4. Change        excluded and counted, with the reason printed.
5. Lesson        Float64 has two edges, not one. The overflow edge is famous; the subnormal edge is
                 where a value still exists but its precision has quietly gone.
```

## Problem 10 — a benchmark that measured the garbage collector

```
1. Expected      the Yukawa/Coulomb cost ratio to be a constant slightly above 1.
2. Observed      1.11, 0.74, 1.00, 1.08 -- including a ratio BELOW 1, i.e. screening apparently
                 cheaper than not screening.
3. Suspicious    0.74 is not a plausible physical result, so the measurement was suspect.
4. Hypothesis    a single sweep takes milliseconds; a garbage collection lands in about one run in four
                 and dominates the mean.
5. Change        report the MINIMUM of 15 repetitions -- the run that was not interrupted.
6. Result        1.13, 1.12, 1.12, 1.08. Flat, as the O(N) claim requires.
7. Lesson        do not report a single-shot timing as a benchmark. A ratio below 1 was the tell;
                 without it the noisy numbers would have been written down as fact.
```

## Smaller issues, recorded so they are not rediscovered

`Radial.Grid(true; rnt=...)` takes only `printout` (use `Radial.Grid(Radial.Grid(false); ...)`);
`Radial.Grid` has no `rbox` property (use `grid.r[end]`); `Basics.EOLField` and `CoulombBreit` have no
zero-argument constructors; Julia's `let a, b = f()` binds only `b`; Julia soft-scope warnings inside `for`
loops at top level, fixed with `Ref`.

---

# 13. Validation Strategy

## 13.1 The layered hierarchy

```
   Bessel functions          ihat_k, khat_k
          |                  isolates: NORMALISATION (the pi/2 trap), overflow, series truncation
          v
   Yukawa kernel             U_k(r_<, r_>)
          |                  isolates: the EXPANSION itself, the (2k+1) factor, the Coulomb limit, sign
          v
   radial integrals          V_k(r), the O(N) sweep
          |                  isolates: the SWEEP -- separability, the two-region split, the tail model
          v
   Slater integrals          R^k(abcd)
          |                  isolates: the CONTRACTION against rho_ac, the mtp/extent bookkeeping
          v
   two-electron ME           X^L(abcd)
          |                  isolates: the ANGULAR x RADIAL junction; direct vs exchange slots
          v
   SCF / direct potential    V^dir(r)
          |                  isolates: PATH C -- the second, independent Coulomb implementation
          v
   atomic energy             E(mu)
          |                  isolates: the WHOLE, against a parameter-free analytic law
          v
   full regression           83 / 83
                             isolates: that nothing else moved
```

## 13.2 Why each layer isolates a different failure mode

`[INFERENCE]` The layering is not ceremony; each level catches a class of error that the level above cannot
see, and — importantly — each level *below* the atomic energy can pass while the physics is still wrong.

| Layer | The failure it and only it catches |
|---|---|
| Bessel | A wrong **normalisation**. The $\pi/2$ error is a smooth constant: every higher layer would still converge, still show the right *shape*, still approach the Coulomb limit in the right direction. Only an absolute check against a closed form sees it. |
| Kernel | A wrong **expansion** — the $(2k+1)$, the $\mu$ prefactor, an $r_<$/$r_>$ swap. The angular projection catches these because it never uses a Bessel function; it integrates the actual $e^{-\mu r_{12}}/r_{12}$. |
| Radial integrals | A wrong **sweep**: a dropped cell, a wrong recurrence, a Coulomb tail left attached to screened data. Caught by comparing the $O(N)$ sweep against an explicit $O(N^2)$ double sum over the same kernel — two constructions sharing nothing but the kernel. |
| Slater | Extent (`mtp`) bookkeeping — the "dropped tail" class that once moved an approved reference by 100 %. |
| Two-electron ME | The **angular/radial junction**, and specifically whether exchange is reached. A direct-only implementation passes every layer above this one. |
| SCF | **Path C.** An implementation that screens only the Slater integrals passes every layer above and is still wrong by up to 9 % (§10.5). This is the layer the architecture analysis existed to find. |
| Atomic energy | The **absolute scale and sign**, via a law that involves no orbital, grid or basis. |
| Regression | That the unscreened path is untouched — the one claim no screened test can make. |

## 13.3 Two mathematically independent routes

`[INFERENCE]` A single reference can be wrong, and in this work three of them were (§12, problems 4, 5, 6). The
validation therefore uses references that share **nothing**:

* **Route 1 — special functions.** Library Bessel functions (`besselix`/`besselkx`, with the $\pi/2$
  conversion), and a 256-bit arbitrary-precision re-implementation using a different algorithm.
* **Route 2 — no Bessel functions at all.** The angular projection
  $U_k = \frac{2k+1}{2}\int_{-1}^{1}P_k(w)\,e^{-\mu r_{12}}/r_{12}\,dw$, computed by adaptive quadrature with a
  recurrence-built $P_k$; and the exact $k=0$ closed form $\sinh(\mu r_<)e^{-\mu r_>}/(\mu r_< r_>)$.

Route 2 is the stronger of the two, because it validates the **expansion itself** rather than the evaluation of
its terms. A code that had adopted the $\pi/2$ convention would pass every check in Route 1 that used the same
convention, and fail Route 2 immediately.

---

# 14. Validation Results

Every number in this section was produced by executing code. The commands are in §18. Where a figure is a
timing rather than a computed quantity, it is marked as such and qualified.

## 14.1 Existing regression

`[NUMERICAL VERIFICATION]` Run as `cd test && julia --project=.. runtests.jl` (§18 explains why not
`Pkg.test()`).

| Stage | Result | Runtime | Approved references changed |
|---|---|---|---|
| Baseline, before any edit (`75f963ac`) | **82 / 82** | 9m00.2 s | — |
| After the implementation, screening off | **82 / 82** | 7m42.6 s | none |
| After adding the regression test | **83 / 83** | 7m42.8 s | none |
| **After the final review** (this report) | **83 / 83** | **7m42.1 s** | **none** |
| **The same, re-run to confirm** | **83 / 83** | **7m42.4 s** | **none** |

Zero `Fail`, zero `Error During Test`, zero `Broken` at every stage. The last two rows are two independent
runs of the identical source, made because the report must not rest on a stale measurement: **the latest
modification time under `src/` and `test/` was 02:27:04 and the first of those runs began at 02:35:29**, so the
tree stood still for the whole of it and the result applies to the final source state. The count rose from 82 to 83 because
exactly one test was **added**. `git status --porcelain test/` reports only `test/runtests.jl` — this work's
deliberate four-line addition. **No `test/approved/*.sum` reference was re-approved to make anything pass.**

## 14.2 Bessel and kernel validation

`[NUMERICAL VERIFICATION]`

| Check | Reference | Shares with implementation | Result |
|---|---|---|---|
| $\hat\imath_k$, $\hat k_k$ | `SpecialFunctions` `besselix`/`besselkx` with the $\pi/2$ conversion, $k = 0,1,2,4,6,8$, $x = 10^{-8} \ldots 10^3$ | mathematics only | **2.56e−14** ($\hat\imath$), **8.93e−16** ($\hat k$) |
| $U_k$, 210 combinations, $k = 0,1,2,4,8,12,16$, $\mu r_> = 10^{-12} \ldots 10^{5}$ | 256-bit `BigFloat`, different algorithm | nothing | worst **7.58e−15** (at $k{=}16$, $x{=}35.1$); 21 further combinations underflowed to an exact 0, which is correct in Float64 |
| $U_k$ vs the angular Legendre projection | adaptive quadrature, **no Bessel function at all** | nothing | **~1e−15** |
| $U_0$ vs $\sinh(\mu r_<)e^{-\mu r_>}/(\mu r_< r_>)$ | exact closed form, no Bessel function | nothing | **2.06e−16** |
| Branch continuity at $x = 35$, $k = 0 \ldots 8$ | 256-bit, each branch at its own $x$ | nothing | series **≤ 8.81e−16**, closed form **≤ 3.28e−16**; no step |

The third and fourth rows are the ones that pin the **normalisation** and therefore exclude the $2/\pi$ error
of §2.3, because neither uses a Bessel function of any convention.

## 14.3 Coulomb-limit validation

`[NUMERICAL VERIFICATION]` **At the kernel level**, `yukawaKernel(k, r_<, r_>, 0.0)` is bit-identical to
$r_<^k/r_>^{k+1}$ for $k = 0 \ldots 4$ (tested with `!=`, not a tolerance). The *observed order* of approach,
from the ratio of relative errors at $\mu = 10^{-4}$ and $10^{-5}$:

| $k$ | observed order | expected (§2.4) |
|---|---|---|
| 0 | 0.996 – 0.999 | **1** (the correction is the constant $-\mu$) |
| 1, 2, 3, 4 | 1.995 – 2.000 | **2** |

`[NUMERICAL VERIFICATION]` **At the production level**, on real `perform()` runs of Be $1s^2 2s^2$ ($N = 4$, so
the first-order shift is $-6\mu$):

| $\mu$ [$a_0^{-1}$] | $E$ [Ha] | $\Delta E$ [Ha] | $\Delta E / (-6\mu)$ | residual $\Delta E + 6\mu$ |
|---|---|---|---|---|
| 0 (Coulomb) | −14.570977349023 | — | — | — |
| 1.0e−02 | −14.630234506802 | −5.9257157779e−02 | 0.987619296 | 7.428e−04 |
| 1.0e−03 | −14.576969843895 | −5.9924948720e−03 | 0.998749145 | 7.505e−06 |
| 1.0e−04 | −14.571577273895 | −5.9992487239e−04 | 0.999874787 | 7.513e−08 |
| 1.0e−05 | −14.571037348248 | −5.9999224945e−05 | 0.999987082 | 7.751e−10 |
| 1.0e−06 | −14.570983349030 | −6.0000073940e−06 | 1.000001232 | −7.394e−12 |

The residual falls by **exactly two decades per decade of $\mu$** — second order, as §2.5 requires — until
$\mu = 10^{-6}$, where it reaches the SCF convergence floor and changes sign. The Coulomb limit is recovered.

## 14.4 Weak-screening limit

`[NUMERICAL VERIFICATION]` $\Delta E_{ee}$ against $-\mu N(N-1)/2$ at $\mu = 10^{-5}$, all from real SCF runs:

| System | $N$ | pairs | $\Delta E / [-\mu N(N-1)/2]$ | sign |
|---|---|---|---|---|
| He $1s^2$ | 2 | 1 | **0.999982** | negative ✓ |
| Be $1s^2 2s^2$ | 4 | 6 | **0.999987** | negative ✓ |
| Ne $1s^2 2s^2 2p^6$ | 10 | 45 | **0.999994** | negative ✓ |

Agreement to five significant figures, across a factor 45 in the pair count — which is what tests that the
counting, and not merely the magnitude, is right.

## 14.5 Production atomic calculations

`[NUMERICAL VERIFICATION]` Run through the user-facing route `Atomic.Computation` → `Basics.perform`, grid 1582
points, $r_{\text{box}} = 25\,a_0$. $E$ is the ground level.

| System | Configuration | $N$ | $\mu$ [$a_0^{-1}$] | $E$(Coulomb) [Ha] | $E$(screened) [Ha] | $\Delta E$ [Ha] | $-\mu N(N{-}1)/2$ | $\Delta E$ / prediction |
|---|---|---|---|---|---|---|---|---|
| He | $1s^2$ | 2 | 0.20 | −2.857967131 | −3.034819576 | −0.176852 | −0.200000 | 0.884 |
| Be | $1s^2 2s^2$ | 4 | 0.20 | −14.570977349 | −15.533391190 | −0.962414 | −1.200000 | 0.802 |
| Ne | $1s^2 2s^2 2p^6$ | 10 | 0.20 | −128.672259270 | −136.690308921 | −8.018050 | −9.000000 | 0.891 |
| He | $1s^2$ | 2 | 1.00 | −2.857967131 | −3.444526208 | −0.586559 | −1.000000 | 0.587 |
| Be | $1s^2 2s^2$ | 4 | 1.00 | −14.570977349 | −17.390701542 | −2.819724 | −6.000000 | 0.470 |
| Ne | $1s^2 2s^2 2p^6$ | 10 | 1.00 | −128.672259270 | −156.961366968 | −28.289108 | −45.000000 | 0.629 |

**How to read the last column.** `[PHYSICS]` The weak-screening law is the $\mu \to 0$ limit and is quoted here
as a *scale*, not a tolerance. At $\mu = 0.2$–$1.0\,a_0^{-1}$ the typical $r_{12}$ is of order $1\,a_0$, so
$\mu r_{12}$ is not small and the higher orders of §2.5 are meant to be visible. That the ratio falls
monotonically as $\mu$ grows (0.88, 0.80, 0.89 at $\mu = 0.2$; 0.59, 0.47, 0.63 at $\mu = 1.0$) is the expected
behaviour of a series in $\mu r_{12}$. The quantitative agreement is demonstrated in §14.3 and §14.4, where
$\mu$ is actually small.

**`[INFERENCE]` A caution against over-interpreting these shifts.** Every energy above is a *total* energy, and
the shift is dominated by the trivial $-\mu N(N-1)/2$ piece which is common to all levels. Excitation energies
and transition energies are *differences*, in which that common piece cancels, so they move far less and are
correspondingly more sensitive to grid, basis and convergence. Nothing in this report validates a plasma shift
of a *transition* energy, and §17 records that as untested.

## 14.6 Exchange validation

`[NUMERICAL VERIFICATION]` The decisive calculation is the He $1s2s$ singlet–triplet splitting, because
`[PHYSICS]` to leading order it equals $2G^0(1s,2s)$ — **twice the exchange integral, with no direct
contribution at all**. An implementation that screened only $F^k$ would leave it exactly unchanged.

| Quantity | Coulomb | Screened | Moved? |
|---|---|---|---|
| He $1s2s$ singlet–triplet splitting [Ha] | 3.4558e−02 | 7.3395e−02 | **yes, by 112 %** |
| $F^0(1s,2s)$ | | | yes |
| $G^0(1s,2s)$ | | | yes |
| $G^1(1s,2p_{3/2})$ | | | yes |
| $R^0(1s,2s,1s,2s)$, screened vs unscreened | 4.208057097553249 | 3.346223116076047 | yes, and **smaller**, as screening requires |

The splitting was verified to equal $2G^0(1s,2s)$ exactly at every $\mu$ tested, so the quantity being tracked
is the exchange integral itself and not an accidental combination.

## 14.7 SCF validation

`[NUMERICAL VERIFICATION]` See the table in §10.5. The evidence that the orbitals are genuinely re-optimised is
that screening everything differs from screening only the CI matrix on Coulomb-converged orbitals — by
−1.504e−05 Ha at $\mu = 0.01$ (0.025 % of the shift) rising to −0.2553 Ha at $\mu = 1$ (**9.05 %**). If the
implementation were frozen-orbital, that difference would be exactly zero.

`[NUMERICAL VERIFICATION]` All five screened mean fields respond, and respond in the correct direction — added
to the suite test during the final review, because a production run exercises only whichever one its
`AsfSettings` names:

| Mean field | Route to the screening | Responds with correct sign |
|---|---|---|
| `DFSField(level)` | inline `yukawaKernel(0, …)` loop | ✓ |
| `DFSField(basis)` | inline `yukawaKernel(0, …)` loop | ✓ |
| `HSField` | **`RadialIntegrals.Yk_ab(0, …, mu)`** — a different code path | ✓ |
| `KSField` | inline `yukawaKernel(0, …)` loop | ✓ |
| `CHField` | inline `yukawaKernel(0, …)` loop | ✓ |

## 14.8 Stability

`[NUMERICAL VERIFICATION]`

| Quantity | Value |
|---|---|
| Max relative error vs 256-bit reference | **7.58e−15** (210 combinations, $k \le 16$) |
| Combinations swept for finiteness/sign | **13 544** ($k = 0 \ldots 16$, $\mu r_> = 10^{-30} \ldots 10^{6}$, six ratios) |
| **NaN count** | **0** |
| **Inf count** | **0** |
| Negative values | **0** |
| Values exceeding the Coulomb kernel | **0** |
| Largest $\mu r$ tested | $10^{6}$ (accuracy checked to $10^{5}$) |
| Smallest $\mu r$ tested | $10^{-30}$ |
| Excluded: Coulomb kernel itself non-finite | 1 098 — `(1e-30)^11` underflows, so $r_<^k/r_>^{k+1}$ is `Inf` screened or not |
| Excluded: $\hat\imath_k(x_<)$ subnormal | 148 — e.g. $\hat\imath_{10}(5.6\text{e-}32) = 2.5\text{e-}323$, about two significant bits |
| Underflow to exact 0 | 906 — $e^{-(x_> - x_<)}$ below $10^{-323}$; **0 is the correct Float64 answer** |
| Branch continuity at $x = 35$ | series ≤ 8.81e−16, closed form ≤ 3.28e−16 |
| Rank robustness at vanishing $\mu$ | finite and equal to Coulomb for $k = 0 \ldots 30$ (was `Inf` at $k{=}12$, `NaN` at $k{\ge}14$ before §12 problem 8) |

**Both exclusions are documented in the diagnostic's own output**, with counts, so that a future reader sees
what was excluded and why rather than a bare "0 violations". Neither is a widened tolerance: in both, Float64
cannot represent the quantity at all, and the unscreened Coulomb kernel fails identically at the same points.

`[NUMERICAL VERIFICATION]` Cache safety, from the same run:

```
XLCache, one object across a change of mu:
  mu = 0.0  (fresh)        X^0 = 4.208057097553249
  mu = 0.5  (same cache)   X^0 = 3.346223116076047      differs: true
  mu = 0.0  (same cache)   X^0 = 4.208057097553249      identical (===): true
vkCache, within one mu:
  mu = 0.0 :  uncached 2.104028548776624   cached 2.104028548776624   re-read ...   exact: true
  mu = 0.7 :  uncached 1.533638006593006   cached 1.533638006593006   re-read ...   exact: true
Tail beyond the source extent (mu = 0.8):
  power-law exponent -11.085 (Coulomb would be -1.000);  chord decay -0.879 vs -0.879 predicted
```

## 14.9 Grid convergence

`[NUMERICAL VERIFICATION]` Following CLAUDE.md Rule 12, a compact and a deliberately diffuse case:

| System | Box varied | Converged to |
|---|---|---|
| Be $1s^2 2s^2$ (compact) | $r_{\text{box}}$ 25 → 40 $a_0$ | **~1e−9 Ha** |
| Li $1s^2 3s$ (diffuse) | $r_{\text{box}}$ 70 → 110 $a_0$ | **~1e−7 Ha** |

`[INFERENCE]` The diffuse case is included deliberately: screening acts most strongly on the outermost, most
weakly bound electron, so a diffuse state is where a box-size artefact would most easily be mistaken for a
plasma effect. Rule 12 records that a too-small box corrupts *progressively* and therefore produces a
plausible-looking trend rather than an obvious failure.

## 14.10 Performance

See §11.5 for the full table and the measurement method. Summary: the screened potential sweep costs
**1.08–1.13×** the Coulomb sweep, and that ratio is **flat** over a 5.4× range in grid size — so no
$O(N_r^2)$ has been introduced. The end-to-end SCF cost is roughly **1.7×** (single runs, ±10 %).

## 14.11 The tested path IS the production path

`[NUMERICAL VERIFICATION]` Two independent instruments, neither of which changes any number.

**(a) A caller census.** `Defaults.eeScreeningMu()` was temporarily replaced by a version returning the
*identical* value and additionally recording its caller from `stacktrace()`. Over six production runs
(He / Be / Ne at two Debye lengths, through `Atomic.Computation` -> `Basics.perform`):

```
  XL_Coulomb                        148 calls      <- the CI Hamiltonian       (Path A)
  computePotential                  108 calls      <- the SCF direct potential (Path C)
  #buildScreenedPotential#1          88 calls      <- the Slater radial kernel (Path A)
  checkEeScreeningIsConsistent       12 calls      <- the consistency guard
```

Both paths of §3.2 are reached, from the ordinary user entry point, in a real calculation.

**(b) A sampling profile of one screened `Basics.perform` run** (Be, two configurations, mu = 0.2), which
proves the Yukawa code is not merely *dispatched to* but actually *executed*:

```
  computePotential               1264 samples        SlaterRkKinkAware             195
  yukawaKernel                   1202 samples        buildScreenedPotential        195
  iScaledYukawa                   977 samples        buildScreenedPotentialYukawa  189
  kScaledYukawa                    63 samples
```

`[CODE EVIDENCE]` Together with the source chain
`Basics.perform` (`module-BasicsAZ-inc-perform.jl:15`) -> `SelfConsistent.performSCF` (line 21 of the same
file) -> {`Basics.computePotential`; `Hamiltonian.performCI` -> `InteractionStrength.XL_Coulomb` ->
`RadialIntegrals.SlaterRkKinkAware` -> `buildScreenedPotential`}, the chain is closed at both ends: the source
says the call path exists, and the instrumentation says it was taken.

---

# 15. Verification Status

| Requirement | Status | Evidence |
|---|---|---|
| Yukawa interaction implemented | **ANALYTICALLY + NUMERICALLY VERIFIED** | derivation §2.3, §9.2; kernel vs an angular projection using no Bessel function, ~1e−15 (§14.2) |
| Normalisation ($i_0 = \sinh x/x$, $k_0 = e^{-x}/x$; **not** the $\pi/2$ convention) | **ANALYTICALLY + NUMERICALLY VERIFIED** | closed form $\sinh(\mu r_<)e^{-\mu r_>}/(\mu r_< r_>)$ to 2.06e−16 (§14.2) |
| Coulomb limit | **ANALYTICALLY + NUMERICALLY VERIFIED** | exact derivation §2.4; bit-identical at $\mu=0$; observed orders 1 ($k{=}0$) and 2 ($k{\ge}1$); production $\Delta E/(-6\mu) \to 1.000001$ (§14.3) |
| Direct terms ($F^k$) screened | **SOURCE-VERIFIED + NUMERICALLY VERIFIED** | one shared radial call in `XL_Coulomb` (§9.4); $F^0(1s,2s)$ moves (§14.6) |
| Exchange terms ($G^k$) screened | **SOURCE-VERIFIED + NUMERICALLY VERIFIED** | same call; He singlet–triplet splitting $= 2G^0$ moves 3.4558e−2 → 7.3395e−2 Ha (§14.6) |
| SCF consistency | **SOURCE-VERIFIED + NUMERICALLY VERIFIED** | `:meanFieldIteration` re-optimises (§10.1); full vs frozen differs by 9.05 % at $\mu = 1$ (§10.5); all five mean fields respond (§14.7) |
| The tested path is the production path | **NUMERICALLY VERIFIED** | caller census + sampling profile of a real `Basics.perform` run (§3.3, §14 below) |
| Cache safety | **SOURCE-VERIFIED + NUMERICALLY VERIFIED** | full audit §11.3; round trip across $\mu$ returns `===` the original (§14.8) |
| No Coulomb power-law tail under screening | **SOURCE-VERIFIED + NUMERICALLY VERIFIED** | separate tail branch (§7.3); measured exponential, −0.879 vs −0.879 predicted (§11.4) |
| Small-$x$ stability | **NUMERICALLY VERIFIED** | 0 NaN / 0 Inf down to $\mu r = 10^{-30}$; rank-aware floor verified $k = 0 \ldots 30$ (§14.8) |
| Large-$x$ stability | **NUMERICALLY VERIFIED** | scaled formulation, no clamp; verified to $\mu r = 10^{6}$; underflow to 0 beyond 745 Debye lengths is correct (§14.8) |
| Production calculation | **NUMERICALLY VERIFIED** | He / Be / Ne through `Atomic.Computation` → `perform` at two $\mu$ (§14.5) |
| Unscreened regression | **NUMERICALLY VERIFIED** | 83 / 83, no approved reference re-approved (§14.1); unscreened path is byte-identical by construction (§6.4) |
| Electron–nucleus potential untouched | **SOURCE-VERIFIED** | no edit outside the five e–e `computePotential` methods; at $\mu = 200$ the Be energy saturates at −19.8897 Ha, the bare hydrogenic four-electron sum |
| Invalid parameters refused | **NUMERICALLY VERIFIED** | $\lambda_D \in \{0, -1, \text{NaN}, \text{Inf}\}$ and a mis-ordered kernel call all raise (§14.2 checks, §6.5) |
| Plasma shift of a *transition* energy | **NOT VERIFIED** | see §17; only total energies were validated |
| Agreement with an external code's published Debye results | **NOT VERIFIED** | no cross-code comparison of final energies was performed; only the kernel convention and the tail hazard were cross-checked (§2.3, §11.4) |

---

# 16. What Was NOT Modified

`[CODE EVIDENCE]` Each of the following was examined and deliberately left alone.

| Item | Why not modified |
|---|---|
| **Angular coefficients** (`SpinAngular`, `AngularMomentum`) | The substitution replaces the coefficient of $P_k(\cos\omega)$ and leaves $P_k(\cos\omega)$ itself. Every Wigner coefficient, reduced matrix element and selection rule depends only on rank and quantum numbers, never on the radial kernel (§9.5). Modifying them would have been a *bug*. |
| **The electron–nucleus potential** | Explicitly out of scope. Nuclear Debye screening already exists separately as `Nuclear.nuclearPotentialDH` and is a different physical effect (§1.3). |
| **Breit / Gaunt** (`InteractionStrength.XL_Breit`) | The transverse photon-exchange part of the e–e interaction has its own kernels. A Debye-screened Coulomb term beside a bare Breit term is not a defined approximation, so the combination is **refused** rather than silently half-screened (§17). Breit screening was also explicitly out of scope. |
| **QED corrections** | Out of scope; untouched. |
| **Transition operators** (`PhotoEmission`, `MultipoleMoment`, …) | They are one-body electromagnetic operators; the plasma screening of the e–e interaction does not enter them. They see the screening only through the orbitals and energies, which is correct. |
| **`AaDFSField`, `AaHSField`** | The average-atom plasma line — a separate model with its own chemical potential, out of scope. A patch that touched `AaDFSField` was written, detected and reverted (§12, problem 3). |
| **`ThomasFermiField`** | A statistical *starting* potential built from $Z$ and the electron number, not a Hartree term from the density. Thomas–Fermi screening was explicitly out of scope. |
| **`RadialIntegrals.SlaterRk`** (the naive $O(N^2)$ form) | Used only by `module-AtomicFeatures.jl` to build feature vectors for the `DeepLearning` module. It takes no part in any Hamiltonian. Recorded in §17 rather than changed. |
| **`Plasma.LineShiftScheme`'s design** | It still screens only the CI matrix and not the SCF. That is its own documented approximation and changing it would be new physics, not this task. Its *numerics* did change — see §17, item 5. |
| **The `lambda` argument name of `SlaterRkDebyeHueckel`** | It is a public signature and has always in fact carried $\mu$. Renaming it is an unrelated API change; it is now documented instead (§6.3). |
| **Trailing whitespace flagged by `git diff --check`** | Two docstring lines. This is JAC **house style** — 107 such lines already exist in `module-Defaults.jl` and 113 in `module-RadialIntegrals.jl` — and the two trailing spaces are a Markdown hard line break, so removing them would change how the docstring renders and make the diff inconsistent with its own file. Classified **DO NOT CHANGE**. |
| **The `rtol` keyword** | Accepted but unused by `buildScreenedPotential` **before** this work; the Yukawa counterpart mirrors it for signature symmetry. Removing it is unrelated cleanup. |
| **Exporting `Defaults.eeScreeningMu`** | Adds one name to the `using` surface. Classified **OPTIONAL CLEANUP** and not done: it is a legitimate query, it is documented, and un-exporting functioning code for aesthetics is outside this task's remit. |

---

# 17. Limitations and Open Issues

## 17.1 Known implementation limitations

**(1) `ALField` and `EOLField` cannot be screened — and are REFUSED.**
`[CODE EVIDENCE]` They optimise their orbitals through the B-spline tensor line
(`buildScreenedPotentialPair` / `buildScreenedPotentialCache` → `XL_CoulombTensor`), a construction separate
from the tabulated-orbital sweep the switch reaches. Their SCF would run on the bare Coulomb interaction while
the Hamiltonian built from it was screened. `SelfConsistent.checkEeScreeningIsConsistent` raises a naming error
instead. **This is deliberate**: a partly screened calculation converges and returns plausible numbers with
nothing in the output to say half the interaction was bare. `[NUMERICAL VERIFICATION]` Both refusals fire, and
both fields are still accepted with screening off. Supported: `DFSField` (the default), `HSField`, `KSField`,
`CHField`.

**(2) Breit and Gaunt are refused in combination with screening.** `eeInteraction` or `eeInteractionCI` set to
`BreitInteraction`, `CoulombBreit` or `CoulombGaunt` raises. Same reasoning; same verification.

**(3) The screening parameter is a session global.** Two concurrent calculations in one Julia session cannot
use different Debye lengths, and a threaded CI matrix reads one value. This follows from the design decision of
§6.1 and is the price of guaranteeing that the SCF and the Hamiltonian agree.

**(4) `Basics.AtomicFeatures` is unscreened.** It computes $F^k$/$G^k$ with the naive
`RadialIntegrals.SlaterRk`, which writes the Coulomb kernel inline and is not reached by the switch. It builds
feature vectors for `DeepLearning` and takes no part in any Hamiltonian — but a feature vector computed during
a screened session is a **Coulomb** one. Documentation point, not an inconsistency in any energy.

**(5) `Plasma.LineShiftScheme` and `AutoIonization.computeLinesPlasma` NOW RETURN DIFFERENT NUMBERS.**
`[NUMERICAL VERIFICATION]` They reach `SlaterRkDebyeHueckel`, whose body was replaced. The old body truncated
the $i_k$ series at $p \le 2$ and was in error by 1.7e−4 at $\mu r_< = 1$, 7.4e−3 at 2, 0.15 at 4 and
**0.908 at $\mu r_< = 10$**. **The new numbers are the correct ones.** No approved reference is affected —
`test/approved/test-PlasmaShift-approved.sum` is referenced by nothing (`[CODE EVIDENCE]`, grep over `src/` and
`test/`), an orphan of a test deleted on 09-Aug-2026. Those two schemes are otherwise unchanged and in
particular still screen only the CI matrix, which is their own documented design.

## 17.2 Known physical-model limitations

**(6) $\mu$ is constant in $r$.** This is the Debye–Hückel model as specified. Ion-sphere, Stewart–Pyatt,
Ecker–Kröll, dynamic (frequency-dependent) and Thomas–Fermi screening were all explicitly out of scope and none
was attempted.

**(7) The electron–nucleus interaction is not screened here.** `[PHYSICS]` A physically complete Debye plasma
model screens both the e–e repulsion and the e–n attraction, and they shift the energy in *opposite*
directions. JAC has `Nuclear.nuclearPotentialDH` for the second; combining them is a modelling decision for the
user and was not made here.

**(8) Linear response, weak coupling.** `[PHYSICS]` Debye–Hückel is a linearised Poisson–Boltzmann result. At
strong coupling ($\lambda_D$ comparable to or smaller than the orbital radii) it is being used outside its own
derivation. The implementation remains *numerically* sound there — `[NUMERICAL VERIFICATION]` finite energies
with no NaN or Inf to $\mu = 200\,a_0^{-1}$ ($\lambda_D = 0.005\,a_0$), saturating at −19.8897 Ha, which is the
bare hydrogenic four-electron sum for Be — but the *model* is not, and the saturation is the physical statement
that the e–e interaction has been switched off entirely.

## 17.3 Known numerical limitations

**(9) Two Float64 edges, both documented and both shared with the unscreened code.** Below
$\mu(r_> - r_<) \approx -745$ the kernel underflows to an exact 0, which is correct. At radii below about
$10^{-24}\,a_0$ with rank $\ge 10$, $\hat\imath_k$ becomes subnormal and loses precision; below
$10^{-30}\,a_0$ the *unscreened* Coulomb kernel is itself `Inf`. Neither regime is reachable on any JAC grid
(`rnt` $\approx 2\times10^{-6}$).

**(10) Above rank 16 the small-$\mu$ Coulomb substitution becomes approximate.** `[ANALYTICAL VERIFICATION]`
`yukawaXFloor(k)` exceeds the Float64 epsilon for $k > 16$, so below that floor the substitution carries a
relative error of order the floor rather than being exact. It is a graceful degradation replacing a `NaN`, and
$k > 16$ requires $l > 8$, which no atomic Slater integral in practice reaches.

## 17.4 Untested features and unsupported modes

| Item | Status |
|---|---|
| Plasma shift of a **transition** or excitation energy | **NOT VERIFIED.** Only total energies were validated. Differences cancel the dominant $-\mu N(N-1)/2$ term and are correspondingly more sensitive to grid and basis (§14.5). |
| Comparison against an external code's published Debye energies | **NOT VERIFIED.** The kernel convention and the tail hazard were cross-checked against an independent implementation; final energies were not. |
| Screening in `Cascade`, `PhotoIonization`, `DielectronicRecombination`, … | **NOT VERIFIED.** These reach the same `performSCF` and the same Slater machinery, so `[INFERENCE]` they should inherit the screening — but no such calculation was run screened. |
| Continuum orbitals under screening | **NOT VERIFIED.** Not exercised. |
| Screening with `scField = NuclearField()` or a one-electron system | Degenerate by construction: there is no e–e interaction to screen. |
| Thread safety of the global under `Threads.@threads` CI assembly | `[INFERENCE]` reads only, never written during a calculation, so it should be safe; **NOT VERIFIED** by a threaded run. |

## 17.5 Future work, recorded and deliberately not done

* a Yukawa counterpart of the B-spline tensor line, which would lift limitation (1);
* screening the Breit kernels, which would lift limitation (2);
* a `Plasma.Computation` scheme that screens the SCF as well as the CI matrix, now that the machinery exists;
* a cross-code validation of final energies against published Debye–Hückel results;
* deleting the orphaned `test/approved/test-PlasmaShift-approved.sum` — an editorial act for the maintainer.

---

# 18. Reproducible Workflow

## 18.1 Environment

```
repository   JenaAtomicCalculator.jl
base commit  75f963acac367e7365eed27854199a572a465d2c   (branch master)
Julia        1.12.7
platform     Linux 6.8.0-138-generic, x86_64
```

```bash
git clone <the JAC repository>
cd JenaAtomicCalculator.jl
git checkout 75f963acac367e7365eed27854199a572a465d2c
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

`Manifest.toml` is machine- and Julia-specific and must **not** be copied between machines; if the environment
misbehaves, `rm Manifest.toml` and instantiate again.

## 18.2 Running the test suite

```bash
cd <JAC-root>/test
julia --project=.. runtests.jl
```

**This form, and not `Pkg.test()` or a run from the repository root.** `[CODE EVIDENCE]` CLAUDE.md Rule 11
requires it: every `TestFrames` include opens its summary with a *relative* path, so the working directory
decides where `test-*.sum` and `zzz-*.sum` land. Running from the root scatters them across the repository.
The Rule 11 form runs the identical `runtests.jl` in the identical project.

Expected: `Test Summary: | Pass Total` → `83 83`, about 7m45 s.

## 18.3 The screening API

```julia
using JenaAtomicCalculator

Defaults.setDefaults("e-e screening", Basics.DebyeHueckelModel(5.0))   # lambda_D = 5 a_o
Defaults.setDefaults("e-e screening", Basics.NoPlasmaModel())          # off (the default)

Defaults.getDefaults("e-e screening")   # the model in force
Defaults.eeScreeningMu()                # mu = 1/lambda_D, exactly 0.0 when unscreened
```

## 18.4 A complete example calculation

```julia
using JenaAtomicCalculator

grid    = Radial.Grid(Radial.Grid(false); rbox = 25.0, printout = false)
configs = [Configuration("1s^2 2s^2")]

function groundEnergy(lambdaD)                     # lambdaD = nothing  ->  unscreened
    Defaults.setDefaults("e-e screening",
        isnothing(lambdaD) ? Basics.NoPlasmaModel() : Basics.DebyeHueckelModel(lambdaD))
    try
        comp = Atomic.Computation(Atomic.Computation(), name = "Be", grid = grid,
                                  nuclearModel = Nuclear.Model(4.), configs = configs)
        mp = Basics.perform(comp; output = true)["multiplet:"]
        return  minimum(lv.energy for lv in mp.levels)
    finally
        Defaults.setDefaults("e-e screening", Basics.NoPlasmaModel())   # ALWAYS restore
    end
end

e0 = groundEnergy(nothing)      # -14.570977349  Ha
e1 = groundEnergy(5.0)          # -15.533391190  Ha   (mu = 0.2 a_o^-1)
println("Delta E = ", e1 - e0, " Ha;   -mu N(N-1)/2 = ", -0.2*4*3/2)
```

**The `try`/`finally` is not optional.** The switch is a session global; leaving it set would silently screen
every later calculation in the same session.

## 18.5 Validation commands

All five diagnostics live in `tools/` and are run from the repository root. None of them is required to use the
feature; they exist so that every number in this report can be regenerated.

```bash
julia --project=. tools/diag-debye-screening-kernel.jl        # Levels 1-4: Bessel, kernel, angular projection,
                                                              #   the k=0 closed form, O(N) vs O(N^2)
julia --project=. tools/diag-debye-screening-atomic.jl        # Levels 5-7: exchange, SCF consistency, the
                                                              #   -mu N(N-1)/2 law, grid convergence, guards
julia --project=. tools/diag-debye-screening-production.jl    # the production path: perform() runs, the caller
                                                              #   census, the profile, the Coulomb limit
julia --project=. tools/diag-debye-screening-stability.jl     # 256-bit reference, 13544-point sweep, cache
                                                              #   audit, tail model, timings
julia --project=. tools/diag-debye-screening-performance.jl   # the retired truncated kernel, isolated and timed
```

To reproduce the regression comparison specifically:

```bash
git stash                                   # set the working tree back to the base commit
cd test && julia --project=.. runtests.jl   # expect 82 / 82
cd .. && git stash pop
cd test && julia --project=.. runtests.jl   # expect 83 / 83
```

---

# 19. File-by-File Change Summary

## 19.1 Production files — these belong in the commit

| File | Changes | Reason | Production? | Test coverage |
|---|---|---|---|---|
| `src/module-Defaults.jl` | +74 / −1 | the screening state, its validating setter, its getter, and `eeScreeningMu()` | **yes** | suite checks 5, 8 (round trip, four refusals) |
| `src/module-RadialIntegrals.jl` | +351 / −35* | Bessel primitives, kernel, rank-aware floor, $O(N)$ Yukawa sweep, `SlaterRkYukawa`, the **one branch point**, `Yk_ab(…,mu)`, `SlaterRkDebyeHueckel` body | **yes** | suite checks 1–4b, 5, 7; all five diagnostics |
| `src/module-BasicsAZ-inc-compute.jl` | +65 / −6 | the five e–e mean-field potentials (Path C) | **yes** | suite checks 10, 11; §14.7 |
| `src/module-InteractionStrength.jl` | +12 / −2 | $\mu$ in both `XLCache` keys | **yes** | §11.1, §14.8 |
| `src/module-SelfConsistent.jl` | +58 | `checkEeScreeningIsConsistent` + the screened-run notice | **yes** | suite check 8; guard diagnostic |
| `src/module-TestFrames-inc-plasma.jl` | +285 | `testModule_PlasmaScreening` — 11 checks, no approved data | **yes** (test suite) | is the coverage |
| `test/runtests.jl` | +4 | wires it into the previously empty `@testset "JAC plasma"` | **yes** (test suite) | — |

Totals, from `git diff --numstat` at the close of this work: **7 files, 849 insertions, 44 deletions.**

\* 35 of the 44 deletions are the retired truncated `ul_DH` body of `SlaterRkDebyeHueckel` (§5.1c); the other
nine are the six `max()` lines replaced by branches in `computePotential`, two cache-key lines, and one
`export` line. **No line was deleted that this work did not replace.**

## 19.2 Documentation — belongs in the commit

| File | What |
|---|---|
| `docs/debye-electron-electron-screening.md` | **this report** — the permanent technical record |
| `progress.md` | the chronological development log (currently untracked; see §20) |

## 19.3 Diagnostic files — `tools/`, and a maintainer decision

`[CODE EVIDENCE]` CLAUDE.md Rule 6 designates `tools/` as "Assistant-maintained diagnostics and cross-code
comparison harnesses", **tracked**, and states that the maintainer is not expected to read or maintain it. So
these are eligible for the commit by the repository's own rules, and the reason to keep them is that they are
what regenerates every number in §14.

| File | What it validates | Runtime |
|---|---|---|
| `tools/diag-debye-screening-kernel.jl` | Levels 1–4: Bessel vs `SpecialFunctions`; the angular Legendre projection; the $k=0$ closed form; $O(N)$ vs $O(N^2)$ | seconds |
| `tools/diag-debye-screening-atomic.jl` | Levels 5–7: exchange, SCF consistency, $-\mu N(N-1)/2$, grid convergence, stronger screening, 7 guards | minutes |
| `tools/diag-debye-screening-production.jl` | the production path: `perform()` runs, the **caller census**, the **profile**, the Coulomb limit, SCF-vs-frozen | minutes |
| `tools/diag-debye-screening-stability.jl` | the 256-bit reference, the 13 544-point sweep, the cache audit, the tail model, the timings | minutes |
| `tools/diag-debye-screening-performance.jl` | reproduces the **retired truncated kernel** verbatim and quantifies what changed | seconds |

**Files that must NOT be committed.** `[CODE EVIDENCE]` These are present in the working tree and are not part
of this work:

| File | Why not |
|---|---|
| `CLAUDE.local.md` | the user's private project instructions, deliberately untracked |
| `OVERNIGHT_TASK.md`, `progress_before_resume.md` | session scaffolding, not repository content |
| `launch_claude.py`, `automation_logs/` | local automation, unrelated to JAC |

These were present before this work began and have not been touched.

---

# 20. Final Git Diff Review

Commands and their output at the close of this work are recorded in `progress.md`; the summary is:

```bash
git diff --stat     # 7 files changed
git diff --check    # 2 trailing-whitespace hits, both JAC house style -- see below
git status          # 7 modified, plus untracked docs/diagnostics/scaffolding
```

**What the working tree contains:**

| Category | Contents | Should be committed? |
|---|---|---|
| Intended implementation changes | the 5 `src/` modules of §19.1 | **yes** |
| Test changes | `src/module-TestFrames-inc-plasma.jl`, `test/runtests.jl` | **yes** |
| Documentation | `docs/debye-electron-electron-screening.md`, `progress.md` | **yes** |
| Diagnostic scripts | the 5 `tools/diag-debye-screening-*.jl` | **yes** — Rule 6 tracks `tools/` |
| Unrelated files | `CLAUDE.local.md`, `OVERNIGHT_TASK.md`, `progress_before_resume.md`, `launch_claude.py`, `automation_logs/` | **no** — pre-existing, not this work's |

**On `git diff --check`.** It reports trailing whitespace on two lines:
`src/module-Defaults.jl:696` and `src/module-RadialIntegrals.jl:1420`, both of the form
`` `Module.function(args)`  `` at the head of a docstring. `[CODE EVIDENCE]` This is **JAC house style**: 107
such lines already exist in `module-Defaults.jl` and 113 in `module-RadialIntegrals.jl`, and the two trailing
spaces are a Markdown hard line break, so removing them would change how the docstring renders and make the new
code inconsistent with the file around it. Classified **DO NOT CHANGE**.

**No commit has been made.** CLAUDE.md requires explicit user approval before `git add`/`git commit`, and
separately before any push.

---

# 21. Final Assessment

```
Implementation status:   COMPLETE
Validation status:       VALIDATED at seven layers, against two mathematically independent references,
                         with the production path proved by instrumentation of a real run.
                         Four items remain NOT VERIFIED and are named in §17.4 -- chief among them the
                         plasma shift of a TRANSITION energy, and any cross-code comparison of final
                         energies.
Regression status:       CLEAN.  83 / 83 (82 baseline + 1 added), no approved reference re-approved,
                         and the unscreened path is byte-identical by construction rather than by
                         tolerance.
Production readiness:    READY for the calculation modes it supports -- scField in {DFSField (default),
                         HSField, KSField, CHField} with a Coulomb e-e interaction, through
                         Atomic.Computation / perform.  NOT ready, and explicitly REFUSES, for
                         ALField, EOLField, and any Breit or Gaunt interaction.
Known limitations:       10 items in §17, of which the three that matter for a user are:
                         (1) ALField/EOLField are refused;
                         (2) Breit/Gaunt are refused;
                         (3) the parameter is a session global.
                         Plus one behavioural change to flag: Plasma.LineShiftScheme and
                         AutoIonization.computeLinesPlasma now return DIFFERENT -- and correct --
                         numbers, because the truncated kernel they used was replaced (§17.1 item 5).
Recommended next action: MAINTAINER REVIEW OF THE DIFF, then a commit of the seven source/test files
                         plus this report, plus progress.md and the five tools/ diagnostics.
                         Do NOT commit CLAUDE.local.md, OVERNIGHT_TASK.md, progress_before_resume.md,
                         launch_claude.py or automation_logs/.
                         Two decisions are the maintainer's alone and are not taken here:
                           - whether the orphaned test/approved/test-PlasmaShift-approved.sum should be
                             deleted, now that nothing references it;
                           - whether the changed Plasma.LineShiftScheme numbers warrant a note in
                             docs/src/news.md before the next release.
```

## A closing note on what this report is for

`[INFERENCE]` The most valuable single fact in this document is not any of the measured numbers. It is that
JAC contains **two independent implementations of the electron–electron Coulomb interaction** — the Slater
machinery and the mean-field potential — and that they are reached by different call chains, share no code, and
must be changed together. That was discovered by tracing the source, not by reading documentation, and an
implementation that missed it would have converged, returned plausible energies, and been wrong by up to 9 %
with nothing in the output to say so.

The second most valuable is §12: three of the ten problems recorded there were failures of a *reference* or a
*test criterion* rather than of the implementation, and in each the tempting repair was to widen a tolerance.
Anyone reimplementing this should expect the same, and should treat a disagreement between an implementation
and its reference as an open question with two suspects.
