
## Functions in this file cover: plasma shift and plasma environment effects.
## Alphabetical order within this file.


## PlasmaShift: no test exists and none is planned here. The former testModule_PlasmaShift was deleted
## on 09-Aug-2026 as obsolete: PlasmaShift is no longer an Atomic property but a Plasma.Computation
## scheme, so a test for it belongs with the Plasma module and has to be written afresh.


"""
`TestFrames.testModule_PlasmaScreening(; short::Bool=true)`
    ... tests the Debye-Hueckel (Yukawa) screening of the ELECTRON-ELECTRON interaction, which is switched on by
        Defaults.setDefaults("e-e screening", Basics.DebyeHueckelModel(lambda_D)) and replaces 1/r_12 by
        exp(-r_12/lambda_D)/r_12 in the Slater integrals of the Hamiltonian and in the direct potential of the
        mean-field SCF alike.  A success::Bool is returned.

        NO APPROVED DATA IS USED and none can be: every check below is either an exact closed form, an identity
        between two independent constructions, or the mu -> 0 limit, so this test cannot pass on a stale stored
        reference.  The one thing it does NOT check is a value against an independent code; that belongs with
        the diagnostics in tools/, which compare the kernel against an angular Legendre projection.

        The checks, in the order they are made:
          (1) the NORMALISATION of the two Bessel primitives against i_0(x) = sinh(x)/x and k_0(x) = exp(-x)/x,
              on BOTH branches of iScaledYukawa.  This is the check worth having first: the spherical k_n of GSL
              and of SpecialFunctions are (pi/2) times the k_k used here, and a code that adopted the library
              convention would still converge, still look plausible, and be wrong by 2/pi in every screened
              matrix element.
          (2) the rank-0 kernel against its closed form sinh(mu r_<) exp(-mu r_>) / (mu r_< r_>), which involves
              no Bessel function at all;
          (3) the COULOMB LIMIT: the kernel at mu = 0 is bit-identical to r_<^k/r_>^(k+1), and the relative
              error at small mu falls with the documented order -- 1 for k = 0, where the leading correction is
              the constant -mu, and 2 for k >= 1;
          (4) the SIGN: screening only ever weakens the repulsion, at every rank and radius;
          (5) that the screening is OFF by default and that switching it on and off again restores the
              production Slater integral BIT FOR BIT, not merely to within a tolerance;
          (6) that the DIRECT and the EXCHANGE slot assignment both respond -- one shared kernel serves both,
              and this is what says so;
          (7) that the O(N) forward/backward sweep of buildScreenedPotentialYukawa reproduces the definition,
              by contracting the same orbital pair against the kernel directly;
          (8) that an invalid Debye length and a wrongly ordered kernel call are REFUSED rather than absorbed;
          (9) a PRODUCTION regression: a real SelfConsistent.performSCF of helium, screened and unscreened,
              against the first-order shift  Delta E = -mu N(N-1)/2, which is exact as mu -> 0 and independent
              of the wavefunction because the operator is then a constant;
         (10) that ALL FIVE screened mean fields respond, and respond in the right DIRECTION -- a production
              run exercises only the one its AsfSettings names, so without this the other four would be
              covered by inspection alone.
"""
function testModule_PlasmaScreening(; short::Bool=true)
    Defaults.setDefaults("print summary: open", "test-PlasmaScreening-new.sum")
    printstyled("\n\nTest the module  Plasma:  Debye-Hueckel screening of the e-e interaction  ... \n", color=:cyan)
    success = true;    printTest, iostream = Defaults.getDefaults("test flag/stream")
    # The switch is a module global, so a failure part-way through this test would leave every LATER test in the
    # suite running under screening. Everything below therefore sits inside a try/finally that restores it.
    setMu(mu) = redirect_stdout(devnull) do
        Defaults.setDefaults("e-e screening", mu == 0.  ?  Basics.NoPlasmaModel()  :  Basics.DebyeHueckelModel(1/mu))
    end

    try
        # (1) normalisation of the Bessel primitives; x = 50 exercises the closed-form branch of iScaledYukawa,
        #     the others its ascending series. exp(-x) i_0(x) = (1 - exp(-2x))/(2x) is written in the scaled form
        #     on purpose: sinh(50) alone is 2.6e21 and would compare two overflow-prone numbers. It is then
        #     written with expm1 rather than as 1 - exp(-2x), which at x = 1e-6 subtracts two numbers agreeing in
        #     their first eleven digits and leaves the REFERENCE wrong by 5e-12 -- comfortably enough to fail a
        #     correct implementation, and exactly the way a bad tolerance gets adopted to hide it.
        for  x  in  [1.0e-6, 1.0e-2, 0.5, 3.0, 20.0, 50.0, 400.0]
            wa = RadialIntegrals.iScaledYukawa(0, x);    wb = -expm1(-2x) / (2x)
            wc = RadialIntegrals.kScaledYukawa(0, x);    wd = 1.0 / x
            if  abs(wa - wb) > 1.0e-12 * abs(wb)  ||  abs(wc - wd) > 1.0e-12 * abs(wd)
                success = false
                if printTest   info(iostream, "Bessel normalisation at x = $x: ihat_0 = $wa against $wb, " *
                                              "khat_0 = $wc against $wd")   end
            end
        end
        # (2) the rank-0 kernel against its closed form. This is the same statement as (1) seen from the kernel,
        #     and it is the one that pins the (pi/2): a kernel built from the library convention fails here by
        #     exactly 2/pi.
        for  (rs, rl, mu)  in  [(0.1, 0.5, 1.0), (0.5, 2.0, 0.3), (1.0, 1.0, 2.0), (0.02, 30.0, 0.5),
                               (2.0, 4.0, 1.0e-3)]
            wa = RadialIntegrals.yukawaKernel(0, rs, rl, mu)
            wb = sinh(mu*rs) * exp(-mu*rl) / (mu * rs * rl)
            if  abs(wa - wb) > 1.0e-11 * abs(wb)
                success = false
                if printTest   info(iostream, "U_0($rs, $rl; mu = $mu) = $wa, closed form $wb")   end
            end
        end
        # (3) the Coulomb limit. At mu = 0 the kernel must be the Coulomb one EXACTLY -- an unscreened run has to
        #     take the old code path, not a screened one that agrees with it -- and the approach to it must carry
        #     the documented order. The k = 0 correction is the constant -mu, so its RELATIVE error is O(mu r_>)
        #     and the observed order is 1; for k >= 1 the linear term cancels and the order is 2. Getting this
        #     backwards is easy and would hide a genuine first-order error at higher rank.
        for  k  in  0:4
            wa = RadialIntegrals.yukawaKernel(k, 0.7, 2.5, 0.0)
            if  wa != 0.7^k / 2.5^(k+1)
                success = false
                if printTest   info(iostream, "U_$k at mu = 0 is $wa, must be exactly $(0.7^k/2.5^(k+1))")   end
            end
            exact = 0.7^k / 2.5^(k+1)
            e1    = abs(RadialIntegrals.yukawaKernel(k, 0.7, 2.5, 1.0e-4) - exact) / exact
            e2    = abs(RadialIntegrals.yukawaKernel(k, 0.7, 2.5, 1.0e-5) - exact) / exact
            order = log10(e1 / e2)
            if  abs(order - (k == 0  ?  1.0  :  2.0)) > 0.05
                success = false
                if printTest   info(iostream, "the Coulomb limit of U_$k is approached with observed order " *
                                              "$order, expected $(k == 0 ? 1 : 2)")   end
            end
        end
        # (4) the sign. Screening removes interaction; it never adds any.
        for  k  in  0:4,  mu  in  [0.05, 0.5, 5.0],  (rs, rl)  in  [(0.05, 0.05), (0.3, 1.2), (1.0, 40.0)]
            wa = RadialIntegrals.yukawaKernel(k, rs, rl, mu)
            if  !(0. < wa < rs^k / rl^(k+1))
                success = false
                if printTest   info(iostream, "U_$k($rs, $rl; mu = $mu) = $wa is not in (0, " *
                                              "$(rs^k/rl^(k+1))), the Coulomb value")   end
            end
        end
        # (4b) HIGH RANK at vanishing mu, which is where khat_k is formed closest to its overflow. Until
        #      05-Sep-2026 the floor below which the Coulomb kernel is substituted was a single constant,
        #      1e-25, and that is safe only to k = 10: the kernel returned Inf at k = 12 and NaN at k >= 14,
        #      where the correct answer is the Coulomb value. The floor is now rank-aware, and the ranks
        #      beyond 10 are checked here BECAUSE the earlier sweep stopped at 8 and so never saw it.
        for  k  in  [0, 8, 10, 11, 12, 14, 16, 20],  mu  in  [1.0e-24, 1.0e-18, 1.0e-8, 1.0e-2]
            wa = RadialIntegrals.yukawaKernel(k, 1.0, 1.0, mu)
            if  !isfinite(wa)  ||  !(0. < wa <= 1.0 + 1.0e-12)
                success = false
                if printTest   info(iostream, "U_$k(1, 1; mu = $mu) = $wa, which must be finite and in " *
                                              "(0, 1]; the Coulomb value is 1")   end
            end
        end

        grid = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
        nm   = Nuclear.Model(10.)
        orb  = Dict{String,Orbital}()
        for  s  in  ["1s_1/2", "2s_1/2", "2p_1/2", "2p_3/2"]    orb[s] = HydrogenicIon.orbital(Subshell(s), nm, grid)    end

        # (5) the switch is off by default, and a round trip restores the production integral BIT FOR BIT. The
        #     comparison is `!=` and not a tolerance: the whole design rests on an unscreened run re-entering the
        #     original code path rather than a screened path that happens to agree with it.
        if  Defaults.eeScreeningMu() != 0.
            success = false
            if printTest   info(iostream, "the e-e screening is not off at the start of this test")   end
        end
        wBefore = RadialIntegrals.SlaterRkKinkAware(0, orb["1s_1/2"], orb["2s_1/2"], orb["1s_1/2"], orb["2s_1/2"], grid)
        setMu(0.5)
        wOn     = RadialIntegrals.SlaterRkKinkAware(0, orb["1s_1/2"], orb["2s_1/2"], orb["1s_1/2"], orb["2s_1/2"], grid)
        setMu(0.0)
        wAfter  = RadialIntegrals.SlaterRkKinkAware(0, orb["1s_1/2"], orb["2s_1/2"], orb["1s_1/2"], orb["2s_1/2"], grid)
        if  wBefore != wAfter
            success = false
            if printTest   info(iostream, "R^0(1s,2s,1s,2s) before screening $wBefore, after switching it off " *
                                          "again $wAfter; these must be identical")   end
        end
        if  !(0. < wOn < wBefore)
            success = false
            if printTest   info(iostream, "R^0(1s,2s,1s,2s) is $wOn under screening and $wBefore without; the " *
                                          "screened value must be smaller and positive")   end
        end

        # (6) DIRECT and EXCHANGE both respond. The two differ only in the slot assignment the angular
        #     coefficient supplies -- (a,b,a,b) is the direct F^k, (a,b,b,a) the exchange G^k -- and both are
        #     built from the one kernel, so a change that reached only one of them would be a defect this is the
        #     cheapest way to see.
        for  (label, L, sa, sb, sc, sd)  in  [("direct   F^0(1s,2s)", 0, "1s_1/2","2s_1/2","1s_1/2","2s_1/2"),
                                              ("exchange G^0(1s,2s)", 0, "1s_1/2","2s_1/2","2s_1/2","1s_1/2"),
                                              ("exchange G^1(1s,2p)", 1, "1s_1/2","2p_3/2","2p_3/2","1s_1/2")]
            setMu(0.0);    w0 = InteractionStrength.XL_CoulombKinkAware(L, orb[sa], orb[sb], orb[sc], orb[sd], grid)
            setMu(0.5);    w1 = InteractionStrength.XL_CoulombKinkAware(L, orb[sa], orb[sb], orb[sc], orb[sd], grid)
            setMu(0.0)
            if  abs(w0) < 1.0e-10  ||  abs(w1 - w0) < 1.0e-6 * abs(w0)
                success = false
                if printTest   info(iostream, "$label did not respond to the screening: $w0 -> $w1")   end
            end
        end

        # (7) the O(N) sweep against the definition. buildScreenedPotentialYukawa accumulates two running moments
        #     in exponentially scaled form; the double sum below contracts the very same orbital pair against
        #     RadialIntegrals.yukawaKernel with no sweep at all, so the two share nothing but the kernel. They
        #     are two quadratures of one integral and agree only to quadrature accuracy, which is what the 1e-3
        #     allows for; a wrong sweep is out by tens of percent, not by 0.1 %.
        for  (mu, k, sb, sd)  in  [(0.5, 0, "1s_1/2", "1s_1/2"), (1.0, 0, "2s_1/2", "2s_1/2"),
                                   (0.5, 1, "1s_1/2", "2p_3/2")]
            b = orb[sb];    d = orb[sd];    mtp = min(size(b.P, 1), size(d.P, 1))
            Vk = RadialIntegrals.buildScreenedPotentialYukawa(k, b, d, grid, mu; mtpOut = mtp)
            for  i  in  [max(2, mtp ÷ 8), mtp ÷ 3, mtp ÷ 2]
                r  = grid.r[i];    wa = 0.
                for  j = 2:mtp
                    rl = min(r, grid.r[j]);   rg = max(r, grid.r[j])
                    wa = wa + (b.P[j]*d.P[j] + b.Q[j]*d.Q[j]) * RadialIntegrals.yukawaKernel(k, rl, rg, mu) * grid.wr[j]
                end
                if  abs(Vk[i] - wa) > 1.0e-3 * max(abs(wa), 1.0e-30)
                    success = false
                    if printTest   info(iostream, "V^$k(r = $r) for ($sb,$sd) at mu = $mu: sweep $(Vk[i]), " *
                                                  "direct contraction $wa")   end
                end
            end
        end

        # (8) refusals. An invalid Debye length must be rejected where it is set, not carried into the kernels as
        #     an Inf or a NaN, and the kernel is not symmetric in its two radii so a wrongly ordered call is a
        #     caller error rather than something to repair silently.
        for  lambda  in  [0.0, -1.0, NaN, Inf]
            try
                redirect_stdout(devnull) do
                    Defaults.setDefaults("e-e screening", Basics.DebyeHueckelModel(lambda))
                end
                success = false
                if printTest   info(iostream, "lambda_D = $lambda was accepted; it must be refused")   end
            catch    # the expected outcome
            end
        end
        setMu(0.0)
        try
            RadialIntegrals.yukawaKernel(0, 2.0, 1.0, 0.5)
            success = false
            if printTest   info(iostream, "yukawaKernel accepted rSmall > rLarge; it must refuse")   end
        catch    # the expected outcome
        end

        # (9) THE PRODUCTION REGRESSION. A real self-consistent calculation of helium, run through the same
        #     SelfConsistent.performSCF that Basics.perform(::Atomic.Computation) calls, screened and unscreened.
        #     exp(-mu r_12)/r_12 = 1/r_12 - mu + O(mu^2 r_12), so at first order every electron PAIR contributes
        #     -mu and Delta E = -mu N(N-1)/2 exactly, whatever the orbitals are; the orbital relaxation is second
        #     order. mu = 1e-4 is small enough for the law to hold to ~1e-4 relative and large enough that the
        #     shift, 1e-4 Ha, sits far above the SCF convergence noise.
        heGrid  = Radial.Grid(Radial.Grid(false), rnt = 2.0e-6, h = 5.0e-2, hp = 2.0e-2, rbox = 20.0)
        heConf  = [Configuration("1s^2")];    muHe = 1.0e-4;    nHe = 2
        setMu(0.0);    eCb = minimum(lv.energy for lv in
                            SelfConsistent.performSCF(heConf, Nuclear.Model(2.), heGrid, AsfSettings(); printout=false).levels)
        setMu(muHe);   eDh = minimum(lv.energy for lv in
                            SelfConsistent.performSCF(heConf, Nuclear.Model(2.), heGrid, AsfSettings(); printout=false).levels)
        setMu(0.0)
        predicted = -muHe * nHe * (nHe - 1) / 2
        if  eDh >= eCb  ||  abs((eDh - eCb) - predicted) > 1.0e-3 * abs(predicted)
            success = false
            if printTest   info(iostream, "He 1s^2 at mu = $muHe: E = $eDh against $eCb unscreened, i.e. " *
                                          "Delta E = $(eDh - eCb) where -mu N(N-1)/2 = $predicted")   end
        end

        # (10) ALL FIVE screened mean fields, not only the default one. Basics.computePotential carries the
        #      screening in five methods, and a production run exercises whichever ONE the AsfSettings names --
        #      so DFSField alone was being tested and the other four only by inspection. The check is the
        #      physical one and needs no reference: these potentials are stored as Z(r) = -r V(r) with V the
        #      attractive total, so removing e-e repulsion by screening makes the electronic part smaller and
        #      Z(r) LARGER at every radius. HSField is the one that matters most to include, because it alone
        #      reaches the screening through RadialIntegrals.Yk_ab rather than through yukawaKernel.
        heLevel = redirect_stdout(devnull) do
                      SelfConsistent.performSCF(heConf, Nuclear.Model(2.), heGrid, AsfSettings();
                                                printout=false).levels[1]
                  end
        for  scF  in  [Basics.DFSField(), Basics.HSField(), Basics.KSField(), Basics.CHField()]
            setMu(0.0);    p0 = Basics.computePotential(scF, heGrid, heLevel)
            setMu(0.5);    p1 = Basics.computePotential(scF, heGrid, heLevel)
            setMu(0.0)
            nMoved = 0;    nWrong = 0
            for  i = 2:min(length(p0.Zr), length(p1.Zr))
                d = p1.Zr[i] - p0.Zr[i]
                if  d < -1.0e-10                        nWrong = nWrong + 1    end
                if  abs(d) > 1.0e-6 * max(abs(p0.Zr[i]), 1.0e-3)   nMoved = nMoved + 1    end
            end
            if  nWrong > 0  ||  nMoved == 0
                success = false
                if printTest   info(iostream, "$(nameof(typeof(scF))): screening moved Z(r) at $nMoved points " *
                                              "and moved it the WRONG way at $nWrong")   end
            end
        end
        # (11) the DFSField(basis) method, which is a separate method from DFSField(level) and is the one the
        #      mean-field SCF iterates against.
        setMu(0.0);    q0 = Basics.computePotential(Basics.DFSField(), heGrid, heLevel.basis)
        setMu(0.5);    q1 = Basics.computePotential(Basics.DFSField(), heGrid, heLevel.basis)
        setMu(0.0)
        if  !any(q1.Zr[i] - q0.Zr[i] > 1.0e-6 for i = 2:min(length(q0.Zr), length(q1.Zr)))
            success = false
            if printTest   info(iostream, "computePotential(DFSField, grid, basis) did not respond to the " *
                                          "screening")   end
        end

    finally
        redirect_stdout(devnull) do
            Defaults.setDefaults("e-e screening", Basics.NoPlasmaModel())
        end
    end

    println(iostream, "Plasma screening: the Debye-Hueckel e-e kernel against the closed forms i_0 = sinh(x)/x, " *
                      "k_0 = exp(-x)/x and U_0 = sinh(mu r_<) exp(-mu r_>)/(mu r_< r_>); its Coulomb limit and "  *
                      "the order in which it is approached; the sign; a bit-identical round trip of the switch; " *
                      "the response of the direct AND the exchange integral; the O(N) sweep against a direct "    *
                      "contraction of the same kernel; the finiteness of ranks 0..20 at vanishing mu; four "    *
                      "refusals; a self-consistent helium calculation against Delta E = -mu N(N-1)/2; and the "   *
                      "response, with the right sign, of all five screened mean fields. No approved data is used.")
    Defaults.setDefaults("print summary: close", "")
    testPrint("testModule_PlasmaScreening()::", success)
    return( success )
end
