MODULE reactions_module

! RATE SETS: Prefers most recent first, then continues
! down the list when #defines are set to zero.

! Use rates from 1992 code:
!   Smith, Kawano, Malaney, Astrophys. J. Supp. Ser. 85, 219-247 (1993)
#define rates_1992 1
! Use rates from NACRE database:
!   C. Angulo et al., Nucl. Phys. A 656, 3-187 (1999)
#define rates_1999 1
! Use rates of Cyburt review:
!   R. H. Cyburt, Phys. Rev. D 70, 023505 (2004)
#define rates_2004C 1
! Use rates of Descouvemont R-matrix analysis:
!   P. Descouvemont, A. Adahchour, C. Angulo, A. Coc, E. Vangioni-Flam,
!   At. Data. Nucl. Data Tables 88, 203 (2004)
#define rates_2004D 0

! Rate for p(n,g)d
!   Ando et al. Phys. Rev. C 74, 025809 (2006)
#define rate12_2006 1

! Rate for He3(a,g)Be7
!   Cyburt and Davids, arXiv:0809.3240 [nucl-ex]
! Refit from S(E) using mathematica
#define rate27_2008 0
#define rate27_2008_refit 1

! Integrate the cross-section explicitly for each value of T
#define integrate_rates 0
#if integrate_rates
    use rate_integrator
!  1    p(n,g)d
#define integrate_12 0
!  2    d(p,g)3He
#define integrate_20 0
!  3    d(d,n)3He
#define integrate_28 0
!  4    d(d,p)t
#define integrate_29 0
!  5  3He(n,p)t
#define integrate_16 0
!  6    t(d,n)4He
#define integrate_30 1
!  7  3He(d,p)4He
#define integrate_31 1
!  8  3He(a,g)7Be
#define integrate_27 0
!  9    t(a,g)7Li
#define integrate_26 0
! 10  7Be(n,p)7Li
#define integrate_17 0
! 11  7Li(p,a)4He
#define integrate_24 0
#endif

#define linear_q_variation 0
#define include_gamow_shift 1

    use physics_parameters
    use computational_parameters
    use variables

    REAL, SAVE:: lnvar

    REAL, PARAMETER :: T9LIMIT = 100.0  !Default range of reaction rate fits
    REAL, SAVE, DIMENSION(nrec) :: &
            t9max           !Maximum temp for which the reaction rate is valid

    !Reaction rates.
    REAL, SAVE, DIMENSION(nrec) :: &
            f, &            !Forward reaction rate coefficients.
            r, &            !Reverse reaction rate coefficients.
            fmax            !Store f(t9max)

    !Change in deutron binding energy
    REAL, PRIVATE, SAVE :: DeutronBindingVariation

    !Change in element binding energy
    LOGICAL, PRIVATE, SAVE :: BindingVariation

CONTAINS

    SUBROUTINE InitialiseReactions
        ! Read in reaction parameters and set rates to zero
        INTEGER :: i
        REAL :: GamowPrefactor, newQ
        INTEGER :: nuc1, nuc2

        f = 0.                 !Forward rate coeff.
        r = 0.                 !Reverse rate coeff.
        DeutronBindingVariation = 0.0
        BindingVariation = .false.
#if integrate_rates
        CALL InitialiseRateIntegrator
#endif
        !Gamow energy Eg = 2 pi^2 alpha^2 mu A (Z_1 Z_2)^2
        GamowPrefactor = 2. * pi**2 * alpha**2 * amu_MeV / kB
        EgIn  = 0.0
        EgOut = 0.0

        ! Set reaction rates to default values
        DO i  = 1,nrec
            iform(i) = int(reacpr(i,2))!Reaction type.
            ii(i)    = int(reacpr(i,3))!Incoming nuclide type.
            jj(i)    = int(reacpr(i,4))!Incoming nuclide type.
            kk(i)    = int(reacpr(i,5))!Outgoing nuclide type.
            ll(i)    = int(reacpr(i,6))!Outgoing nuclide type.
            rev(i)   = reacpr(i,7)     !Reverse reaction coefficient.
            q9(i)    = reacpr(i,8)     !Energy released.

            ! EgIn
            nuc1 = 0
            nuc2 = 0
            select case (iform(i))
            case (2, 3, 5, 9, 10)
                nuc1 = ii(i)
                nuc2 = jj(i)
            case (6, 11)
                nuc1 = ii(i)
                nuc2 = ii(i)
            end select
            if((nuc1.ne.0) .and. (nuc2.ne.0)) then
                EgIn(i) = GamowPrefactor * (zm(nuc1) * zm(nuc2))**2 / &
                    (1./(am(nuc1)+dm(nuc1)) + 1./(am(nuc2)+dm(nuc2)))
            end if

            ! EgOut
            nuc1 = 0
            nuc2 = 0
            select case (iform(i))
            case (3, 6)
                nuc1 = kk(i)
                nuc2 = ll(i)
            case (4, 5)
                nuc1 = ll(i)
                nuc2 = ll(i)
            end select
            if((nuc1.ne.0) .and. (nuc2.ne.0)) then
                EgOut(i) = GamowPrefactor * (zm(nuc1) * zm(nuc2))**2 / &
                    (1./(am(nuc1)+dm(nuc1)) + 1./(am(nuc2)+dm(nuc2)))
!                print *, i, EgOut(i)/q9(i)
            end if

            ! Check q values
            newQ = si(iform(i))*km(ii(i)) - sl(iform(i))*km(ll(i))
            if(sj(iform(i)).ne.0) newQ = newQ + sj(iform(i))*km(jj(i))
            if(sk(iform(i)).ne.0) newQ = newQ - sk(iform(i))*km(kk(i))
            newQ = newQ /1000.0 / kB
            
!            if(abs(newQ - q9(i)) .gt. 0.01) then
!                print *, i, q9(i), (q9(i)-newQ), (q9(i)-newQ)/newQ
!            end if
        END DO

        CALL InitialiseMaxima
    END SUBROUTINE InitialiseReactions


    SUBROUTINE InitialiseMaxima
        ! Set maximum temperatures and reaction rates
        INTEGER :: i
        t9max = T9LIMIT
        fmax  = 0.0
        ! Set maximum temperatures for important reactions
#if rates_2004C
        t9max(20) = 3.9
        t9max(28) = 12.5
        t9max(29) = 5.8
        t9max(16) = 3.9
        t9max(30) = 2.3
        t9max(31) = 2.3
        t9max(27) = 7.8
        t9max(26) = 3.9
        t9max(17) = 11.7
        t9max(24) = 3.9
#elif rates_2004D
        t9max(20) = 0.8
        t9max(28) = 3.0
        t9max(29) = 3.0
        t9max(30) = 0.5
        t9max(26) = 8.0
        t9max(16) = 3.0
        t9max(31) = 2.0
        t9max(27) = 8.0
        t9max(24) = 7.0
        t9max(17) = 0.2
#endif
        ! Unset for rates that are integrated explicitly
#if integrate_rates
#if integrate_12
        t9max(12) = T9LIMIT
#endif
#if integrate_20
        t9max(20) = T9LIMIT
#endif
#if integrate_28
        t9max(28) = T9LIMIT
#endif
#if integrate_29
        t9max(29) = T9LIMIT
#endif
#if integrate_16
        t9max(16) = T9LIMIT
#endif
#if integrate_30
        t9max(30) = T9LIMIT
#endif
#if integrate_31
        t9max(31) = T9LIMIT
#endif
#if integrate_27
        t9max(27) = T9LIMIT
#endif
#if integrate_26
        t9max(26) = T9LIMIT
#endif
#if integrate_17
        t9max(17) = T9LIMIT
#endif
#if integrate_24
        t9max(24) = T9LIMIT
#endif
#endif
        ! Calculate all maxima
        DO i = 12, jsize
            IF(t9max(i) .lt. T9LIMIT) THEN
                CALL Rate2(t9max(i))
                fmax(i) = f(i)
            ENDIF
        END DO
        DO i = 35, jsize
            IF(t9max(i) .lt. T9LIMIT) THEN
                CALL Rate3(t9max(i))
                fmax(i) = f(i)
            ENDIF
        END DO
        DO i = 64, jsize
            IF(t9max(i) .lt. T9LIMIT) THEN
                CALL Rate4(t9max(i))
                fmax(i) = f(i)
            ENDIF
        END DO
        f = 0.0     ! Reset forward rate coefficients
    END SUBROUTINE


    SUBROUTINE SetDeutronBindingVariation(dQ)
        ! Set change in deutron binding energy
        REAL, INTENT(IN) :: dQ
        INTEGER :: reaction, element

        element = 3
        DeutronBindingVariation = dQ
        BindingVariation = .true.

        do reaction = 12, jsize
            if(ii(reaction).eq.element) then
                q9(reaction) = q9(reaction) &
                    - si(iform(reaction)) * DeutronBindingVariation
            else if(jj(reaction).eq.element) then
                q9(reaction) = q9(reaction) &
                    - sj(iform(reaction)) * DeutronBindingVariation
            else if(kk(reaction).eq.element) then
                q9(reaction) = q9(reaction) &
                    + sk(iform(reaction)) * DeutronBindingVariation
            else if(ll(reaction).eq.element) then
                q9(reaction) = q9(reaction) &
                    + sl(iform(reaction)) * DeutronBindingVariation
            end if
        end do
        !q9(12) = reacpr(12, 8) + DeutronBindingVariation
    END SUBROUTINE


    SUBROUTINE SetBindingVariation(element, dQ)
        ! Set change in binding energy of element
        INTEGER, INTENT(IN) :: element
        REAL, INTENT(IN) :: dQ
        INTEGER :: reaction

        BindingVariation = .true.

        do reaction = 12, jsize
            if(ii(reaction).eq.element) then
                q9(reaction) = q9(reaction) &
                    - si(iform(reaction)) * dQ
            else if(jj(reaction).eq.element) then
                q9(reaction) = q9(reaction) &
                    - sj(iform(reaction)) * dQ
            else if(kk(reaction).eq.element) then
                q9(reaction) = q9(reaction) &
                    + sk(iform(reaction)) * dQ
            else if(ll(reaction).eq.element) then
                q9(reaction) = q9(reaction) &
                    + sl(iform(reaction)) * dQ
            end if
            ! Print relative change in Q
            if(element.eq.0 .and. q9(reaction).ne.reacpr(reaction, 8)) then
                select case (iform(reaction))
                case (3, 4, 5, 6)
                    print *, reaction, reacpr(reaction, 8), (q9(reaction)/reacpr(reaction, 8)-1.)*100. &
                           , sqrt(EgOut(reaction)/reacpr(reaction, 8))
                case (1, 2, 7, 8)
                    print *, reaction, reacpr(reaction, 8), (q9(reaction)/reacpr(reaction, 8)-1.)*100.
                end select
            end if
        end do

        ! Print percentage changes in forward reaction rates
        IF (element.eq.0) then
        DO reaction = 12, 34
            IF (q9(reaction) /= reacpr(reaction, 8)) THEN
                select case (iform(reaction))
                case (3, 4, 5, 6)
#if include_gamow_shift
                    print *, reaction, 0.5 * (1. + sqrt(EgOut(reaction)/reacpr(reaction, 8))) &
                                  * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8) * 100.
#else
                    print *, reaction, 0.5 * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8) * 100.
#endif
                case (1, 2, 7, 8)
                    print *, reaction, 3.0 * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8) * 100.
                end select
            END IF
        END DO
        END IF
    END SUBROUTINE


    SUBROUTINE SetHe5ResonancePositions(dEr)
        ! Shift resonances due to He5
        REAL, INTENT(IN) :: dEr
!  6    t(d,n)4He
#if integrate_30
        call VaryResonance(6, dEr, 0.0)
#endif
    END SUBROUTINE

    SUBROUTINE SetLi5ResonancePositions(dEr)
        ! Shift resonances due to Li5
        REAL, INTENT(IN) :: dEr
!  7  3He(d,p)4He
#if integrate_31
        call VaryResonance(7, dEr, 0.0)
#endif
    END SUBROUTINE

    SUBROUTINE SetHe5ResonanceWidth(dG)
        ! Shift resonances due to He5
        REAL, INTENT(IN) :: dG
!  6    t(d,n)4He
#if integrate_30
        call VaryResonance(6, 0.0, dG)
#endif
    END SUBROUTINE

    SUBROUTINE SetLi5ResonanceWidth(dG)
        ! Shift resonances due to Li5
        REAL, INTENT(IN) :: dG
!  7  3He(d,p)4He
#if integrate_31
        call VaryResonance(7, 0.0, dG)
#endif
    END SUBROUTINE


    SUBROUTINE PrintQValues
        ! Print Q and dQ for all reactions
        INTEGER :: reaction
        do reaction = 12, jsize
            if(q9(reaction).ne.reacpr(reaction, 8)) &
                write(*, "(i2,': ',f6.2,' -> ',f6.2,';  dQ = ',f6.2,' = ',f6.1,'%')") &
                    reaction, reacpr(reaction, 8), q9(reaction), q9(reaction)-reacpr(reaction, 8), &
                    (q9(reaction)-reacpr(reaction, 8))/reacpr(reaction, 8) * 100.0
        end do
    END SUBROUTINE


    SUBROUTINE Rate2(t9)

        !-------LINKAGES.
        !     CALLED BY - [subroutine] derivs

        !-------REMARKS.
        !     Generates rate coefficients for reactions involving nuclides
        !     up to A = 9.

        REAL, INTENT(IN) :: t9
        INTEGER :: reaction
        REAL :: temp

        !-------TEMPERATURE FACTORS--------------------------------

        REAL :: t92, t93, t913, t923, t943, t953, t912, t932, &
                t9m1, t9m23, t9m32, t9a, t9a32, t9b, t9b32, &
                t9c, t9c13, t9c56, t9d, t9d13, t9d56, &
                t9e, t9e13, t9e56, t9f, t9f13, t9f56

        t92   = t9*t9
        t93   = t9*t92
        t913  = t9**(.33333333)      !t9**(1/3)
        t923  = t913*t913            !t9**(2/3)
        t943  = t923*t923            !t9**(4/3)
        t953  = t9*t923              !t9**(5/3)
        t912  = sqrt(t9)             !t9**(1/2)
        t932  = t9*t912              !t9**(3/2)
        t9m1  = 1/t9                 !t9**(-1)
        t9m23 = 1.0/t923             !t9**(-2/3)
        t9m32 = 1.0/t932             !t9**(-3/2)
        t9a   = t9/(1.0+13.076*t9)   !For reaction 17.
        t9a32 = t9a**(1.5)           !t9a**(3/2)
        t9b   = t9/(1.+49.18*t9)     !For reaction 18.
        t9b32 = t9b**(1.5)           !t9b**(3/2)
        IF (t9.gt.10.) THEN          !For reaction 22.
            t9c = 1.
        ELSE
            t9c = t9/(1.-9.69e-2*t9+2.84e-2*t953/(1.-9.69e-2*t9)**(2./3.))
        END IF
        t9c13 = t9c**(.3333333)      !t9c**(1/3)
        t9c56 = t9c**(.8333333)      !t9c**(5/6)
        t9d   = t9/(1.+0.759*t9)     !For reaction 24.
        t9d13 = t9d**(.3333333)      !t9d**(1/3)
        t9d56 = t9d**(.8333333)      !t9d**(5/6)
        t9e   = t9/(1.+0.1378*t9)    !For reaction 26.
        t9e13 = t9e**(.3333333)      !t9e**(1/3)
        t9e56 = t9e**(.8333333)      !t9e**(5/6)
        t9f   = t9/(1.+0.1071*t9)    !For reaction 27.
        t9f13 = t9f**(.3333333)      !t9f**(1/3)
        t9f56 = t9f**(.8333333)      !t9f**(5/6)

        !-------NEUTRON, PHOTON REACTIONS-----------------------------

#if integrate_12
        !.......H(n,g)H2...................
        f(12)  = GetRate(1, t9)
#elif rate12_2006
        f(12)  = 44216.0 * (1. +3.75191*t9 +1.92934*t92 +0.746503*t93 &
                            +1.97023e-2*t92*t92 +3.00491e-6*t92*t93) &
                    / (1. +5.4678*t9 +5.62395*t92 +0.489312*t93 +0.00747806*t92*t92)
#elif rates_2004C
        f(12)  = 4.40654e4 * (1. +.0457518*t912 -2.47101*t9 +4.17185*t932 -3.44553*t92 &
                    +1.72766*t9*t932 -.546196*t93 +.106066*t912*t93 -.0115306*t92*t92 &
                    +.536436e-3*t932*t93)
#elif rates_1992
        !.......H(n,g)H2...................(Smith-Kawano-Malaney 1992)
        f(12)  = 4.742e+4*(1.-.8504*t912+.4895*t9-.09623*t932 &
                              +8.471e-3*t9*t9-2.80e-4*t9*t932)
#endif

        !.......H2(n,g)H3..................(Wagoner 1969)
        f(13)  = 6.62e+1*(1.+18.9*t9)

        !.......He3(n,g)He4................(Wagoner 1969)
        f(14)  = 6.62e+0*(1.+905.*t9)

        !.......Li6(n,g)Li7................(Malaney-Fowler 1989)
        f(15)  = 5.10e+3

        !-------NEUTRON, PROTON REACTIONS-----------------------------

        !.......He3(n,p)H3.................
#if integrate_16
        f(16)  = GetRate(5, t9)
#elif rates_2004D
        f(16)  = 6.505e+8 * (1. -0.655*t9 +0.445*t92 -0.082*t93)
#elif rates_2004C
        f(16)  = 6.84713e+8 * (1. -.0171094*t912 -2.66179*t9 +8.27463*t932 -14.3898*t92 &
                    +15.6385*t932*t9 -10.3337*t93 +3.80177*t912*t93 -.599790*t92*t92 &
                    -0.0139213*t932*t93 +.0140311*t92*t93 -.00106709*t932*t92*t92 &
                    +1.06709e-6*t93*t93)
#elif rates_1999
        !.......He3(n,p)H3.................(Cyburt-Fields-Olive 2001)
        f(16)  = 7.3546e+8*(1.-.7757*t912+.5376*t9-.1018*t932)
#elif rates_1992
        !.......He3(n,p)H3.................(Smith-Kawano-Malaney 1992)
        f(16)  = 7.21e+8*(1.-.508*t912+.228*t9)
#endif

        !.......Be7(n,p)Li7................
#if integrate_17
        f(17)  = GetRate(10, t9)
#elif rates_2004D
        f(17)  = 4.609e+9 * (1. -7.518*t9 +53.093*t92 -135.953*t93)
#elif rates_2004C
        f(17)  = 5.17900e9 * (1. -1.44587*t912 +1.12925*t9 -0.493526*t932 +0.126269*t92 &
                    -0.0194265*t932*t9 +.00177188*t93 -.883411e-4*t912*t93 +.185551e-5*t92*t92) &
                 + 4.2994e9*t9m32*exp(-3.713442/t9) &
                 + 1.36949e11*t9m32*exp(-31.332167/t9)
#elif rates_1992
        !.......Be7(n,p)Li7................(Smith-Kawano-Malaney 1992)
        f(17)  = 2.675e+9*(1.-.560*t912+.179*t9-.0283*t932 &
                 + 2.214e-3*t9*t9-6.851e-5*t9*t932) &
                 + 9.391e+8*t9a32*t9m32 &
                 + 4.467e+7*t9m32*exp(-0.07486/t9)
#endif

        !-------NEUTRON, ALPHA REACTIONS------------------------------

        !.......Li6(n,a)H3.................(Caughlan-Fowler 1988)
        f(18)  = 2.54e+9*t9m32*exp(-2.39/t9) &
                  + 1.68e+8*(1.-.261*t9b32/t932)

        !.......Be7(n,a)He4................(Wagoner 1969)
        f(19)  = 2.05e+4*(1.+3760.*t9)

        !-------PROTON, PHOTON REACTIONS------------------------------

#if integrate_20
        !.......H2(p,g)He3.................
        f(20)  = GetRate(2, t9)
#elif rates_2004D
        f(20)  = 2.173e3 * t9m23*exp(-3.7208/t913) * (1. +6.899*t9 -4.442*t92 +3.134*t93)
#elif rates_2004C
        f(20)  = 7.30909e3 * t9m23* exp(-3.7209/t913) * (1. -10.3497*t913 +63.4315*t923 -209.780*t9 &
                    +432.557*t943 -571.937*t953 +497.303*t92 -284.936*t943*t9 +106.863*t953*t9 &
                    -25.7496*t93 +3.81387*t913*t93 -.313823*t923*t93 +.0108908*t92*t92)
#elif rates_1999
        !.......H2(p,g)He3.................(NACRE)
        if(t9.gt.0.11) then
            f(20) = 2.58E+03*t9m23*exp(-3.721/t913)*(1.+3.96*t9+0.116*t9**2)
        else
            f(20) = 1.81e+3*t9m23*exp(-3.721/t913) &
                    *(1.+14.3*t9 -90.5*t9**2 + 395.*t9**3)
        end if
#elif rates_1992
        !.......H2(p,g)He3.................(Smith-Kawano-Malaney 1992)
        f(20)  = 2.65e+3*t9m23*exp(-3.720/t913) &
                  *(1.+.112*t913+1.99*t923+1.56*t9+.162*t943+.324*t953)
#endif

        !.......H3(p,g)He4.................(Caughlan-Fowler 1988)
        f(21)  = 2.20e+4*t9m23*exp(-3.869/t913) &
                  *(1.+.108*t913+1.68*t923+1.26*t9+.551*t943+1.06*t953)

#if rates_1999
        !.......Li6(p,g)Be7................(NACRE)
        f(22)  = 1.25e+6*t9m23*exp(-8.415/t913) &
                 *(1.-0.252*t9+5.19e-2*t9**2-2.92e-3*t9**3)
#elif rates_1992
        !.......Li6(p,g)Be7................(Caughlan-Fowler 1988)
        f(22)  = 6.69e+5*t9c56*t9m32*exp(-8.413/t9c13)
#endif

        !-------PROTON, ALPHA REACTIONS-------------------------------

#if rates_1999
        !.......Li6(p,a)He3................(NACRE)
        f(23)  = 3.54E+10*t9m23*exp(-8.415/t913) &
                 *(1.-0.137*t9 + 2.41E-02*t9**2-1.28E-03*t9**3)
#elif rates_1992
        !.......Li6(p,a)He3................(Caughlan-Fowler 1988)
        f(23)  = 3.73e+10*t9m23*exp(-8.413/t913-(t9/5.50)**2) &
                  *(1.+.050*t913-.061*t923-.021*t9+.006*t943+.005*t953) &
                  + 1.33e+10*t9m32*exp(-17.763/t9) &
                  + 1.29e+09*t9m1*exp(-21.820/t9)
#endif

        !.......Li7(p,a)He4................
#if integrate_24
        f(24)  = GetRate(11, t9)
#elif rates_2004D
        f(24)  = 8.309e8 * t9m23*exp(-8.4727/t913) * (1. + 0.278*t9 -0.018*t92 +0.005*t93)
#elif rates_2004C
        f(24)  = 9.19322e8 * t9m23*exp(-8.4730/t913) * (1. -2.26222*t913 +11.3224*t923 -27.3071*t9 &
                    +41.1901*t943 -37.4242*t953 +18.3941*t92 -3.72281*t9*t943 +2.58125e-2*t953*t9)
#elif rates_1999
        !.......Li7(p,a)He4................(NACRE)
        f(24)  = 7.20E+08*t9m23*exp(-8.473/t913-(t9/6.5)**2) &
                 *(1. + 1.05*t9 -0.653*t9**2 &
                   + 0.185*t9**3 -2.12E-02*t9**4 + 9.30E-04*t9**5) &
                 + 9.85E+06*t9**0.576*exp(-10.415/t9)
#elif rates_1992
        !.......Li7(p,a)He4................(Smith-Kawano-Malaney 1992)
        f(24)  = 1.096e+9*t9m23*exp(-8.472/t913) &
                  - 4.830e+8*t9d56*t9m32*exp(-8.472/t9d13) &
                  + 1.06e+10*t9m32*exp(-30.442/t9) &
                  + 1.56e+5*t9m23*exp((-8.472/t913)-(t9/1.696)**2) &
                    *(1.+.049*t913-2.498*t923+.860*t9+3.518*t943+3.08*t953) &
                  + 1.55e+6*t9m32*exp(-4.478/t9)
#endif

        !-------ALPHA, PHOTON REACTIONS-------------------------------

#if rates_1999
        !.......H2(a,g)Li6.................(NACRE)
        f(25)  = 14.82*t9m23*exp(-7.435/t913) &
                 *(1. + 6.572*t9 + 0.076*t9**2 + 0.0248*t9**3) &
                 + 82.8/t932*exp(-7.904/t9)
!        f(25) = f(25) * 1.e3
#elif rates_1992
        !.......H2(a,g)Li6.................(Caughlan-Fowler 1988)
        f(25)  = 3.01e+01*t9m23*exp(-7.423/t913) &
                  *(1.+.056*t913-4.85*t923+8.85*t9-.585*t943-.584*t953) &
                  + 8.55e+1*t9m32*exp(-8.228/t9)
#endif

        !.......H3(a,g)Li7.................
#if integrate_26
        f(26)  = GetRate(9, t9)
#elif rates_2004D
        f(26)  = 7.717e5 * t9m23*exp(-8.0805/t913) * (1. -0.268*t9 +0.068*t92 -0.004*t93)
#elif rates_2004C
        f(26)  = 4.65494351e6 * t9m23*exp(-8.0808/t913) * (1. -12.3956341*t913 +76.2717899*t923 &
                    -250.678479*t9 +446.413119*t943 -289.008201*t953 -474.786707*t92 +1346.42142*t9*t943 &
                    -1503.09444*t9*t953 +923.138882*t93 -306.14089*t913*t93 +42.9886919*t923*t93)
#elif rates_1999
        !.......H3(a,g)Li7.................(NACRE)
        f(26)  = 8.20E+5*t9m23*exp(-8.081/t913) &
                 * (1.-0.389*t9+0.134*t9**2-0.0181*t9**3+9.23E-4*t9**4)
#elif rates_1992
        !.......H3(a,g)Li7.................(Smith-Kawano-Malaney 1992)
        f(26)  = 3.032e+5*t9m23*exp(-8.090/t913) &
                  *(1.+.0516*t913+.0229*t923+8.28e-3*t9 &
                      -3.28e-4*t943-3.01e-4*t953) &
                  + 5.109e+5*t9e56*t9m32*exp(-8.068/t9e13)
#endif

#if integrate_27
        !.......He3(a,g)Be7................(Cyburt-Davids 2008)
        f(27)  = GetRate(8, t9)
#elif rate27_2008
        f(27)  = 6015602. * t9m23 * exp(-12.82707707/t913) &
                 *(1. -2.0478E-2*t923 +0.211995*t943) &
                 /(1. +0.255059*t923 + 0.338573*t943)
#elif rate27_2008_refit
        f(27)  = 6015602. * t9m23 * exp(-12.82707707/t913) &
                 *(1. +0.454877*t923 +0.244011*t943) &
                 /(1. +0.727006*t923 +0.527836*t943)
#elif rates_2004D
        f(27)  = 5.216e6 * t9m23*exp(-12.827/t913) * (1. -0.235*t9 +0.041*t92 -0.002*t93)
#elif rates_2004C
        f(27)  = 3.94207e6 * t9m23*exp(-12.8274/t913) * (1. +.185267*t913 -.837432*t923 +7.23019*t9 &
                    -26.1976*t943 +41.6914*t953 +19.4465*t92 -215.248*t9*t943 +422.548*t9*t953 &
                    -412.866*t93 +176.691*t913*t93 +45.8891*t923*t93 -100.644*t92*t92 &
                    +54.8984*t943*t93 -14.1903*t953*t93 +1.48464*t92*t93)
#elif rates_1999
        !.......He3(a,g)Be7................(NACRE)
        f(27)  = 5.46E+06*t9m23*exp(-12.827/t913) &
                 *(1.-0.307*t9+8.81E-2*t9**2-1.06E-2*t9**3+4.46E-4*t9**4)
#elif rates_1992
        !.......He3(a,g)Be7................(Smith-Kawano-Malaney 1992)
        f(27)  = 4.817e+6*t9m23*exp(-14.964/t913) &
                  *(1.+.0325*t913-1.04e-3*t923-2.37e-4*t9 &
                      -8.11e-5*t943-4.69e-5*t953) &
                  + 5.938e+6*t9f56*t9m32*exp(-12.859/t9f13)
#endif

        !-------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------

        !.......H2(d,n)He3.................
#if integrate_28
        f(28)  = GetRate(3, t9)
#elif rates_2004D
        f(28)  = 4.371e8 * t9m23*exp(-4.2586/t913) * (1. +1.737*t9 -0.633*t92 +0.109*t93)
#elif rates_2004C
        f(28)  = 1.00749e9 * t9m23*exp(-4.2586/t913) * (1. -9.59015*t913 +65.2448*t923 -247.756*t9 &
                    +596.231*t943 -941.064*t953 +980.076*t92 -643.032*t9*t943 +211.982*t9*t953 &
                    +29.0491*t93 -66.1847*t913*t93 +31.6452*t923*t93 -7.15147*t92*t92 &
                    +.372749*t943*t93 +.208645*t953*t93 -.0545129*t92*t93 +.00536216*t913*t92*t93 &
                    -.000157984*t923*t92*t93 -.457514e-5*t93*t93 +2.123592e-9*t913*t93*t93)
#elif rates_1999
        !.......H2(d,n)He3.................(NACRE)
        f(28)  = 4.67E+8*t9m23*exp(-4.259/t913) &
                 * (1. + 1.079*t9 - 0.1124*t9**2 + 5.68E-3*t9**3)
#elif rates_1992
        !.......H2(d,n)He3.................(Smith-Kawano-Malaney 1992)
        f(28)  = 3.95e+8*t9m23*exp(-4.259/t913) &
                  *(1.+.098*t913+.765*t923+.525*t9+9.61e-3*t943+.0167*t953)
#endif

        !.......H2(d,p)H3..................
#if integrate_29
        f(29)  = GetRate(4, t9)
#elif rates_2004D
        f(29)  = 4.682e8 * t9m23*exp(-4.2586/t913) * (1. +0.745*t9 -0.065*t92 +0.003*t93)
#elif rates_2004C
        f(29)  = 3.91889e8 * t9m23* exp(-4.2586/t913) * (1. +.309233*t913 -.337260*t923 &
                    +2.51922*t9 -2.79097*t943 +2.16082*t953 -.976181*t92 +.210883*t943*t9 &
                    -.0169027*t953*t9 +7.845538e-6*t93)
#elif rates_1999
        !.......H2(d,p)H3..................(NACRE)
        f(29)  = 4.66E+8*t9m23*exp(-4.259/t913) &
                 * (1. + 0.759*t9 - 0.0612*t9**2 + 2.78E-3*t9**3)
#elif rates_1992
        !.......H2(d,p)H3..................(Smith-Kawano-Malaney 1992)
        f(29)  = 4.17e+8*t9m23*exp(-4.258/t913) &
                 *(1.+.098*t913+.518*t923+.355*t9-.010*t943-.018*t953)
#endif

        !.......H3(d,n)He4.................
#if integrate_30
        f(30)  = GetRate(6, t9)
#elif rates_2004D
        f(30)  = 8.656e10 * t9m23*exp(-4.5244/t913) * (1. +14.002*t9 -59.683*t92 +64.236*t93)
#elif rates_2004C
        f(30)  = 1.78988e12 * t9m23*exp(-4.5245/t913) &
                / (1. + ((0.129964*t923 -0.0482)/(0.5*0.0806))**2) &
                * (1. -14.3137899*t913 +92.4325675*t923 -314.645738*t9 +641.100355*t943 -844.106855*t953 &
                    +752.418564*t92 -465.820564*t943*t9 +202.276143*t953*t9 -61.3172473*t93 &
                    +12.6913874*t913*t93 -1.707344*t923*t93 +.134399048*t92*t92 &
                    -.00469341945*t943*t93)
#elif rates_1999
        !.......H3(d,n)He4.................(NACRE)
        f(30)  = 8.29E+10*t9m23*exp(-4.524/t913 - (t9/0.08)**2) &
                 *(1.+17.2*t9+175.*t9**2) &
                 + 8.12E+8*t9**-0.712*exp(-0.506/t9)
#elif rates_1992
        !.......H3(d,n)He4.................(Smith-Kawano-Malaney 1992)
        f(30)  = 1.063e+11*t9m23*exp(-4.559/t913-(t9/.0754)**2) &
                  *(1.+.092*t913-.375*t923-.242*t9+33.82*t943+55.42*t953) &
                  + 8.047e+8*t9m23*exp(-0.4857/t9)
#endif

        !.......He3(d,p)He4................
#if integrate_31
        f(31)  = GetRate(7, t9)
#elif rates_2004D
        f(31)  = 5.477e10 * t9m23*exp(-7.1820/t913) * (1. +4.367*t9 -4.329*t92 +1.115*t93)
#elif rates_2004C
        f(31)  = 5.67897e12 * t9m23*exp(-7.1840/t913) &
                / (1. + ((0.206357*t923-0.183)/(0.5*0.256))**2) &
                * (1. -8.59410908*t913 +31.1979775*t923 -61.2218616*t9 +72.0331037*t943 -52.8696341*t953 &
                    +23.7371543*t92 -5.4569107*t943*t9 -.226478266*t953*t9 +.583380161*t93 &
                    -.190978484*t913*t93 +.031949394*t923*t93 &
                    -.00284146599*t92*t92 +.106749198e-3*t943*t93)
#elif rates_1992
        !.......He3(d,p)He4................(Smith-Kawano-Malaney 1992)
        f(31)  = 5.021e+10*t9m23*exp(-7.144/t913-(t9/.270)**2) &
                  *(1.+.058*t913+.603*t923+.245*t9+6.97*t943+7.19*t953) &
                  + 5.212e+8/t912*exp(-1.762/t9)
#endif

        !-------THREE PARTICLE REACTIONS------------------------------

#if rates_1999
        !.......He3(He3,2p)He4.............(NACRE)
        f(32)  = 5.59E+10*t9m23*exp(-12.277/t913) &
                 *(1. - 0.135*t9 + 2.54E-2*t9**2 -1.29E-03*t9**3)
#elif rates_1992
        !.......He3(He3,2p)He4.............(Caughlan-Fowler 1988)
        f(32)  = 6.04e+10*t9m23*exp(-12.276/t913) &
                  *(1.+.034*t913-.522*t923-.124*t9+.353*t943+.213*t953)
#endif

        !.......Li7(d,na)He4...............(Caughlan-Fowler 1988)
        f(33)  = 2.92e+11*t9m23*exp(-10.259/t913)

        !.......Be7(d,pa)He4...............(Caughlan-Fowler 1988)
        f(34)  = 1.07e+12*t9m23*exp(-12.428/t913)

#if linear_q_variation
        IF (DeutronBindingVariation /= 0.0) THEN
            f(12) = f(12) * (1.0 + (2.5 + sqrt(reacpr(12, 8)/0.812)) * DeutronBindingVariation/reacpr(12, 8))
!            f(12) = f(12) * (1.0 + 2.5  * DeutronBindingVariation/reacpr(12, 8))
        END IF
#else
        IF (DeutronBindingVariation /= 0.0) THEN
            temp = q9(12)/reacpr(12, 8)
            f(12) = f(12) * temp * temp * sqrt(temp)
            f(12) = f(12) * (1.0 + sqrt(reacpr(12, 8)/0.812) * DeutronBindingVariation/reacpr(12, 8))
        END IF
#endif

#if linear_q_variation
        IF (BindingVariation) THEN
        DO reaction = 13, 34
            IF (q9(reaction) /= reacpr(reaction, 8)) THEN
                select case (iform(reaction))
                case (3, 4, 5, 6)
#if include_gamow_shift
                    f(reaction) = f(reaction) * (1.0 + 0.5 * (1. + sqrt(EgOut(reaction)/reacpr(reaction, 8))) &
                                  * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8))
#else
                    f(reaction) = f(reaction) * (1.0 + 0.5 * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8))
#endif
                case (1, 2, 7, 8)
                    f(reaction) = f(reaction) * (1.0 + 3.0 * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8))
                end select
            END IF
        END DO
        END IF
#else
        IF (BindingVariation) THEN
        DO reaction = 13, 34
            IF (q9(reaction) .le. 0.0) THEN
                f(reaction) = 0.0
            ELSE IF (q9(reaction) /= reacpr(reaction, 8)) THEN
                select case (iform(reaction))
                case (3, 4, 5, 6)
#if include_gamow_shift
                    f(reaction) = f(reaction) * sqrt(q9(reaction)/reacpr(reaction, 8)) &
                                  * exp(-sqrt(EgOut(reaction)/q9(reaction)) + sqrt(EgOut(reaction)/reacpr(reaction, 8)))
#else
                    f(reaction) = f(reaction) * sqrt(q9(reaction)/reacpr(reaction, 8))
#endif
                case (1, 2, 7, 8)
                    temp = q9(reaction)/reacpr(reaction, 8)
                    f(reaction) = f(reaction) * temp * temp * temp
                end select
            END IF
        END DO
        END IF
#endif
        RETURN

    !-------REFERENCES--------------------------------------
    !     Smith, M., Kawano, L.H., and Malaney, R.A., 1992, submitted to Ap. J.
    !     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5.
    !     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data 
    !       Tables, 40, 283.
    !     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.
    !     NACRE    http://pntpm3.ulb.ac.be/Nacre/
    !     Cyburt R.H., Fields B.D., Olive K.A., 2001, New Astronomy, 6, 215

    END SUBROUTINE Rate2


    SUBROUTINE Rate3(t9)

        !-------LINKAGES.
        !     CALLED BY - [subroutine] derivs

        !-------REMARKS.
        !     Generates rate coefficients for reactions involving nuclides 
        !     up to A = 18.

        REAL, INTENT(IN) :: t9
        INTEGER :: reaction
        REAL :: temp

        !-------TEMPERATURE FACTORS--------------------------------

        REAL :: t913, t923, t943, t953, t912, t932, t915, t954, &
                t9m1, t9m23, t9m32, t9m34, t9m15, t9m54, &
                t9a, t9a13, t9a56, fff

        t913  = t9**(.33333333)      !t9**(1/3)
        t923  = t913*t913            !t9**(2/3)
        t943  = t923*t923            !t9**(4/3)
        t953  = t9*t923              !t9**(5/3)
        t912  = sqrt(t9)             !t9**(1/2)
        t932  = t9*t912              !t9**(3/2)
        t915  = t9**(.2)             !t9**(1/5)
        t954  = t9**(1.25)           !t9**(5/4)
        t9m1  = 1.0/t9               !t9**(-1)
        t9m23 = 1.0/t923             !t9**(-2/3)
        t9m32 = 1.0/t932             !t9**(-3/2)
        t9m34 = sqrt(t9m32)          !t9**(-3/4)
        t9m15 = 1.0/t915             !t9**(-1/5)
        t9m54 = 1.0/t954             !t9**(-5/4)
        t9a   = t9/(1.+t9/15.1)      !For reaction 53.
        t9a13 = t9a**(.3333333)      !t9a**(1/3)
        t9a56 = t9a**(.8333333)      !t9a**(5/6)

        !-------NEUTRON, PHOTON REACTIONS-----------------------------

        !.......Li7(n,g)Li8................(Wagoner 1969)
        f(35)  = 4.90e+3 + 9.96e+3*t9m32*exp(-2.62/t9)

        !.......B10(n,g)B11................(Wagoner 1969)
        f(36)  = 6.62e+4

        !.......B11(n,g)B12................(Malaney-Fowler 1989)
        f(37)  = 7.29e+2 + 2.40e+3*t9m32*exp(-0.223/t9)

        !-------NEUTRON, PROTON REACTIONS-----------------------------

        !.......C11(n,p)B11................(Caughlan-Fowler 1988)
        f(38)  = 1.69e+8*(1.-.048*t912+.010*t9)

        !-------NEUTRON, ALPHA REACTIONS------------------------------

        !.......B10(n,a)Li7................(Caughlan-Fowler 1988)
        f(39)  = 5.07e+8

        !-------PROTON, PHOTON REACTIONS------------------------------

#if rates_1999
        !.......Be7(p,g)B8.................(NACRE)
        f(40)  = 2.61e+5*t9m23*exp(-10.264/t913) &
                 *(1.-5.11e-2*t9+4.68e-2*t9**2-6.6e-3*t9**3+3.12e-4*t9**4) &
                 +2.05e+3/t932*exp(-7.345/t9)
#elif rates_1992
        !.......Be7(p,g)B8.................(Caughlan-Fowler 1988)
        f(40)  = 3.11e+5*t9m23*exp(-10.262/t913) &
                  + 2.53e+3*t9m32*exp(-7.306/t9)
#endif

#if rates_1999
        !.......Be9(p,g)B10................(NACRE)
        f(41)  = 1.36e+7*t9m23*exp(-10.361/t913 -(t9/1.5)**2) &
                 * (1. + 2.71*t9 - 1.95*t9**2 + 0.594*t9**3) &
                 + 4.8e+3/t932*exp(-3.102/t9)+2.75e+6/t932*exp(-10.615/t9)
#elif rates_1992
        !.......Be9(p,g)B10................(Caughlan-Fowler 1988)
        f(41)  = 1.33e+7*t9m23*exp(-10.359/t913-(t9/.846)**2) &
                  *(1.+.040*t913+1.52*t923+.428*t9+2.15*t943+1.54*t953) &
                  + 9.64e+4*t9m32*exp(-3.445/t9) &
                  + 2.72e+6*t9m32*exp(-10.620/t9)
#endif

#if rates_1999
        !.......B10(p,g)C11................(NACRE)
        f(42)  = 1.68e+6*t9m23*exp(-12.064/t913) &
                 * (1.+ 0.977*t9+1.87*t9**2-0.272*t9**3+1.3e-2*t9**4) &
                 / ((t923-0.0273)**2 + 4.69e-4)
#elif rates_1992
        !.......B10(p,g)C11................(Caughlan-Fowler 1988)
        f(42)  = 4.61e+5*t9m23*exp(-12.062/t913-(t9/4.402)**2) &
                  *(1.+.035*t913+.426*t923+.103*t9+.281*t943+.173*t953) &
                  + 1.93e+5*t9m32*exp(-12.041/t9) &
                  + 1.14e+4*t9m32*exp(-16.164/t9)
#endif

#if rates_1999
        !.......B11(p,g)C12................(NACRE)
        f(43)  = 4.58e+7*t9m23*exp(-12.097/t913-(t9/0.6)**2) &
                 *(1.+0.353*t9+0.842*t9**2) &
                 + 6.82e+3/t932*exp(-1.738/t9) &
                 + 2.8e+4*t9**0.104*exp(-3.892/t9)
#elif rates_1992
        !.......B11(p,g)C12................(Caughlan-Fowler 1988)
        f(43)  = 4.62e+7*t9m23*exp(-12.095/t913-(t9/.239)**2) &
                  *(1.+.035*t913+3.00*t923+.723*t9+9.91*t943+6.07*t953) &
                  + 7.89e+3*t9m32*exp(-1.733/t9) &
                  + 9.68e+4*t9m15*exp(-5.617/t9)
#endif

        !.......C11(p,g)N12................(Caughlan-Fowler 1988)
        f(44)  = 4.24e+4*t9m23*exp(-13.658/t913-(t9/1.627)**2) &
                  *(1.+.031*t913+3.11*t923+.665*t9+4.61*t943+2.50*t953) &
                  + 8.84e+3*t9m32*exp(-7.021/t9)

        !-------PROTON, NEUTRON REACTIONS-----------------------------

        !.......B12(p,n)C12................(Wagoner 1969)
        f(45)  = 4.02e+11*t9m23*exp(-12.12/t913)

        !-------PROTON, ALPHA REACTIONS-------------------------------

#if rates_1999
        !.......Be9(p,a)Li6................(NACRE)
        f(46)  = 2.11e+11*t9m23*exp(-10.361/t913-(t9/0.4)**2) &
                 *(1.-0.189*t9+35.2*t9**2) &
                 + 5.24e+8/t932*exp(-3.446/t9) &
                 + 4.65e+8*t9**(-0.293)*exp(-4.396/t9)
#elif rates_1992
        !.......Be9(p,a)Li6................(Caughlan-Fowler 1988)
        f(46)  = 2.11e+11*t9m23*exp(-10.359/t913-(t9/.520)**2) &
                  *(1.+.040*t913+1.09*t923+.307*t9+3.21*t943+2.30*t953) &
                  + 4.51e+8*t9m1*exp(-3.046/t9) &
                  + 6.70e+8*t9m34*exp(-5.160/t9)
#endif

#if rates_1999
        !.......B10(p,a)Be7................(NACRE)
        if(t9.gt.0.8) then
            f(47)  = 1.01e+10*t9m23*exp(-12.064/t913) &
                     * (-1.+15.8*t9-2.6*t9**2+0.125*t9**3)
        else
            f(47)  = 2.56e+10*t9m23*exp(-12.064/t913) &
                     * (1.+5.95*t9+29.2*t9**2-316.*t9**3+914.*t9**4 &
                        -1085.*t9**5+465.*t9**6)/((t923-0.026)**2+4.7e-4)
        end if
#elif rates_1992
        !.......B10(p,a)Be7................(Caughlan-Fowler 1988)
        f(47)  = 1.26e+11*t9m23*exp(-12.062/t913-(t9/4.402)**2) &
                  *(1.+.035*t913-.498*t923-.121*t9+.300*t943+.184*t953) &
                  + 2.59e+9*t9m1*exp(-12.260/t9)
#endif

        !.......B12(p,a)Be9................(Wagoner 1969)
        f(48)  = 2.01e+11*t9m23*exp(-12.12/t913)

        !-------ALPHA, PHOTON REACTIONS-------------------------------

        !.......Li6(a,g)B10................(Caughlan-Fowler 1988)
        f(49)  = 4.06e+6*t9m23*exp(-18.790/t913-(t9/1.326)**2) &
                  *(1.+.022*t913+1.54*t923+.239*t9+2.20*t943+.869*t953) &
                  + 1.91e+3*t9m32*exp(-3.484/t9) &
                  + 1.01e+4*t9m1*exp(-7.269/t9)

#if rates_1999
        !.......Li7(a,g)B11................(NACRE)
        f(50)  = 9.72e+7*t9m23*exp(-19.163/t913-(t9/0.4)**2) &
                 *(1.+2.84*t9-7.89*t9**2) &
                 + 3.35e+2/t932*exp(-2.959/t9) &
                 + 1.04e+4*t9**(-0.023)*exp(-4.922/t9)
#elif rates_1992
        !.......Li7(a,g)B11................(Caughlan-Fowler 1988)
        f(50)  = 3.55e+7*t9m23*exp(-19.161/t913-(t9/4.195)**2) &
                  *(1.+.022*t913+.775*t923+.118*t9+.884*t943+.342*t953) &
                  + 3.33e+2*t9m32*exp(-2.977/t9) &
                  + 4.10e+4*t9m1*exp(-6.227/t9)
#endif

#if rates_1999
!.......Be7(a,g)C11................(NACRE)
        if(t9.gt.2.) then
            f(51)  = 1.41e+3*t9**0.636*exp(-3.015/t9)
        else
            f(51)  = 1.29e+10*t9m23*exp(-23.214/t913-(t9/0.8)**2) &
                     * (1.-6.47*t9+19.5*t9**2-19.3*t9**3) &
                     + 1.25e+4/t932*exp(-6.498/t9)+1.44e+5/t932*exp(-10.177/t9) &
                     + 1.63e+4*t9**0.178*exp(-15.281/t9)
        end if
#elif rates_1992
        !.......Be7(a,g)C11................(Caughlan-Fowler 1988)
        f(51)  = 8.45e+7*t9m23*exp(-23.212/t913-(t9/4.769)**2) &
                  *(1.+.018*t913+.488*t923+.061*t9+.296*t943+.095*t953) &
                  + 1.25e+4*t9m32*exp(-6.510/t9) &
                  + 1.29e+5*t9m54*exp(-10.039/t9)
#endif

        !-------ALPHA, PROTON REACTIONS-------------------------------

        !.......B8(a,p)C11.................(Wagoner 1969)
        f(52)  = 1.08e+15*t9m23*exp(-27.36/t913)

        !----------ALPHA, NEUTRON REACTIONS------------------------------

        !.......Li8(a,n)B11................(Malaney-Fowler 1989)
        f(53)  = 8.62e+13*t9a56*t9m32*exp(-19.461/t9a13)

#if rates_1999
        !.......Be9(a,n)C12................(NACRE)
        f(54)  = 5.e+13*t9m23*exp(-23.872/t913-(t9/0.154)**2) &
                 * (1.+27.3*t9+1632.*t9**2) &
                 + 0.7/t932*exp(-1.832/t9)+1.77e+5/t932*exp(-4.385/t9) &
                 + 4.12e+7*t9**0.65*exp(-10.06/t9)
#elif rates_1992
        !.......Be9(a,n)C12................(Caughlan-Fowler 1988)
        f(54)  = 4.62e+13*t9m23*exp(-23.870/t913-(t9/.049)**2) &
                  *(1.+.017*t913+8.57*t923+1.05*t9+74.51*t943+23.15*t953) &
                  + 7.34e-5*t9m32*exp(-1.184/t9) &
                  + 2.27e-1*t9m32*exp(-1.834/t9) &
                  + 1.26e+5*t9m32*exp(-4.179/t9) &
                  + 2.40e+8*exp(-12.732/t9)
#endif

        !----------DEUTERIUM, NEUTRON AND DEUTERIUM, PROTON REACTIONS-------------

        !.......Be9(d,n)B10................(original Wagoner code)
        f(55)  = 7.16e+8*t9m23*exp(6.44-12.6/t913)

        !.......B10(d,p)B11................(original Wagoner code)
        f(56)  = 9.53e+8*t9m23*exp(7.30-14.8/t913)

        !.......B11(d,n)C12................(original Wagoner code)
        f(57)  = 1.41e+9*t9m23*exp(7.40-14.8/t913)

        !----------THREE PARTICLE REACTIONS------------------------------

#if rates_1999
        !.......He4(an,g)Be9...............(NACRE)
        fff  = 2.43e+9*t9m23*exp(-13.49/t913-(t9/0.15)**2) &
               *(1.+74.5*t9)+6.09e+5/t932*exp(-1.054/t9)
        if(t9.gt.0.03) then
            f(58) = fff*2.42e-12*(1-1.52*log10(t9)+0.448*log10(t9)**2 &
                    + 0.435*log10(t9)**3)
        else
            f(58) = fff*6.69e-12*(1.-192.*t9+2.48e+4*t9**2-1.5e+6*t9**3 &
                    + 4.13e+7*t9**4 - 3.9e+8*t9**5)
        end if
#elif rates_1992
        !.......He4(an,g)Be9...............(Caughlan-Fowler 1988)
        f(58)  = (2.59e-6/((1.+.344*t9)*t9**2))*exp(-1.062/t9)
#endif

#if rates_1999
        !.......He4(2a,g)C12...............(NACRE)
        f(59)  = 2.76e+7*t9m23*exp(-23.57/t913-(t9/0.4)**2) &
                 * (1.+5.47*t9+326*t9**2)+130.7/t932*exp(-3.338/t9) &
                 + 2.51e+4/t932*exp(-20.307/t9)
        if(t9.gt.0.03) then     ! fff defined in rate 58
            f(59) = fff*f(59)*3.44e-16*(1.+0.0158*t9**(-0.65))
        else
            f(59) = fff*f(59)*3.07e-16*(1.-29.1*t9 +1308.*t9**2)
        end if
#elif rates_1992
        !.......He4(2a,g)C12...............(Caughlan-Fowler 1988)
        f(59)  = 2.79e-8*t9m32*t9m32*exp(-4.4027/t9) &
                  + 1.35e-8*t9m32*exp(-24.811/t9)
#endif

        !.......Li8(p,na)He4...............(original Wagoner code)
        f(60)  = 8.65e+9*t9m23*exp(-8.52/t913-(t9/2.53)**2) &
                  + 2.31e+9*t9m32*exp(-4.64/t9)

        !.......B8(n,pa)He4................(original Wagoner code)
        f(61)  = 4.02e+8

#if rates_1999
        !.......Be9(p,da)He4...............(NACRE)
        f(62)  = 2.18e+11*t9m23*exp(-10.361/t913-(t9/0.42)**2) &
                 * (1.-0.427*t9+34.055*t9**2)+6.24e+8/t932*exp(-3.446/t9) &
                 + 3.53e+8*t9**(-0.205)*exp(-3.889/t9)
#elif rates_1992
        !.......Be9(p,da)He4...............(Caughlan-Fowler 1988)
        f(62)  = 2.11e+11*t9m23*exp(-10.359/t913-(t9/.520)**2) &
                  *(1.+.040*t913+1.09*t923+.307*t9+3.21*t943+2.30*t953) &
                  + 5.79e+8*t9m1*exp(-3.046/t9) &
                  + 8.50e+8*t9m34*exp(-5.800/t9)
#endif

#if rates_1999
        !.......B11(p,2a)He4...............(NACRE)
        if(t9.gt.2.) then
            f(63)  = 5.84e+11*t9m23*exp(-12.097/t913) &
                     *(-1.+0.883*t9+0.012*t9**2)/((t923-1.47)**2+0.187)
        else
            f(63)  = 2.68e+12*t9m23*exp(-12.097/t913) &
                     *(1.+1.62*t9-1.31*t9**2+0.26*t9**3) &
                     + 2.12e+6/t932*exp(-1.724/t9)
        end if
#elif rates_1992
        !.......B11(p,2a)He4...............(Caughlan-Fowler 1988)
        f(63)  = 2.20e+12*t9m23*exp(-12.095/t913-(t9/1.644)**2) &
                  *(1.+.034*t913+.140*t923+.034*t9+.190*t943+.116*t953) &
                  + 4.03e+6*t9m32*exp(-1.734/t9) &
                  + 6.73e+9*t9m32*exp(-6.262/t9) &
                  + 3.88e+9*t9m1*exp(-14.154/t9)
#endif

        !.......C11(n,2a)He4...............(Wagoner 1969)
        f(64)  = 1.58e+8

#if linear_q_variation
        IF (BindingVariation) THEN
        DO reaction = 35, 64
            IF (q9(reaction) /= reacpr(reaction, 8)) THEN
                IF((iform(reaction).eq.3).or.(iform(reaction).eq.5).or.(iform(reaction).eq.6)) THEN
                    f(reaction) = f(reaction) * (1.0 + 0.5 * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8))
                ELSE IF((iform(reaction).eq.2).or. &
                        (iform(reaction).eq.7).or.(iform(reaction).eq.8)) THEN
                    f(reaction) = f(reaction) * (1.0 + 3.0 * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8))
                END IF
            END IF
        END DO
        END IF
#else
        IF (BindingVariation) THEN
        DO reaction = 35, 64
            IF (q9(reaction) .le. 0.0) THEN
                f(reaction) = 0.0;
            ELSE IF (q9(reaction) /= reacpr(reaction, 8)) THEN
                IF((iform(reaction).eq.3).or.(iform(reaction).eq.5).or.(iform(reaction).eq.6)) THEN
                    f(reaction) = f(reaction) * sqrt(q9(reaction)/reacpr(reaction, 8))
                ELSE IF((iform(reaction).eq.2).or. &
                        (iform(reaction).eq.7).or.(iform(reaction).eq.8)) THEN
                    temp = q9(reaction)/reacpr(reaction, 8)
                    f(reaction) = f(reaction) * temp * temp * temp
                END IF
            END IF
        END DO
        END IF
#endif
        RETURN

    !-------REFERENCES--------------------------------------
    !     Malaney, R.A., and Fowler, W.A., 1989, Astrophys. J., 345, L5.
    !     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data 
    !       Tables, 40, 283.
    !     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.
    !     NACRE    http://pntpm3.ulb.ac.be/Nacre/

    END SUBROUTINE Rate3


    SUBROUTINE Rate4(t9)

    !-------LINKAGES.
    !     CALLED BY - [subroutine] derivs

    !-------REMARKS.
    !     Generates rate coefficients for rest of reactions.

        REAL, INTENT(IN) :: t9
        INTEGER :: reaction
        REAL :: temp

        !-------TEMPERATURE FACTORS--------------------------------

        REAL :: t913, t923, t943, t953, t912, t932, t935, t965, t938, &
                t9m13, t9m23, t9m32, t9m65, t9a, t9a13, t9a56, &
                t9b, t9b13, t9b56

        t913  = t9**(.33333333)      !t9**(1/3)
        t923  = t913*t913            !t9**(2/3)
        t943  = t923*t923            !t9**(4/3)
        t953  = t9*t923              !t9**(5/3)
        t912  = sqrt(t9)             !t9**(1/2)
        t932  = t9*t912              !t9**(3/2)
        t935  = t9**(.6)             !t9**(3/5)
        t965  = t9**(1.2)            !t9**(6/5)
        t938  = t9**(.375)           !t9**(3/8)
        t9m13 = 1.0/t913             !t9**(1/3)
        t9m23 = 1.0/t923             !t9**(-2/3)
        t9m32 = 1.0/t932             !t9**(-3/2)
        t9m65 = 1.0/t965             !t9**(-6/5)
        t9a   = t9 &                 !For reaction 82.
                   /(1.+4.78e-2*t9+7.56e-3*t953/(1.+4.78e-2*t9)**(2./3.))  
        t9a13 = t9a**(.33333333)     !t9a**(1/3)
        t9a56 = t9a**(.83333333)     !t9a**(5/6)
        t9b   = t9 &                 !For reaction 84.
                   /(1.+7.76e-2*t9+2.64e-2*t953/(1.+7.76e-2*t9)**(2./3.))
        t9b13 = t9b**(.33333333)     !t9b**(1/3)
        t9b56 = t9b**(.83333333)     !t9b**(5/6)

        !-------NEUTRON, PHOTON REACTIONS-----------------------------

        !.......C12(n,g)C13................(Wagoner 1969)
        f(65)  = 4.50e+2

        !.......C13(n,g)C14................(Wagoner 1969)
        f(66)  = 1.19e+2 + 2.38e+5*t9m32*exp(-1.67/t9)

        !.......N14(n,g)N15................(Wagoner 1969)
        f(67)  = 9.94e+3

        !-------NEUTRON, PROTON REACTIONS-----------------------------

        !.......N13(n,p)C13................(Caughlan-Fowler 1988)
        f(68)  = 1.88e+8*(1.-.167*t912+.037*t9)

        !.......N14(n,p)C14................(Caughlan-Fowler 1988)
        f(69)  = 2.39e+5*(1.+.361*t912+.502*t9) &
                  + 1.112e+8/t912*exp(-4.983/t9)

        !.......O15(n,p)N15................(Caughlan-Fowler 1988)
        f(70)  = 3.50e+8*(1.+.452*t912-.191*t9)

        !-------NEUTRON, ALPHA REACTIONS------------------------------

        !.......O15(n,a)C12................(Caughlan-Fowler 1988)
        f(71)  = 3.50e+7*(1.+.188*t912+.015*t9)

        !-------PROTON, PHOTON REACTIONS------------------------------

#if rates_1999
        !.......C12(p,g)N13................(NACRE)
        f(72)  = 2.0e+7*t9m23*exp(-13.692/t913-(t9/0.46)**2) &
                 * (1.+9.89*t9-59.8*t9**2+266.*t9**3) &
                 + 1.e+5/t932*exp(-4.913/t9)+4.24e+5/t932*exp(-21.62/t9)
#elif rates_1992
        !.......C12(p,g)N13................(Caughlan-Fowler 1988)
        f(72)  = 2.04e+7*t9m23*exp(-13.690/t913-(t9/1.500)**2) &
                  *(1.+.030*t913+1.19*t923+.254*t9+2.06*t943+1.12*t953) &
                  + 1.08e+5*t9m32*exp(-4.925/t9) &
                  + 2.15e+5*t9m32*exp(-18.179/t9)
#endif

        !.......C13(p,g)N14................(Caughlan-Fowler 1988)
        f(73)  = 8.01e+7*t9m23*exp(-13.717/t913-(t9/2.000)**2) &
                  *(1.+.030*t913+.958*t923+.204*t9+1.39*t943+.753*t953) &
                  + 1.21e+6*t9m65*exp(-5.701/t9)

        !.......C14(p,g)N15................(Caughlan-Fowler 1988)
        f(74)  = 6.80e+6*t9m23*exp(-13.741/t913-(t9/5.721)**2) &
                  *(1.+.030*t913+.503*t923+.107*t9+.213*t943+.115*t953) &
                  + 5.36e+3*t9m32*exp(-3.811/t9) &
                  + 9.82e+4*t9m13*exp(-4.739/t9)

        !.......N13(p,g)O14................(Caughlan-Fowler 1988)
        f(75)  = 4.04e+7*t9m23*exp(-15.202/t913-(t9/1.191)**2) &
                  *(1.+.027*t913-.803*t923-.154*t9+5.00*t943+2.44*t953) &
                  + 2.43e+5*t9m32*exp(-6.348/t9)

        !.......N14(p,g)O15................(Caughlan-Fowler 1988)
        f(76)  = 4.90e+7*t9m23*exp(-15.228/t913-(t9/3.294)**2) &
                  *(1.+.027*t913-.778*t923-.149*t9+.261*t943+.127*t953) &
                  + 2.37e+3*t9m32*exp(-3.011/t9) &
                  + 2.19e+4*exp(-12.530/t9)

        !.......N15(p,g)O16................(Caughlan-Fowler 1988)
        f(77)  = 9.78e+8*t9m23*exp(-15.251/t913-(t9/.450)**2) &
                  *(1.+.027*t913+.219*t923+.042*t9+6.83*t943+3.32*t953) &
                  + 1.11e+4*t9m32*exp(-3.328/t9) &
                  + 1.49e+4*t9m32*exp(-4.665/t9) &
                  + 3.80e+6*t9m32*exp(-11.048/t9)

        !-------PROTON, ALPHA REACTIONS-------------------------------

        !.......N15(p,a)C12................(Caughlan-Fowler 1988)
        f(78)  = 1.08e+12*t9m23*exp(-15.251/t913-(t9/.522)**2) &
                  *(1.+.027*t913+2.62*t923+.501*t9+5.36*t943+2.60*t953) &
                  + 1.19e+8*t9m32*exp(-3.676/t9) &
                  + 5.41e+8/t912*exp(-8.926/t9) &
                  + 4.72e+7*t9m32*exp(-7.721/t9) &
                  + 2.20e+8*t9m32*exp(-11.418/t9)

        !-------ALPHA, PHOTON REACTIONS-------------------------------

        !.......C12(a,g)O16................(Caughlan-Fowler 1988)
        f(79)  = 1.04e+8/t9**2*exp(-32.120/t913-(t9/3.496)**2) &
                  /(1.+.0489*t9m23)**2 &
                  + 1.76e+8/(t9)**2/(1.+.2654*t9m23)**2*exp(-32.120/t913) &
                  + 1.25e+3*t9m32*exp(-27.499/t9) &
                  + 1.43e-2*(t9)**5*exp(-15.541/t9)

        !-------ALPHA, PROTON REACTIONS-------------------------------

        !.......B10(a,p)C13................(Wagoner 1969)
        f(80)  = 9.60e+14*t9m23*exp(-27.99/t913)

        !.......B11(a,p)C14................(Caughlan-Fowler 1988)
        f(81)  = 5.37e+11*t9m23*exp(-28.234/t913-(t9/0.347)**2) &
                  *(1.+.015*t913+5.575*t923+.576*t9+15.888*t943+4.174*t953) &
                  + 5.44e-3*t9m32*exp(-2.827/t9) &
                  + 3.36e+2*t9m32*exp(-5.178/t9) &
                  + 5.32e+6/t938*exp(-11.617/t9)

        !.......C11(a,p)N14................(Caughlan-Fowler 1988)
        f(82)  = 7.15e+15*t9a56*t9m32*exp(-31.883/t9a13)

        !.......N12(a,p)O15................(Caughlan-Fowler 1988)
        f(83)  = 5.59e+16*t9m23*exp(-35.60/t913)

        !.......N13(a,p)O16................(Caughlan-Fowler 1988)
        f(84)  = 3.23e+17*t9b56*t9m32*exp(-35.829/t9b13)

        !-------ALPHA, NEUTRON REACTIONS------------------------------

        !.......B10(a,n)N13................(Caughlan-Fowler 1988)
        f(85)  = 1.20e+13*t9m23*exp(-27.989/t913-(t9/9.589)**2)

        !.......B11(a,n)N14................(Caughlan-Fowler 1988)
        f(86)  = 6.97e+12*t9m23*exp(-28.234/t913-(t9/0.140)**2) &
                  *(1.+.015*t913+8.115*t923+.838*t9+39.804*t943 &
                      +10.456*t953) &
                  + 1.79e+0*t9m32*exp(-2.827/t9) &
                  + 1.71e+3*t9m32*exp(-5.178/t9) &
                  + 4.49e+6*t935*exp(-8.596/t9)

        !.......B12(a,n)N15................(Wagoner 1969)
        f(87)  = 3.04e+15*t9m23*exp(-28.45/t913)

        !.......C13(a,n)O16................(Caughlan-Fowler 1988)
        f(88)  = 6.77e+15*t9m23*exp(-32.329/t913-(t9/1.284)**2) &
                  *(1.+.013*t913+2.04*t923+.184*t9) &
                  + 3.82e+5*t9m32*exp(-9.373/t9) &
                  + 1.41e+6*t9m32*exp(-11.873/t9) &
                  + 2.00e+9*t9m32*exp(-20.409/t9) &
                  + 2.92e+9*t9m32*exp(-29.283/t9)

#if linear_q_variation
        IF (BindingVariation) THEN
        DO reaction = 65, 88
            IF (q9(reaction) /= reacpr(reaction, 8)) THEN
                IF((iform(reaction).eq.3).or.(iform(reaction).eq.5).or.(iform(reaction).eq.6)) THEN
                    f(reaction) = f(reaction) * (1.0 + 0.5 * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8))
                ELSE IF((iform(reaction).eq.2).or. &
                        (iform(reaction).eq.7).or.(iform(reaction).eq.8)) THEN
                    f(reaction) = f(reaction) * (1.0 + 3.0 * (q9(reaction) - reacpr(reaction, 8))/reacpr(reaction, 8))
                END IF
            END IF
        END DO
        END IF
#else
        IF (BindingVariation) THEN
        DO reaction = 65, 88
            IF (q9(reaction) .le. 0.0) THEN
                f(reaction) = 0.0;
            ELSE IF (q9(reaction) /= reacpr(reaction, 8)) THEN
                IF((iform(reaction).eq.3).or.(iform(reaction).eq.5).or.(iform(reaction).eq.6)) THEN
                    f(reaction) = f(reaction) * sqrt(q9(reaction)/reacpr(reaction, 8))
                ELSE IF((iform(reaction).eq.2).or. &
                        (iform(reaction).eq.7).or.(iform(reaction).eq.8)) THEN
                    temp = q9(reaction)/reacpr(reaction, 8)
                    f(reaction) = f(reaction) * temp * temp * temp
                END IF
            END IF
        END DO
        END IF
#endif

        RETURN

    !-------REFERENCES--------------------------------------
    !     Caughlan, G.R., and Fowler, W.A., 1988, Atomic Data and Nuclear Data
    !       Tables, 40, 283.
    !     Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247.

    END SUBROUTINE Rate4

    SUBROUTINE CheckMaxTemps(t9)
        ! Set maxima on f(i) when t9 > t9max(i)
        REAL, INTENT(IN) :: t9
        INTEGER :: i

        DO i = 1, jsize
            if(t9 .gt. t9max(i)) then
                f(i) = fmax(i)
            end if
        END DO
    END SUBROUTINE

END MODULE reactions_module
