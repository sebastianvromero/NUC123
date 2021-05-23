MODULE physics_parameters

    REAL, SAVE :: pi
    REAL, SAVE :: alpha, &  ! Fine-structure constant
                  kB, &     ! Boltzmann's constant (MeV / 10^9 K) = 1/11.604505
                  NA, &     ! Avagadro's number
                  c_cms, &  ! Speed of light (cm/s)
                  amu_MeV   ! Atomic mass unit (MeV)

    !Early universe model parameters
    REAL, SAVE :: g                     !Gravitational constant.
    REAL, SAVE :: tau                   !Neutron lifetime (sec).
    REAL, SAVE :: xnu                   !Number of neutrino species.
    REAL, SAVE, DIMENSION(3) :: c
                        !c(1) is variation of gravitational constant.
                        !c(2) is neutron half-life (min).
                        !c(3) is number of neutrino species.
    REAL, SAVE :: cosmo     !Cosmological constant.
    REAL, SAVE, DIMENSION(3) :: xi
                        !xi(1) is e neutrino degeneracy parameter.
                        !xi(2) is m neutrino degeneracy parameter.
                        !xi(3) is t neutrino degeneracy parameter.

    REAL, SAVE :: eta   !Baryon-to-photon ratio.

    !Neutrino parameters.
    REAL, SAVE :: t9mev, &  !Temperature (in units of MeV).
            tnmev, &        !Neutrino temperature (in units of MeV).
            tnu, &          !Neutrino temperature.
            cnorm, &        !Normalizing constant.
            rhonu           !Neutrino energy density.
    INTEGER, SAVE :: nu     !Type of neutrino.

    INTEGER, PARAMETER :: nnuc = 26     !Number of nuclides in calculation.

    !Nuclide data
    REAL, SAVE, DIMENSION(nnuc) :: &
            am, &       !Atomic number of nuclide.
            zm, &       !Charge of nuclide.
            dm, &       !Mass excess of nuclide (amu)
            km          !Mass excess of nuclide (keV)

    INTEGER, PARAMETER :: nrec = 88     !Number of nuclear reactions.

    !Reaction rate parameters
    INTEGER, SAVE, DIMENSION(nrec) :: &
            iform, &        !Reaction type code (1-11).
            ii, &           !Incoming nuclide type (1-26).
            jj, &           !Incoming light nuclide type (1-6).
            kk, &           !Outgoing light nuclide type (1-6).
            ll              !Outgoing nuclide type (1-26).
    REAL, SAVE, DIMENSION(nrec) :: &
            rev, &          !Reverse reaction coefficient.
            q9              !Energy released in reaction (in 10**9 K).
    REAL, SAVE, DIMENSION(nrec) :: &
            EgIn, &         !Gamow energy of entrance channel (in 10**9 K).
            EgOut           !Gamow energy of exit channel (in 10**9 K).

    !Initial reaction parameter values.
    REAL, SAVE, DIMENSION(nrec, 8) :: reacpr

    ! Number of nuclides in reaction types 1-11
    REAL, SAVE, DIMENSION(11) :: si, sj, sk, sl

!=====================DATA DIVISION========================

!    Nuclide and corresponding number
!    --------------------
!    1) N         7) Li6      13) B10      19) C13      25) O15
!    2) P         8) Li7      14) B11      20) N13      26) O16
!    3) H2        9) Be7      15) C11      21) C14
!    4) H3       10) Li8      16) B12      22) N14
!    5) He3      11) B8       17) C12      23) O14
!    6) He4      12) Be9      18) N12      24) N15

!--------NUCLIDE DATA.
      DATA am /1.,1.,2.,3.,3.,4.,6.,7.,7.,8.,8.,9.,10.,11.,11.,12., &
               12.,12.,13.,13.,14.,14.,14.,15.,15.,16./
      DATA zm /0.,1.,1.,1.,2.,2.,3.,3.,4.,3.,5.,4.,5.,5.,6.,5., &
               6.,7.,6.,7.,6.,7.,8.,7.,8.,8./
      DATA dm /.008665,.007825,.014102,.016050,.016030,.002603,.015125, &
               .016004,.016929,.022487,.024609,.012186,.012939,.009305, &
               .011432,.014354,.000000,.018641,.003354,.005738,.003242, &
               .003074,.008597,.000108,.003070,-.005085/
      DATA km /8071.31710,7288.97050,13135.72158,14949.80600,14931.21475,2424.91565, &
               14086.793,14908.141,15770.034,20946.844,22921.490,11347.648, &
               12050.731,8667.931,10650.342,13368.899,0.0,17338.082, &
               3125.01129,5345.481,3019.89305,2863.41704,8007.356,101.43805, &
               2855.605,-4737.00141/

      ! Number of nuclides in reaction types 1-11
      DATA si /1.,1.,1.,1.,1.,2.,3.,2.,1.,1.,2./
      DATA sj /0.,1.,1.,0.,1.,0.,0.,1.,1.,1.,0./
      DATA sk /0.,0.,1.,0.,0.,1.,0.,0.,1.,0.,2./
      DATA sl /1.,1.,1.,2.,2.,1.,1.,1.,2.,3.,1./

! Select q-values to use
#define rates_1999 1
#define rates_1992 0

      INTEGER, PRIVATE, SAVE :: i, j
      DATA ((reacpr(i,j),j=1,8),i=1,11) / &
      !        reac# type n1 n2 n3 n4 rev-coeff q-value
      !        ----  ---- -- -- -- -- ------ -------
                   1.,1., 1.,0.,0., 2., 0.0  ,   0.0 , &   !N->P         
                   2.,1., 4.,0.,0., 5., 0.0  ,   0.0 , &   !H3->He3
                   3.,4.,10.,0.,0., 6., 0.0  ,   0.0 , &   !Li8->2He4
                   4.,1.,16.,0.,0.,17., 0.0  ,   0.0 , &   !B12->C12
                   5.,1.,21.,0.,0.,22., 0.0  ,   0.0 , &   !C14->N14
                   6.,4.,11.,0.,0., 6., 0.0  ,   0.0 , &   !B8->2He4
                   7.,1.,15.,0.,0.,14., 0.0  ,   0.0 , &   !C11->B11
                   8.,1.,18.,0.,0.,17., 0.0  ,   0.0 , &   !N12->C12
                   9.,1.,20.,0.,0.,19., 0.0  ,   0.0 , &   !N13->C13
                  10.,1.,23.,0.,0.,22., 0.0  ,   0.0 , &   !O14->N14
                  11.,1.,25.,0.,0.,24., 0.0  ,   0.0 /     !O15->N15
      DATA ((reacpr(i,j),j=1,8),i=12,22) / &
      !        reac# type n1 n2 n3 n4 rev-coeff q-value
      !        ----  ---- -- -- -- -- ------ -------
                  12.,2., 2.,1.,0., 3., 0.471,  25.82, &   !H(n,g)H2
                  13.,2., 3.,1.,0., 4., 1.63 ,  72.62, &   !H2(n,g)H3
                  14.,2., 5.,1.,0., 6., 2.61 , 238.81, &   !He3(n,g)He4
                  15.,2., 7.,1.,0., 8., 1.19 ,  84.17, &   !Li6(n,g)Li7
                  16.,3., 5.,1.,2., 4., 1.002,   8.863, &  !He3(n,p)H3
                  17.,3., 9.,1.,2., 8., 0.998,  19.081, &  !Be7(n,p)Li7
                  18.,3., 7.,1.,4., 6., 1.070,  55.494, &  !Li6(n,a)H3
                  19.,5., 9.,1.,0., 6., 4.70 , 220.39, &   !Be7(n,a)He4
#if rates_1999
                  20.,2., 3.,2.,0., 5., 1.63 ,  63.749, &  !H2(p,g)He3
                  21.,2., 4.,2.,0., 6., 2.61 , 229.932, &  !H3(p,g)He4  (unchanged)
                  22.,2., 7.,2.,0., 9., 1.187,  65.052/    !Li6(p,g)Be7
#elif rates_1992
                  20.,2., 3.,2.,0., 5., 1.63 ,  63.750, &  !H2(p,g)He3
                  21.,2., 4.,2.,0., 6., 2.61 , 229.932, &  !H3(p,g)He4
                  22.,2., 7.,2.,0., 9., 1.19 ,  65.054/    !Li6(p,g)Be7
#endif
      DATA ((reacpr(i,j),j=1,8),i=23,33) / &
      !        reac# type n1 n2 n3 n4 rev-coeff q-value
      !        ----  ---- -- -- -- -- ------ -------
#if rates_1999
                  23.,3., 7.,2.,5., 6., 1.067,  46.641, &  !Li6(p,a)He3
                  24.,5., 8.,2.,0., 6., 4.676, 201.30 , &  !Li7(p,a)He4
                  25.,2., 6.,3.,0., 7., 1.53 ,  17.118, &  !H2(a,g)Li6  (unchanged)
                  26.,2., 6.,4.,0., 8., 1.11 ,  28.640, &  !H3(a,g)Li7  (unchanged)
                  27.,2., 6.,5.,0., 9., 1.11 ,  18.423, &  !He3(a,g)Be7 (unchanged)
                  28.,6., 3.,0.,1., 5., 1.73 ,  37.935, &  !H2(d,n)He3  (unchanged)
                  29.,6., 3.,0.,2., 4., 1.732,  46.797, &  !H2(d,p)H3
                  30.,3., 4.,3.,1., 6., 5.506, 204.12,  &  !H3(d,n)He4
                  31.,3., 5.,3.,2., 6., 5.55 , 212.980, &  !He3(d,p)He4 (unchanged)
                  32.,11.,5.,0.,2., 6., 3.392, 149.230, &  !He3(He3,2p)He4
#elif rates_1992
                  23.,3., 7.,2.,5., 6., 1.07 ,  46.631, &  !Li6(p,a)He3
                  24.,5., 8.,2.,0., 6., 4.69 , 201.291, &  !Li7(p,a)He4
                  25.,2., 6.,3.,0., 7., 1.53 ,  17.118, &  !H2(a,g)Li6
                  26.,2., 6.,4.,0., 8., 1.11 ,  28.640, &  !H3(a,g)Li7
                  27.,2., 6.,5.,0., 9., 1.11 ,  18.423, &  !He3(a,g)Be7
                  28.,6., 3.,0.,1., 5., 1.73 ,  37.935, &  !H2(d,n)He3
                  29.,6., 3.,0.,2., 4., 1.73 ,  46.798, &  !H2(d,p)H3
                  30.,3., 4.,3.,1., 6., 5.54 , 204.117, &  !H3(d,n)He4
                  31.,3., 5.,3.,2., 6., 5.55 , 212.980, &  !He3(d,p)He4
                  32.,11.,5.,0.,2., 6., 3.39 , 149.230, &  !He3(He3,2p)He4
#endif
                  33.,9., 8.,3.,1., 6., 9.95 , 175.476/    !Li7(d,na)He4
      DATA ((reacpr(i,j),j=1,8),i=34,44) / &
      !        reac# type n1 n2 n3 n4 rev-coeff q-value
      !        ----  ---- -- -- -- -- ------ -------
                  34.,9., 9.,3.,2., 6., 9.97 , 194.557, &  !Be7(d,pa)He4
                  35.,2., 8.,1.,0.,10., 1.31 ,  23.59, &   !Li7(n,g)Li8
                  36.,2.,13.,1.,0.,14., 3.04 , 132.95, &   !B10(n,g)B11
                  37.,2.,14.,1.,0.,16., 2.34 ,  39.10, &   !B11(n,g)B12
                  38.,3.,15.,1.,2.,14., 1.002,  32.080, &  !C11(n,p)B11
                  39.,3.,13.,1.,6., 8., 0.758,  32.382, &  !B10(n,a)Li7
                  40.,2., 9.,2.,0.,11., 1.30 ,   1.595, &  !Be7(p,g)B8
#if rates_1999
                  41.,2.,12.,2.,0.,13., 0.9734, 76.424, &  !Be9(p,g)B10
                  42.,2.,13.,2.,0.,15., 3.026 , 100.83, &  !B10(p,g)C11
                  43.,2.,14.,2.,0.,17., 7.004 , 185.17, &  !B11(p,g)C12
#elif rates_1992
                  41.,2.,12.,2.,0.,13., 0.973,  76.427, &  !Be9(p,g)B10
                  42.,2.,13.,2.,0.,15., 3.03 , 100.840, &  !B10(p,g)C11
                  43.,2.,14.,2.,0.,17., 7.01 , 185.173, &  !B11(p,g)C12
#endif
                  44.,2.,15.,2.,0.,18., 2.33 ,   6.975/    !C11(p,g)N12
      DATA ((reacpr(i,j),j=1,8),i=45,55) / &
      !        reac# type n1 n2 n3 n4 rev-coeff q-value
      !        ----  ---- -- -- -- -- ------ -------
                  45.,3.,16.,2.,1.,17., 3.00 , 146.08, &   !B12(p,n)C12
#if rates_1999
                  46.,3.,12.,2.,6., 7., 0.6177, 24.663, &  !Be9(p,a)Li6
                  47.,3.,13.,2.,6., 9., 0.7537, 13.291, &  !B10(p,a)Be7
                  48.,3.,16.,2.,6.,12., 0.292,  79.89, &   !B12(p,a)Be9 (unchanged)
                  49.,2., 7.,6.,0.,13., 1.58 ,  51.753, &  !Li6(a,g)B10 (unchanged)
                  50.,2., 8.,6.,0.,14., 4.015, 100.55, &   !Li7(a,g)B11
                  51.,2., 9.,6.,0.,15., 4.016,  87.541, &  !Be7(a,g)C11
                  52.,3.,11.,6.,2.,15., 3.08 ,  86.00, &   !B8(a,p)C11  (unchanged)
                  53.,3.,10.,6.,1.,14., 3.07 ,  76.96, &   !Li8(a,n)B11 (unchanged)
                  54.,3.,12.,6.,1.,17.,10.28 ,  66.158, &  !Be9(a,n)C12
#elif rates_1992
                  46.,3.,12.,2.,6., 7., 0.618,  24.674, &  !Be9(p,a)Li6
                  47.,3.,13.,2.,6., 9., 0.754,  13.301, &  !B10(p,a)Be7
                  48.,3.,16.,2.,6.,12., 0.292,  79.89, &   !B12(p,a)Be9
                  49.,2., 7.,6.,0.,13., 1.58 ,  51.753, &  !Li6(a,g)B10
                  50.,2., 8.,6.,0.,14., 4.02 , 100.538, &  !Li7(a,g)B11
                  51.,2., 9.,6.,0.,15., 4.02 ,  87.539, &  !Be7(a,g)C11
                  52.,3.,11.,6.,2.,15., 3.08 ,  86.00, &   !B8(a,p)C11
                  53.,3.,10.,6.,1.,14., 3.07 ,  76.96, &   !Li8(a,n)B11   
                  54.,3.,12.,6.,1.,17.,10.3  ,  66.160, &  !Be9(a,n)C12
#endif
                  55.,3.,12.,3.,1.,13., 2.07 ,  50.63/     !Be9(d,n)B10
      DATA ((reacpr(i,j),j=1,8),i=56,66) / &
      !        reac# type n1 n2 n3 n4 rev-coeff q-value
      !        ----  ---- -- -- -- -- ------ -------
                  56.,3.,13.,3.,2.,14., 6.44 , 107.13, &   !B10(d,p)B11
                  57.,3.,14.,3.,1.,17.,14.9  , 159.36, &   !B11(d,n)C12
#if rates_1999
                  58.,8., 6.,1.,0.,12., 0.5844, 18.258, &  !He4(an,g)Be9
                  59.,7., 6.,0.,0.,17., 2.003 , 84.415, &  !He4(2a,g)C12
                  60.,9.,10.,2.,1., 6., 3.58 , 177.73, &   !Li8(p,na)He4 (unchanged)
                  61.,9.,11.,1.,2., 6., 3.58 , 218.82, &   !B8(n,pa)He4  (unchanged)
                  62.,9.,12.,2.,3., 6., 0.8071,  7.556, &  !Be9(p,da)He4
                  63.,10.,14.,2.,0.,6., 3.501 , 100.76, &  !B11(p,2a)He4
#elif rates_1992
                  58.,8., 6.,1.,0.,12., 0.584,  18.260, &  !He4(an,g)Be9
                  59.,7., 6.,0.,0.,17., 2.00 ,  84.420, &  !He4(2a,g)C12
                  60.,9.,10.,2.,1., 6., 3.58 , 177.73, &   !Li8(p,na)He4
                  61.,9.,11.,1.,2., 6., 3.58 , 218.82, &   !B8(n,pa)He4
                  62.,9.,12.,2.,3., 6., 0.807,   7.555, &  !Be9(p,da)He4
                  63.,10.,14.,2.,0.,6., 3.50 , 100.753, &  !B11(p,2a)Be4
#endif
                  64.,10.,15.,1.,0.,6., 3.49 , 132.83, &   !C11(n,2a)He4
                  65.,2.,17.,1.,0.,19., 0.886,  57.41, &   !C12(n,g)C13
                  66.,2.,19.,1.,0.,21., 3.58 ,  94.88/     !C13(n,g)C14
      DATA ((reacpr(i,j),j=1,8),i=67,77) / &
      !        reac# type n1 n2 n3 n4 rev-coeff q-value
      !        ----  ---- -- -- -- -- ------ -------
                  67.,2.,22.,1.,0.,24., 2.71 , 125.74, &   !N14(n,g)N15
                  68.,3.,20.,1.,2.,19., 1.002,  34.846, &  !N13(n,p)C13
                  69.,3.,22.,1.,2.,21., 3.003,   7.263, &  !N14(n,p)C14
                  70.,3.,25.,1.,2.,24., 1.002,  41.037, &  !O15(n,p)N15
                  71.,3.,25.,1.,6.,17., 0.709,  98.661, &  !O15(n,a)C12
#if rates_1999
                  72.,2.,17.,2.,0.,20., 0.8847, 22.553, &  !C12(p,g)N13
                  73.,2.,19.,2.,0.,22., 1.19 ,  87.619, &  !C13(p,g)N14
                  74.,2.,21.,2.,0.,24., 0.900, 118.452, &  !C14(p,g)N15 (unchanged)
                  75.,2.,20.,2.,0.,23., 3.571,  53.705, &  !N13(p,g)O14
                  76.,2.,22.,2.,0.,25., 2.699,  84.677, &  !N14(p,g)O15
                  77.,2.,24.,2.,0.,26., 3.622, 140.73/     !N15(p,g)O16
#elif rates_1992
                  72.,2.,17.,2.,0.,20., 0.884,  22.553, &  !C12(p,g)N13
                  73.,2.,19.,2.,0.,22., 1.19 ,  87.621, &  !C13(p,g)N14
                  74.,2.,21.,2.,0.,24., 0.900, 118.452, &  !C14(p,g)N15
                  75.,2.,20.,2.,0.,23., 3.57 ,  53.706, &  !N13(p,g)O14
                  76.,2.,22.,2.,0.,25., 2.70 ,  84.678, &  !N14(p,g)O15
                  77.,2.,24.,2.,0.,26., 3.62 , 140.734/    !N15(p,g)O16
#endif
      DATA ((reacpr(i,j),j=1,8),i=78,88) / &
      !        reac# type n1 n2 n3 n4 rev-coeff q-value
      !        ----  ---- -- -- -- -- ------ -------
                  78.,3.,24.,2.,6.,17., 0.706,  57.623, &  !N15(p,a)C12
                  79.,2.,17.,6.,0.,26., 5.13 ,  83.111, &  !C12(a,g)O16
                  80.,3.,13.,6.,2.,19., 9.36 ,  47.16, &   !B10(a,p)C13  
                  81.,3.,14.,6.,2.,21.,11.0  ,   9.098, &  !B11(a,p)C14  
                  82.,3.,15.,6.,2.,22., 3.68 ,  33.915, &  !C11(a,p)N14  
                  83.,3.,18.,6.,2.,25., 4.26 , 111.87, &   !N12(a,p)O15  
                  84.,3.,20.,6.,2.,26., 5.81 ,  60.557, &  !N13(a,p)O16  
                  85.,3.,13.,6.,1.,20., 9.34 ,  12.287, &  !B10(a,n)N13  
                  86.,3.,14.,6.,1.,22., 3.67 ,   1.835, &  !B11(a,n)N14  
                  87.,3.,16.,6.,1.,24., 4.25 ,  88.47, &   !B12(a,n)N15  
                  88.,3.,19.,6.,1.,26., 5.79 ,  25.711/    !C13(a,n)O16  

!=====================SUBROUTINES DIVISION=================

CONTAINS

    SUBROUTINE ResetPhysicsParameters

        pi = 2.0 * acos(0.0)
        alpha = 1./137.035999679
        kB = 0.08617343         ! Boltzmann constant (MeV / 10^9 K)
                                !  = 1.0/11.604505
        NA = 6.02214179e23      ! Avagadro's number
        c_cms = 29979245800.    ! Speed of light (cm/s)
        amu_MeV = 931.494028    ! Atomic mass unit (MeV)

        ! Set model parameters to default values
        c(1)  = 1.00            !Variation of gravitational constant.
        c(2)  = 885.7           !Neutron lifetime.
        c(3)  = 3.0             !Number of neutrino species.
        cosmo = 0.0             !Cosmological constant.
        xi(1) = 0.0             !Electron degeneracy parameter.
        xi(2) = 0.0             !Muon degeneray parameter.
        xi(3) = 0.0             !Tauon degeneracy parameter.
        eta   = 6.23e-10        !Baryon-to-photon ratio.
        !   BBN 1992 value: eta = 3.162e-10
        !   WMAP 3 value:   eta = 6.14(25)e-10
        !   WMAP 5 value:   eta = 6.23(17)e-10

    END SUBROUTINE ResetPhysicsParameters


END MODULE physics_parameters
