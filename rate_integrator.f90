MODULE rate_integrator

    use physics_parameters

    ! Reaction types and data
    integer, parameter :: nspecial = 11     ! Number of special reactions

    ! Index of reaction in standard BBN modules
    integer, save, dimension(nspecial) :: rmap

    ! reaction types:
    !    1 = Non-resonant neutron-capture
    !    2 = Non-resonant charge-induced
    !    3 = Resonant charge-induced
    integer, save, dimension(nspecial) :: Stype

    integer, parameter :: maxparam = 12     ! Maximum number of parameters for reactions
    !            |  Type
    ! parameters |   1       2       3       4
    !------------|-----------------------------------
    !          1 |   Q       Q       Q       Q
    !          2 |   mu      mu      mu      mu
    !          3 |  Emax    Emax    Emax    Emax
    !          4 |  R(0)    S(0)    S(0)     .
    !          5 |  c1/2     c1      c1      .
    !          6 |   c1      c2      c2    custom
    !          7 |  c3/2     c3      c3      .
    !          8 |   c2      c4      c4      .
    !          9 |  c5/2     c5      c5      .
    !         10 |   c3      Eg      Eg      .
    !         11 |  c7/2             Er      .
    !         12 |   c4              G       .
    ! where
    !     Q = energy released (MeV)
    !    mu = reduced mass (atomic mass units)
    !  Emax = Maximum extrapolated value of E (MeV)
    !  S(0) = Astrophysical S-factor (MeV b)
    !  c(n) = Coefficient of E^(n) in polynomial expansion
    !    Eg = Gamow factor (MeV)
    !    Er = Resonance energy (MeV)
    !     G = Resonance width (MeV)
    real, save, dimension(nspecial, maxparam) :: Sparameters
    integer, parameter :: pQ    = 1, &
                          pMu   = 2, &
                          pEmax = 3, &
                          pS0   = 4, &
                          pR0   = pS0, &
                          pC    = pS0 + 1, &
                          pEg   = 10, &
                          pEr   = 11, &
                          pG    = 12

    DATA rmap / &
        12, & !  1    p(n,g)d
        20, & !  2    d(p,g)3He
        28, & !  3    d(d,n)3He
        29, & !  4    d(d,p)t
        16, & !  5  3He(n,p)t
        30, & !  6    t(d,n)4He
        31, & !  7  3He(d,p)4He
        27, & !  8  3He(a,g)7Be
        26, & !  9    t(a,g)7Li
        17, & ! 10  7Be(n,p)7Li
        24/   ! 11  7Li(p,a)4He

    real, save, private :: SRatePrefactor
    real, save, private :: RRatePrefactor

    ! These are needed to pass reaction parameters and temperature to function
    integer, save, private :: reaction_number
    real, save, private, dimension(maxparam) :: f_params
    real, save, private :: f_kT

    ! Variation parameters
    logical, private, parameter :: shift_whole_cross_section = .false.
    real, save, private, dimension(nspecial) :: delta_Er

CONTAINS

    subroutine InitialiseRateIntegrator()
    ! Set up parameters
        integer :: i
        integer :: nuc1, nuc2
        real :: E_Gamow

        ! Prefactor for Gamow energy
        E_Gamow = 2. * pi**2 * alpha**2 * amu_MeV
        
        ! Prefactor for rates
        SRatePrefactor = NA * sqrt(8./(pi * amu_MeV * kB**3)) * c_cms * 1.e-24
        RRatePrefactor = 2./sqrt(pi * kB**3)

        Sparameters = 0.0
        !  1    p(n,g)d
        ! TODO: This reaction
        !  2    d(p,g)3He
        Stype(2) = 2
        Sparameters(2, pS0:pS0+3) = (/0.2268e-6, 22.05, 30.77, -9.919/)
        Sparameters(2, pEmax) = 1.0

        !  3    d(d,n)3He
        Stype(3) = 2
        Sparameters(3, pS0:pS0+4) = (/0.05067, 7.534, -4.225, 1.508, -0.2041/)
        Sparameters(3, pEmax) = 3.0

        !  4    d(d,p)t
        Stype(4) = 2
        Sparameters(4, pS0:pS0+2) = (/0.05115, 4.685, -1.021/)
        Sparameters(4, pEmax) = 1.5

        !  5  3He(n,p)t
        Stype(5) = 1
        Sparameters(5, pR0:pR0+8) = (/6.846e8, -0.0464743311577, -20.6566636058, &
            145.303829979, -517.845305322, 1061.59032882, -1232.39931680, &
            748.452414743, -184.417975062/)
        Sparameters(5, pEmax) = 1.0

        !  6    t(d,n)4He
        Stype(6) = 3
        Sparameters(6, pS0:pS0+5) = (/24.19, 3.453, -40.16, 285.6, -596.4, 407.1/)
        Sparameters(6, pEr) = 0.0482
        Sparameters(6, pG)  = 0.0806
        Sparameters(6, pEmax) = 0.6

        !  7  3He(d,p)4He
        Stype(7) = 3
        Sparameters(7, pS0:pS0+5) = (/18.52, -4.697, 39.53, -109.6, 130.7, -54.83/)
        Sparameters(7, pEr) = 0.183
        Sparameters(7, pG)  = 0.256
        Sparameters(7, pEmax) = 0.7

        !  8  3He(a,g)7Be
        Stype(8) = 4
        Sparameters(8, pS0:pS0+4) = (/0.3861e-3, 0.8195, -2.194, 1.419, -0.2780/)
        Sparameters(8, pEmax) = 2.0

        !  9    t(a,g)7Li
        Stype(9) = 2
        Sparameters(9, pS0:pS0+4) = (/0.08656e-3, 0.6442, -7.597, 12.16, -5.336/)
        Sparameters(9, pEmax) = 1.0

        ! 10  7Be(n,p)7Li
        Stype(10) = 4
        ! This reaction is mostly non-resonant neutron capture, so we'll use the same ones
        Sparameters(10, pR0:pR0+8) = (/4.7893e9, -4.12682044152, 3.10200988738, &
            15.8164551655, -45.5822669937, 54.7133921087, -34.7483784037, &
            11.3599443403, -1.49669812741/)
        ! TODO: Resonance here
        Sparameters(10, pEmax) = 1.0

        ! 11  7Li(p,a)4He
        Stype(11) = 2
        Sparameters(11, pS0:pS0+4) = (/0.06068, 3.174, -7.586, 8.539, -3.216/)
        Sparameters(11, pEmax) = 1.0

        do i = 1, nspecial
            ! incoming nuclei
            nuc1 = ii(rmap(i))
            nuc2 = jj(rmap(i))
            if(nuc2 .eq. 0) then    ! e.g. incoming nucleons identical
                nuc2 = nuc1
            end if

            ! Q values
            Sparameters(i, pQ) = q9(rmap(i)) * kB

            ! Reduced mass (mu, or A) in atomic mass units
            Sparameters(i, pMu) = 1. / &
                (1./(am(nuc1)+dm(nuc1)) + 1./(am(nuc2)+dm(nuc2)))

            ! Gamow factor for charge induced reactions
            if(Stype(i) .ne. 1) then
                Sparameters(i, pEg) = &
                    E_Gamow * (zm(nuc1) * zm(nuc2))**2 * Sparameters(i, pMu)
            end if
        end do
        
        ! Variation parameters
        delta_Er = 0.0
    end subroutine InitialiseRateIntegrator
    
    subroutine VaryResonance(reaction, dEr, dG)
    ! Move position of a resonance reaction (type 3)
        integer, intent(in) :: reaction
        real, intent(in) :: dEr, dG
        if(Stype(reaction) .ne. 3) then
            print *, "VaryResonance: reaction is wrong type. Reaction = ", reaction
        else if(shift_whole_cross_section) then
            delta_Er(reaction) = dEr
        else
            Sparameters(reaction, pEr) = Sparameters(reaction, pEr) + dEr
            Sparameters(reaction, pG) = Sparameters(reaction, pG) * (1. + dG)
        end if
    end subroutine

    real function GetRate(reaction, temp)
    ! Integrate function for given temp, use to calculate f
        integer, intent(IN) :: reaction
        real, intent(IN)    :: temp
        
        real :: integral
        real, dimension(maxparam) :: params

        ! Parameters for numerical integration routine
        real :: epsabs, epsrel, abserr
        integer :: neval, ierr, limit, lenw, last
        integer, parameter :: iwork_size = 1000
        integer, dimension(iwork_size) :: iwork
        real, dimension(4 * iwork_size) :: work
        
        ierr = 0
        epsabs = 0.0
        epsrel = 1.e-10
        limit = iwork_size
        lenw = 4 * iwork_size

        ! Do integration and get rate
        reaction_number = reaction
        f_params = Sparameters(reaction,:)
        f_kT = kB * temp

        select case (Stype(reaction))
        case (1)
            call dqagi(f1, 0., 1, epsabs, epsrel, integral, abserr, neval, &
                       ierr, limit, lenw, last, iwork, work)
        case (2)
            call dqagi(f2, 0., 1, epsabs, epsrel, integral, abserr, neval, &
                       ierr, limit, lenw, last, iwork, work)
        case (3)
            call dqagi(f3, 0., 1, epsabs, epsrel, integral, abserr, neval, &
                       ierr, limit, lenw, last, iwork, work)
        case (4)
            if(reaction .eq. 8) then
            call dqagi(f_reaction8, 0., 1, epsabs, epsrel, integral, abserr, neval, &
                       ierr, limit, lenw, last, iwork, work)
            end if
        end select
        if(ierr .gt. 0 .and. ierr .ne. 2) then
            print *, "RateIntegrator::GetRate: dqagi failed"
            print *, "     ierr = ", ierr,     "     t9 = ", temp
            print *, " integral = ", integral, " abserr = ", abserr
            print *, "    neval = ", neval,    "   last = ", last
            print *,  "parameters: reaction = ", reaction
            write(*,"(1p,8(e14.5))") f_params
!            stop
        end if

        if(Stype(reaction) .eq. 1) then
            GetRate = RRatePrefactor / sqrt(temp**3) * integral        
        else
            GetRate = SRatePrefactor / sqrt(Sparameters(reaction, pMu) * temp**3) * integral
        end if
    end function GetRate


    REAL FUNCTION f1(x)
    ! Non-resonant neutron capture
        real, intent(in) :: x
        
        integer :: power
        real :: xeff

        if(x .lt. f_params(pEmax)) then
            xeff = x
        else
            xeff = f_params(pEmax)
        end if
        ! Polynomial part
        f1 = 1.
        do power = 1, 8
            f1 = f1 + f_params(power + pR0) * xeff**(power/2.)
        end do
        f1 = f1 * f_params(pR0) * sqrt(x) * exp(-x/f_kT)
    END FUNCTION

    REAL FUNCTION f2(x)
    ! Non-resonant charge-induced
        real, intent(in) :: x
        
        integer :: power
        real :: xeff

        if(x .lt. f_params(pEmax)) then
            xeff = x
        else
            xeff = f_params(pEmax)
        end if
        ! Polynomial part
        f2 = 1.
        do power = 1, 5
            f2 = f2 + f_params(power + pS0) * xeff**(power)
        end do
        f2 = f2 * f_params(pS0)

        ! Gamow factor uses actual x
        f2 = f2 * exp(-sqrt(f_params(pEg)/x) - x/f_kT)
    END FUNCTION

    REAL FUNCTION f3(x)
    ! Resonant charge-induced
        real, intent(in) :: x
        
        integer :: power
        real :: xeff, xcs

        xeff = x - delta_Er(reaction_number)
        if(xeff .lt. 0.0) then
            xcs = 0.0
        else if(xeff .lt. f_params(pEmax)) then
            xcs = xeff
        else
            xcs = f_params(pEmax)
        end if

        ! Polynomial part
        f3 = 1.
        do power = 1, 5
            f3 = f3 + f_params(power + pS0) * xcs**(power)
        end do
        f3 = f3 * f_params(pS0)

        ! Gamow factor uses actual x
        f3 = f3 * exp(-sqrt(f_params(pEg)/x) - x/f_kT)
        ! Resonant part
        f3 = f3 / (1. + ((xeff - f_params(pEr))/(f_params(pG)/2.))**2)
    END FUNCTION

    REAL FUNCTION f_reaction8(x)
    ! Special function: 3He + 4He -> 7Be + gamma
        real, intent(in) :: x
        
        integer :: power
        real :: xeff, xoverEg

        ! Two dimensions are ground and excited parts
        real, dimension(2) :: total, Q, s0, s2, a
        Q  = (/1.5861, 1.1570/)
        s0 = (/0.406, 0.163/)
        s2 = (/0.007, 0.004/)
        a  = (/-0.203, -0.127/)

        if(x .lt. f_params(pEmax)) then
            xeff = x
        else
            xeff = f_params(pEmax)
        end if
        
        xoverEg = xeff/f_params(pEg)

        total = Q / (xeff + Q)
        total = total * (s0 * (1. + a * xeff)**2 &
                         + s2 * (1. + 4. * pi**2 * xoverEg) * (1. + 16. * pi**2 * xoverEg))

        f_reaction8 = (total(1) + total(2))/1000.

        ! Gamow factor uses actual x
        f_reaction8 = f_reaction8 * exp(-sqrt(f_params(pEg)/x) - x/f_kT)
    END FUNCTION

    REAL FUNCTION f_reaction10(x)
    ! Special function: 7Be + n -> 7Li + p
        real, intent(in) :: x
        
        integer :: power
        real :: xeff, xcs

        ! This reaction is mostly non-resonant neutron capture
        real, dimension(9) :: NonResParams = &
            (/4.7893e9, -4.12682044152, 3.10200988738, &
            15.8164551655, -45.5822669937, 54.7133921087, -34.7483784037, &
            11.3599443403, -1.49669812741/)

        ! Resonance parameters
        real, dimension(2) :: ResPositions = (/0.32, 2.7/)
        real, dimension(2) :: ResWidths    = (/0.20, 1.9/)
        real, dimension(2) :: ResHeights   = (/1.0553e9, 2.0364e9/)

        xeff = x - delta_Er(10)
        if(xeff .lt. 0.0) then
            xcs = 0.0
        else if(xeff .lt. f_params(pEmax)) then
            xcs = xeff
        else
            xcs = f_params(pEmax)
        end if

        ! Non-resonant polynomial part
        f_reaction10 = 1.
        do power = 1, 8
            f_reaction10 = f_reaction10 + NonResParams(power + 1) * xcs**(power/2.)
        end do        
        f_reaction10 = f_reaction10 * NonResParams(1)

        ! Resonances
        f_reaction10 = f_reaction10 + ResHeights(1)/(1. + ((xeff - ResPositions(1))/(ResWidths(1)/2.))**2)
        f_reaction10 = f_reaction10 + ResHeights(2)/(1. + ((xeff - ResPositions(2))/(ResWidths(2)/2.))**2)

        f_reaction10 = f_reaction10 * sqrt(x) * exp(-x/f_kT)
    END FUNCTION


    SUBROUTINE Trapezoidal(func, integral)
        real, external :: func
        real, intent(out) :: integral

        integer :: i, imax
        real :: x, dx, total

        x = 0.
        dx = 0.01
        imax = 3./dx
        
        total = 0.
        do i=1, imax
            x = x + dx
            total = total + dx * func(x)
        end do
        integral = total
    END SUBROUTINE

END MODULE rate_integrator
