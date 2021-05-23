MODULE bessel_module

    !Evaluation of functions bl, bm, bn
    REAL, SAVE, DIMENSION(5) :: &
            bl, &       !Evaluation of function bl(z).
            bm, &       !Evaluation of function bm(z).
            bn          !Evaluation of function bn(z).

    !Evaluation of modified bessel functions (K)
    REAL, SAVE, PRIVATE :: bk0,bk1,bk2,bk3,bk4  !Values k0(r),k1(r),k2(r),k3(r),k4(r).
    
    PRIVATE knux

CONTAINS

    SUBROUTINE bessel(z)

!-------LINKAGES.
!     CALLED BY - [subroutine] start, therm
!     CALLS     - [subroutine] knux

!-------REMARKS.
!     Evaluates functions bl(z), bm(z), and bn(z) using solutions to
!     modified Bessel functions.

!=================DECLARATION DIVISION=====================

!-------LOCAL VARIABLES.
      REAL, INTENT(IN) :: z        !Defined by z = m(electron)*c**2/k*t9.
      REAL :: r                    !Multiples of z.
      INTEGER :: i

!==================PROCEDURE DIVISION======================

      DO i=1,5
        r=i*z                      !Multiples of z.
        CALL knux(r)               !Get k0(r),k1(r),k2(r),k3(r),k4(r),k(5).
        bl(i) = bk2/r                   !Function bl.
        bm(i) = 0.25*(3.*bk3+bk1)/r     !Function bm.
        bn(i) = 0.5*(bk4+bk2)/r         !Function bn.
      END DO
      RETURN

    END SUBROUTINE bessel



!===============IDENTIFICATION DIVISION====================

    SUBROUTINE knux(z)

!-------LINKAGES.
!     CALLED BY - [subroutine] bessel
!     CALLS     - [function] exp

!-------REMARKS.
!     A subroutine for modified bessel functions of the second kind
!     k-nu(z).

      REAL, INTENT(IN) :: z                    !Input variable.

!==================DECLARATION DIVISION====================

!--------MODIFIED BESSEL FUNCTION VALUES.
      REAL    bi0,bi1              !Values i0(z),i1(z).

!--------EXPANSION COEFFICIENTS.
      REAL    ci0(7)               !Expansion coefficients for i0 (z.le.2).
      REAL    ci1(7)               !Expansion coefficients for i1 (z.le.2).
      REAL    ck0(7)               !Expansion coefficients for k0 (z.le.2).
      REAL    ck1(7)               !Expansion coefficients for k1 (z.le.2).
      REAL    c0(7)                !Expansion coefficients for k0 (z.gt.2).
      REAL    c1(7)                !Expansion coefficients for k1 (z.gt.2).

!--------VARIABLES TO BE EVALUATED.
      REAL    y                    !Expansion variable = z/2.
      REAL    t                    !Expansion variable = z/3.75.
      REAL    coeff                !Logrithmic or exponential coefficient.
      
      INTEGER i

!=====================DATA DIVISION========================

!-------EXPANSION COEFFICIENTS.
      DATA ci0 / 1., &
                 3.5156229,      3.0899424,      1.2067492, &
                 0.2659732,      0.0360768,      0.0045813/
      DATA ci1 / 0.5, &
                 0.87890594,     0.51498869,     0.15084934, &
                 0.02658733,     0.00301532,     0.00032411/
      DATA ck0 /-0.57721566, &
                 0.42278420,     0.23069756,     0.03488590, &
                 0.00262698,     0.00010750,     0.00000740/
      DATA ck1 / 1., &
                 0.15443144,    -0.67278579,    -0.18156897, &
                -0.01919402,    -0.00110404,    -0.00004686/
      DATA c0  / 1.25331414, &
                -0.07832358,     0.02189568,    -0.01062446, &
                 0.00587872,    -0.00251540,     0.00053208/
      DATA c1  / 1.25331414, &
                 0.23498619,    -0.03655620,     0.01504268, &
                -0.00780353,     0.00325614,    -0.00068245/


!==================PROCEDURE DIVISION======================

!10-----COMPUTE K0 AND K1----------------------------------

      IF (z.le.2.) THEN            !(Ref. 1).
!..........COMPUTE FACTORS.
        t = (z/3.75)
        y = (z/2)
        coeff = alog(y)
!..........VALUES FOR i0(z) and i1(z).
        bi0 = ci0(1)
        bi1 = ci1(1)
        bk0 = ck0(1)
        bk1 = ck1(1)
        DO i = 2,7
          bi0 = bi0 + ci0(i)*t**(2*(i-1))
          bi1 = bi1 + ci1(i)*t**(2*(i-1))
          bk0 = bk0 + ck0(i)*y**(2*(i-1))
          bk1 = bk1 + ck1(i)*y**(2*(i-1))
        END DO
!..........VALUES FOR k0(z) and k1(z).
        bk0 = -coeff*bi0 + bk0
        bk1 = coeff*bi1*z + bk1/z
      ELSE !(z.le.2.)               !(Ref. 2).
!..........COMPUTE FACTORS.
        y = (2.0/z)
        coeff = (exp(-z)/sqrt(z))
!..........VALUES FOR k0(z) and k1(z).
        bk0 = c0(1)
        bk1 = c1(1)       
        DO i = 2,7
          bk0 = bk0 + c0(i)*y**(i-1)
          bk1 = bk1 + c1(i)*y**(i-1)
        END DO
        bk0 = coeff*bk0
        bk1 = coeff*bk1
      END IF !(z.le.2.) 

!20-----FIND K2, K3, AND K4 BY ITERATION (Ref. 3)-------------------

      bk2 = 2.*(bk1/z) + bk0       !k2(z).
      bk3 = 4.*(bk2/z) + bk1       !k3(z).
      bk4 = 6.*(bk3/z) + bk2       !k4(z).

      RETURN

!-------REFERENCES--------------------------------------
!     Handbook of Mathematical Functions (Abramowitz and Stegun),
!       Dover Publications, Inc., New York
!       1) Polynomial approximations for z.le.2
!         page 378, equations 9.8.1 and 9.8.3.
!         page 379, equations 9.8.5 and 9.8.7.
!       2) Polynomial approximations for z > 2
!         page 379, equations 9.8.6 and 9.8.8.
!       3) Recursion relation from 1st line of 9.6.26, page 376.

      END SUBROUTINE knux

END MODULE bessel_module