MODULE universe_module

    use variables
    use scalarfield_module

    TYPE Universe
        !Evolution parameters and derivatives
        REAL :: t9,   dt9,  &               !Temperature of photons (units of 10**9 K).
                hv,   dhv,  &               !Defined by hv = M(atomic)n(baryon)/t9**3.
                phie, dphie                 !Chemical potential of electron.
        REAL, DIMENSION(nnuc) :: y, dydt    !Relative number abundances.

        !Scalar field and time derivatives.
        TYPE(Scalar) scalarfield
    END TYPE Universe

    !Evolution parameters at end of 1st R-K loop.
    TYPE(Universe), SAVE :: uni, uni0
    
CONTAINS

    SUBROUTINE IncrementUniverse(v, v0, v1, dt, ymin)
        ! v = v0 + dv1 * dt
        ! minimum values of y is ymin
        TYPE(Universe), INTENT(OUT) :: v
        TYPE(Universe), INTENT(IN) :: v0, v1
        REAL, INTENT(IN) :: dt, ymin
        
        INTEGER :: i

        v%t9   = v0%t9   + v1%dt9 * dt
        v%hv   = v0%hv   + v1%dhv * dt
        v%phie = v0%phie + v1%dphie * dt
 
        v%y = v0%y + v1%dydt * dt
        DO i = 1, nnuc
            IF (v%y(i).lt.ymin) v%y(i) = ymin
        END DO

        v%scalarfield%scal = v0%scalarfield%scal + v1%scalarfield%dscal * dt
        v%scalarfield%dscal = v0%scalarfield%dscal + v1%scalarfield%d2scal * dt
    END SUBROUTINE IncrementUniverse

    SUBROUTINE StepRK2(v, v0, dt, ymin)
        ! v = v0 + 0.5 (dv + dv0) * dt
        ! minimum values of y is ymin
        TYPE(Universe), INTENT(INOUT) :: v
        TYPE(Universe), INTENT(IN) :: v0
        REAL, INTENT(IN) :: dt, ymin

        INTEGER :: i

        v%t9 = v0%t9 + 0.5 * (v%dt9 + v0%dt9) * dt
        v%hv = v0%hv + 0.5 * (v%dhv + v0%dhv) * dt
        v%phie = v0%phie + 0.5 * (v%dphie + v0%dphie) * dt

        v%y = v0%y + 0.5 * (v%dydt + v0%dydt) * dt
        DO i = 1, nnuc
            IF (v%y(i).lt.ymin) v%y(i) = ymin  !Set at minimum value.
        END DO
        
        ! Inject additional element
        !CALL InjectElement(v, dt)

        v%scalarfield%scal  = v0%scalarfield%scal &
                                 + .5*(v%scalarfield%dscal + v0%scalarfield%dscal)*dt
        v%scalarfield%dscal = v0%scalarfield%dscal &
                                 + .5*(v%scalarfield%d2scal + v0%scalarfield%d2scal)*dt
    END SUBROUTINE StepRK2

    SUBROUTINE StepRK4(v, v0, v1, v2, v3, dt, ymin)
        ! v = v1 + (1/6*dv0 + 1/3*dv1 + 1/3*dv2 + 1/6*dv3) * dt
        ! minimum values of y is ymin
        TYPE(Universe), INTENT(OUT) :: v
        TYPE(Universe), INTENT(IN) :: v0, v1, v2, v3
        REAL, INTENT(IN) :: dt, ymin

        INTEGER :: i

        v%t9 = v0%t9 + (v0%dt9 + 2.0*v1%dt9 + 2.0*v2%dt9 + v3%dt9)/6.0 * dt
        v%hv = v0%hv + (v0%dhv + 2.0*v1%dhv + 2.0*v2%dhv + v3%dhv)/6.0 * dt
        v%phie = v0%phie + (v0%dphie + 2.0*v1%dphie + 2.0*v2%dphie + v3%dphie)/6.0 * dt

        v%y = v0%y + (v0%dydt + 2.0*v1%dydt + 2.0*v2%dydt + v3%dydt)/6.0 * dt
        DO i = 1, nnuc
            IF (v%y(i).lt.ymin) v%y(i) = ymin  !Set at minimum value.
        END DO

        v%scalarfield%scal  = v0%scalarfield%scal &
              + (v0%scalarfield%dscal + 2.0*v1%scalarfield%dscal + &
                 2.0*v2%scalarfield%dscal + v3%scalarfield%dscal)/6.0 * dt
        v%scalarfield%dscal = v0%scalarfield%dscal &
              + (v0%scalarfield%d2scal + 2.0*v1%scalarfield%d2scal + &
                 2.0*v2%scalarfield%d2scal + v3%scalarfield%d2scal)/6.0 * dt
    END SUBROUTINE StepRK4


    SUBROUTINE InjectElement(uni, dt)
        ! Add element to universe at particular temperatures
        use variables, only : dlt9dt

        TYPE(Universe), INTENT(INOUT) :: uni
        REAL, intent(IN) :: dt
        real :: highT9, lowT9, amount
        integer :: element
        real :: dt9, amount_to_add
        real, save :: amount_added = 0.0

        element = 3
        amount  = 1.e-3
        highT9  = 0.9
        lowT9   = 0.8

        if((highT9 > uni%t9) .and. (uni%t9 > lowT9)) then
            dt9 = uni%t9 * dlt9dt * dt
            amount_to_add = amount * dt9/(lowT9 - highT9)
            uni%y(element) = uni%y(element) + amount_to_add
            amount_added = amount_added + amount_to_add
!            print *, dt9, amount_added
        end if
    END SUBROUTINE

    
END MODULE universe_module