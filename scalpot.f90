MODULE scalarfield_module

    ! Potential parameters
    REAL, SAVE :: AA, Phi0, lambda

    ! Scalar field
    TYPE Scalar
        REAL :: scal, dscal, d2scal     ! Scalar field value, d(scal)/dt, d2(scal)/dt2
    END TYPE Scalar

CONTAINS

    SUBROUTINE Initialise(phi, scal0, dscal0, d2scal0)
        ! Initialise the scalar field. Default values all zero.
        TYPE(Scalar), INTENT(INOUT) :: phi
        REAL, OPTIONAL :: scal0, dscal0, d2scal0

        ! Initialise potential parameters
        AA      = -1.0
        Phi0    = sqrt(1.-AA)
        lambda  = 1.

        ! Initialise scalar field values
        IF (PRESENT(scal0)) THEN
            phi%scal = scal0
        ELSE
            phi%scal = 0.
        ENDIF
        IF (PRESENT(dscal0)) THEN
            phi%dscal = dscal0
        ELSE
            phi%dscal = 0.
        ENDIF
        IF (PRESENT(d2scal0)) THEN
            phi%d2scal = d2scal0
        ELSE
            phi%d2scal = 0.
        ENDIF
    END SUBROUTINE Initialise
    
    REAL FUNCTION ScalarPotential(field)
        ! This is a container for the function of the potential of the scalar field.
        REAL, INTENT(IN) :: field
        
        !ScalarPotential = ((field-Phi0)*(field-Phi0) + AA) * exp(-lambda*field)
        ScalarPotential = 0.0
        RETURN
    END FUNCTION ScalarPotential

    REAL FUNCTION dScalarPotential(field)
        ! This function is the derivative of the potential with respect to the scalar field.
        REAL, INTENT(IN) :: field

        !dScalarPotential = field - Phi0
        !dScalarPotential = 2*dScalarPotential - lambda*(dScalarPotential*dScalarPotential + AA)
        !dScalarPotential = dScalarPotential * exp(-lambda*field)
        dScalarPotential = 0.0
        RETURN
    END FUNCTION dScalarPotential

END MODULE scalarfield_module
