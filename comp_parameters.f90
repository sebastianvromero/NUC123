MODULE computational_parameters

    !Computation parameters and default values
    REAL, SAVE :: cy, &     !Time step limiting constant on abundances.
                  ct, &     !Time step limiting constant on temperature.
                  t9i, &    !Initial temperature (in 10**9 K).
                  t9f, &    !Final temperature (in 10**9 K).
                  ytmin     !Smallest abundances allowed.
    INTEGER, SAVE :: inc    !Accumulation increment.

    REAL, SAVE :: dt1       !Initial time step.

    !Run options.
    INTEGER, SAVE :: irun, &    !Run network size.
            isize, &            !Number of nuclides in computation.
            jsize               !Number of reactions in computation.

CONTAINS

    SUBROUTINE ResetComputationalParameters(irun)
    !Set computational parameters to default values

        use physics_parameters
        integer, intent(in) :: irun

        cy    = 0.300           !Time step limiting constant on abundances.
        ct    = 0.030           !Time step limiting constant on temperature.
        t9i   = 1.00e+02        !Initial temperature.
        t9f   = 1.00e-02        !Final temperature.
        ytmin = 1.00e-25        !Smallest abundances allowed.
        inc   = 30              !Accumulation increment.
        dt1   = 1.00e-04        !Initial time step.

        IF (irun.eq.1) THEN         !Maximal network size.
            isize = nnuc
            jsize = nrec
        ELSE IF (irun.eq.2) THEN    !Abridged network size.
            isize = 18
            jsize = 64
        ELSE IF (irun.eq.3) THEN    !Minimal network size.
            isize = 9
            jsize = 34
        END IF

    END SUBROUTINE ResetComputationalParameters

END MODULE computational_parameters
