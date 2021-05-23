MODULE variables

    use physics_parameters

    INTEGER, PARAMETER :: itmax = 200    !Maximum # of lines to be printed.

    !Time variables
    REAL, SAVE :: t, &          !Time.
            dt, &               !Time step.
            dlt9dt              !(1/t9)*d(t9)/d(t).

    !Dynamic variables
    REAL, SAVE, DIMENSION(14) :: thm    !Thermodynamic variables (energy densities).
    REAL, SAVE :: hubcst        !Expansion rate of the universe.

    !Energy densities.
    REAL, SAVE :: rhone0, &     !Initial electron neutrino energy density.
            rhob0, &            !Initial baryon energy density.
            rhob, &             !Baryon energy density.
            rnb                 !Baryon energy density (ratio to init value).

    !Matrix coefficients for linear equation.
    DOUBLE PRECISION, SAVE, DIMENSION(nnuc, nnuc) :: a
                                !Relates y(t+dt) to y(t).
    REAL, SAVE, DIMENSION(nnuc) :: &
            b, &        !Contains y0 in inverse order.
            yx          !yy in reverse order.

    !Flags, counters
    INTEGER, SAVE:: &
            ltime, &            !Indicates if output buffer printed.
            is, &               !# total iterations for particular model.
            ip, &               !# iterations after outputing a line.
            it, &               !# times accumulated in output buffer.
            mbad                !Indicates if gaussian elimination failed.

    !Computation location for check function
    INTEGER, SAVE :: itime      !Computation location.

    !Output arrays
    REAL, SAVE, DIMENSION(itmax,nnuc) :: xout   !Nuclide mass fractions.
    REAL, SAVE, DIMENSION(itmax, 6) :: thmout   !Thermodynamic variables.
    REAL, SAVE, DIMENSION(itmax) :: &
            t9out, &        !Temperature (in units of 10**9 K).
            tout, &         !Time.
            dtout, &        !Time step.
            etaout, &       !Baryon to photon ratio.
            hubout          !Expansion rate.

    !Output file status.
    INTEGER, SAVE :: nout       !Number of output requests.
    LOGICAL, SAVE :: outfile    !Indicates if output file used.

END MODULE variables