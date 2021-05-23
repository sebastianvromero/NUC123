C===============IDENTIFICATION DIVISION====================

      SUBROUTINE DriverRK2
      
C-------LINKAGES.
C     CALLED BY - [subroutine] run
C     CALLS     - [subroutine] start, derivs, accum

C-------REMARKS.
C     Second order Runge-Kutta computational routine

      use computational_parameters
      use universe_module
      use variables

C-------PARAMETERS.
      PARAMETER (cl=1.e-16)        !Lower limit on size of time step.

C-------TIME AND TIME STEP VARIABLES.
      REAL    dtmin                !Mininum time step.
      REAL    dtl                 !Time step from limitation on abund changes.

C==================PROCEDURE DIVISION======================

!C10-----INPUT INITIALIZATION INFORMATION, RELABEL-------------------

      ltime = 0                    !Set termination indicator to zero.
      CALL start                   !Input initialization information.

!C20-----LOOP ONE----------------------------------------

 200  continue                     !Begin Runge-Kutta looping.
      loop = 1                     !Loop indicator.
C..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
      CALL derivs(loop)
      itime = 4                    !Time = 1st R-K loop.
      CALL check                   !Check interface subroutine.
C..........ACCUMULATE.
      IF ((uni%t9.le.t9f).or.                         !Low temp.
     |    (dt.lt.abs(cl/dlt9dt)).or.              !Small dt.
     |    (ip.eq.inc)) CALL accum                 !Enough iterations.
C..........POSSIBLY TERMINATE COMPUTATION.
      IF (ltime.eq.1) THEN         !Return to run selection.
        RETURN
      END IF
C..........RESET COUNTERS.
      IF (ip.eq.inc) THEN          !Reset iteration counters.
        ip = 0
      END IF
      ip = ip + 1
      is = is + 1
C..........ADJUST TIME STEP.
      IF (is.gt.3) THEN            !Adjust time step after 3 iterations.
        dtmin = abs(1./dlt9dt)*ct  !Trial value for minimum time step (Ref 1).
        DO i = 1,isize             !Go through all abundance changes.
          IF ((uni%dydt(i).ne.0.).and.(uni%y(i).gt.ytmin)) THEN
            dtl = abs(uni%y(i)/uni%dydt(i))*cy
     |            *(1.+(alog10(uni%y(i))/alog10(ytmin))**2)  !(Ref 2).
            IF (dtl.lt.dtmin) dtmin = dtl         !Find smallest time step.
          END IF
        END DO
        IF (dtmin.gt.1.5*dt) dtmin = 1.5*dt       !Limit change in time step.
        dt = dtmin                 !Set new time step.
      END IF
      t = t + dt                   !Increment time.
C..........STORE AND INCREMENT VALUES (Ref 3).

      uni0 = uni      

      CALL IncrementUniverse(uni, uni0, uni0, dt, ytmin)

!C30-----LOOP TWO----------------------------------------

      loop = 2                     !Step up loop counter.
C..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
      CALL derivs(loop)
      itime = 7                    !Time = 2nd R-K loop.
      CALL check                   !Check interface subroutine.
C..........INCREMENT VALUES.
      
      CALL StepRK2(uni, uni0, dt, ytmin)

      GO TO 200

C-------REFERENCES--------------------------------------
C     1)  Constraint on dt from the requirement that 
C                (d(t9)/dt)*(dt/t9) < ct
C         Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 293, equation C6.
C     2)  Constraint on dt from 
C                dtl < y/(dy/dt)*cy*(1+(log(y)/log(ytmin))**2)
C         Wagoner, R.V. 1969, page 293, equation C7 but with log term squared.
C     3)  Wagoner, R.V. 1969, page 292, equations C1, C2.

      END


      SUBROUTINE DriverRK4
      
C-------LINKAGES.
C     CALLED BY - [subroutine] run
C     CALLS     - [subroutine] start, derivs, accum

C-------REMARKS.
C     Fourth order Runge-Kutta computational routine

      use computational_parameters
      use universe_module
      use variables

C-------PARAMETERS.
      PARAMETER (cl=1.e-16)        !Lower limit on size of time step.

C-------Temporary RK variables.
C     Let f(t, y) = dy/dt
C       k0 = f(t0, y0)
C       k1 = f(t0 + dt/2, y0 + k0*dt/2)
C       k2 = f(t0 + dt/2, y0 + k1*dt/2)
C       k3 = f(t0 + dt, y0 + k2*dt)
C     Then
C       y1 = y0 + (k0/6 + k1/3 + k2/3 + k3/6) * dt

      TYPE(Universe) :: uni1, uni2, uni3

C-------TIME AND TIME STEP VARIABLES.
      REAL    dtmin                !Mininum time step.
      REAL    dtl                 !Time step from limitation on abund changes.

C==================PROCEDURE DIVISION======================

!C10-----INPUT INITIALIZATION INFORMATION, RELABEL-------------------

      ltime = 0                    !Set termination indicator to zero.
      CALL start                   !Input initialization information.

!C20-----LOOP ONE----------------------------------------

 200  continue                     !Begin Runge-Kutta looping.
      loop = 1                     !Loop indicator.
C..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
      CALL derivs(loop)
      itime = 4                    !Time = 1st R-K loop.
      CALL check                   !Check interface subroutine.
C..........ACCUMULATE.
      IF ((uni%t9.le.t9f).or.                         !Low temp.
     |    (dt.lt.abs(cl/dlt9dt)).or.              !Small dt.
     |    (ip.eq.inc)) CALL accum                 !Enough iterations.
C..........POSSIBLY TERMINATE COMPUTATION.
      IF (ltime.eq.1) THEN         !Return to run selection.
        RETURN
      END IF
C..........RESET COUNTERS.
      IF (ip.eq.inc) THEN          !Reset iteration counters.
        ip = 0
      END IF
      ip = ip + 1
      is = is + 1
C..........ADJUST TIME STEP.
      IF (is.gt.3) THEN            !Adjust time step after 3 iterations.
        dtmin = abs(1./dlt9dt)*ct  !Trial value for minimum time step (Ref 1).
        DO i = 1,isize             !Go through all abundance changes.
          IF ((uni%dydt(i).ne.0.).and.(uni%y(i).gt.ytmin)) THEN
            dtl = abs(uni%y(i)/uni%dydt(i))*cy
     |            *(1.+(alog10(uni%y(i))/alog10(ytmin))**2)  !(Ref 2).
            IF (dtl.lt.dtmin) dtmin = dtl         !Find smallest time step.
          END IF
        END DO
        IF (dtmin.gt.1.5*dt) dtmin = 1.5*dt       !Limit change in time step.
        dt = dtmin                 !Set new time step.
      END IF
C..........STORE AND INCREMENT VALUES (Ref 3).

      uni0 = uni

      ! Increment time to halfway point
      dt = 0.5 * dt
      t = t + dt
      CALL IncrementUniverse(uni, uni0, uni0, dt, ytmin)

!C30-----LOOP TWO----------------------------------------

      loop = 2                     !Step up loop counter.
C..........COMPUTE DERIVATIVES OF VARIABLES TO BE EVOLVED.
      CALL derivs(loop)
      itime = 7                    !Time = 2nd R-K loop.
      CALL check                   !Check interface subroutine.
C..........INCREMENT VALUES.

      uni1 = uni
      CALL IncrementUniverse(uni, uni0, uni1, dt, ytmin)

      CALL derivs(loop)
      uni2 = uni
      
      ! Increment time to end point
      t = t + dt
      dt = 2. * dt
      CALL IncrementUniverse(uni, uni0, uni2, dt, ytmin)
      
      CALL derivs(loop)
      uni3 = uni
      
      CALL StepRK4(uni, uni0, uni1, uni2, uni3, dt, ytmin)

      GO TO 200

C-------REFERENCES--------------------------------------
C     1)  Constraint on dt from the requirement that 
C                (d(t9)/dt)*(dt/t9) < ct
C         Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 293, equation C6.
C     2)  Constraint on dt from 
C                dtl < y/(dy/dt)*cy*(1+(log(y)/log(ytmin))**2)
C         Wagoner, R.V. 1969, page 293, equation C7 but with log term squared.
C     3)  Wagoner, R.V. 1969, page 292, equations C1, C2.

      END


C===============IDENTIFICATION DIVISION====================

      SUBROUTINE start

C-------LINKAGES.
C     CALLED BY - [subroutine] driver
C     CALLS     - [subroutine] rate1, bessel, rate0
C               - [function] ex

C-------REMARKS.
C     Sets initial conditions.

      use scalarfield_module
      use bessel_module
      use computational_parameters
      use reactions_module
      use universe_module
      use variables

C-------PARAMETERS.
      PARAMETER (const1=0.09615)   !Relation between time and temperature.
      PARAMETER (const2=6.673e-8) !Gravitational constant.

C-------LOCAL VARIABLES.
      REAL    z                  !Defined by z = m(electron)*c**2/k*t9.


C==================PROCEDURE DIVISION======================

!C10-----INITIALIZE FLAGS AND COUNTERS-------------------------

      ltime = 0                    !No output yet.
      is    = 1                    !First iteration coming up.
      ip    = inc                  !Set to maximum allowed # of iteration.
      it    = 0                    !No accumulation yet.
      mbad  = 0                    !No computational errors.

!C20-----SETTINGS----------------------------------------

C..........COMPUTATIONAL SETTINGS.
      uni%t9  = t9i                    !Initial temperature.
      tnu = uni%t9                     !Initial neutrino temperature.
      t   = 1/(const1*uni%t9)**2       !Initial time (Ref 1).
      dt  = dt1                    !Initial time step.
C..........MODEL SETTINGS.
      g   = const2*c(1)            !Modify gravitational constant.
      tau = c(2)                 !Convert n half-life (min) to lifetime (sec).
      tau = tau/0.98               !Coulomb correction (Ref 2).
      xnu = c(3)                   !Number of neutrino species.

      CALL Initialise(uni%scalarfield)

!C30-----COMPUTE INITIAL ABUNDANCES FOR NEUTRON AND PROTON--------------

      IF ((15.011/uni%t9+xi(1)).gt.58.) THEN      !Overabundance of antineutrinos.
        uni%y(1) = 1.e-25              !Very little of neutrons.
        uni%y(2) = 1.                  !Essentially all protons.
      ELSE
        IF ((15.011/uni%t9+xi(1)).lt.-58.) THEN   !Overabundance of neutrinos.
          uni%y(1) = 1.                !Essentially all neutrons.
          uni%y(2) = 1.e-25            !Very little of protons.
        ELSE
          uni%y(1) = 1./(ex(15.011/uni%t9+xi(1))+1.)  !Initial n abundance (Ref 3).
          uni%y(2) = 1./(ex(-15.011/uni%t9-xi(1))+1.) !Initial p abundance (Ref 3).
        END IF     
      END IF
      IF (xi(1).ne.0.) THEN        !Electron neutrino degeneracy.
        cnorm = 1.
        tnu   = .00001             !Low temperature.
        CALL rate1(0.00001)        !Find normalization constant at low temp.
        cnorm = 1/tau/f(1)
      END IF
      uni0%y(1) = uni%y(1)
      uni0%y(2) = uni%y(2)

!C40-----FIND RATIO OF BARYON DENSITY TO TEMPERATURE CUBED--------------

      z      = 5.930/uni%t9            !Inverse of temperature.
      CALL bessel(z)
      uni%hv     = 3.3683e+4*eta*2.75 !(Ref 4 but with final eta).
      uni%phie   = uni%hv*(1.784e-5*uni%y(2))  !Chemical potential of electron (Ref 5).
     |            /(.5*z**3*(bl(1)-2.*bl(2)+3.*bl(3)-4.*bl(4)+5.*bl(5)))
      rhob0  = uni%hv*uni%t9**3            !Baryon density.
      IF ((xi(1).eq.0.).and.(xi(2).eq.0.).and.(xi(3).eq.0)) THEN  !Nondegen.
        rhone0 = 7.366*uni%t9**4       !Electron neutrino density (Ref 6).
      END IF
        
!C50-----SET ABUNDANCES FOR REST OF NUCLIDES----------------------

      uni%y(3)  = uni%y(1)*uni%y(2)*rhob0*ex(q9(12)/uni%t9)/(.471e+10*uni%t9**1.5)  !(Ref 7).
      uni0%y(3) = uni%y(3)
      DO i = 4,isize
        uni%y(i)  = ytmin              !Set rest to minimum abundance.
        uni0%y(i) = uni%y(i)               !Init abundances at beginning of iteration.
      END DO
      CALL rate0                   !Compute weak decay rates.
      RETURN

C-------REFERENCES--------------------------------------
C     1) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 44, equation A15.
C     2) Coulomb correction obtained by dividing by correction factor Fp(t9)
C               Fp(t9) = 1 - 0.5(pi/(137<v>/c)) 
C          Wagoner, R.V. 1973, Ap. J. 179, page 358.
C     3) For the nondegenerate case:
C          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 4, equation 3.
C        For the case with neutrino degeneracy:
C          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 
C          page 417, equation 9.
C     4) Wagoner, R.V. 1969, Ap J. Suppl. No. 162, 18, page 250, equation 4.
C          3.3683e+4 = Mu(ng/t9**3) with Mu the atomic mass, ng the 
C          photon density.  2.75 is for the 11/4 factor difference
C          between the initial and final values of eta.
C     5) Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
C          Kellogg Radiation Lab preprint OAP-714,
C          equation D.2.
C     6) Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 43, equation A4.
C          7.366 is used instead of 14.73 as the latter is the sum total 
C          for 2 neutrino species.
C     7) Initial deuterium abundance from nuclear statistical equilibrium
C          Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 19, equation 17.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE derivs(loop)

C-------LINKAGES.
C     CALLED BY - [subroutine] driver
C     CALLS     - [subroutine] therm, rate1, rate4, rate3, rate2, sol

C-------REMARKS.
C     Computes derivatives of
C       - Temperature
C       - hv
C       - Chemical potential
C       - abundances

      use scalarfield_module
      use computational_parameters
      use universe_module
      use variables
      use reactions_module

C-------SUMS.
      REAL    sumy                 !Sum of abundances.
      REAL    sumzy                !Sum of charge*abundances.
      REAL    sumdy                !Sum of abundance flows.
      REAL    summdy               !Sum of (mass excess)*(abundance flows).
      REAL    sumzdy               !Sum of (charge)*(abundance flows).

C-------DERIVATIVES.
      REAL    dphdt9               !d(phi e)/d(t9).
      REAL    dphdln               !d(phi e)/d(h).
      REAL    dphdzy               !d(phi e)/d(sumzy).
      REAL    dlndt9               !(1/h)*d(h)/d(t9).
      REAL    bar                  !Baryon density and pressure terms.

C-------LOCAL VARIABLES.
      INTEGER loop                 !Counts which Runge-Kutta loop.

!      REAL    rhoscal              !Density of scalar field.
!      REAL    drhoscaldt9          !d(rho scal)/d(t9).
!      REAL    dpscaldt9            !d(p scal)/d(t9).


C==================PROCEDURE DIVISION======================

!C10-----COMPUTE DERIVATIVES FOR ABUNDANCES-----------------------

      rnb    = uni%hv*uni%t9*uni%t9*uni%t9/rhob0   !Baryon mass density (ratio to init value).
C..........VARIOUS THERMODYNAMIC QUANTITIES.
      CALL therm
      hubcst = sqrt((8./3.)*pi*g*(thm(10))+(cosmo/3.))  !Expansion rate.
      rhob   = thm(9)             !Baryon mass density.
C..........COMPUTE REACTION RATE COEFFICIENTS.
      CALL rate1(uni%t9)
      GO TO (100,110,120), irun    !Run network selection.
 100  CONTINUE
        CALL Rate4(uni%t9)         !Forward rate for all of reactions.
 110  CONTINUE
        CALL Rate3(uni%t9)         !Forward rate for reactions with A < 19.
 120  CONTINUE
        CALL Rate2(uni%t9)         !Forward rate for reactions with A < 10.

      CALL CheckMaxTemps(uni%t9)
C..........SOLVE COUPLED DIFFERENTIAL EQUATIONS.

      CALL sol(loop)
      IF (mbad.gt.0) RETURN        !Abort in case matrix not invertible.

!C20-----COMPUTE DERIVATIVES FOR TEMPERATURE, hv, AND CHEMICAL POTENTIAL------

C..........SCALAR FIELD VARIABLES.
      uni%scalarfield%d2scal = -3*hubcst*uni%scalarfield%dscal - dScalarPotential(uni%scalarfield%scal)

!      rhoscal      = 0
!      drhoscaldt9  = 0
!      dpscaldt9    = 0

C..........INITIALIZE SUMS TO ZERO.
      sumy   = 0.
      sumzy  = 0.
      sumdy  = 0.
      summdy = 0.
      sumzdy = 0.
C..........ACCUMULATE TO GET SUM.
      DO i = 1,isize
        sumy   = sumy   + uni%y(i)           !Sum of abundance.
        sumzy  = sumzy  + zm(i)*uni%y(i)     !Sum of charge*abundance.
        sumdy  = sumdy  + uni%dydt(i)        !Sum of abundance flow.
        summdy = summdy + dm(i)*uni%dydt(i)  !Sum of (mass excess)*(abundance flow).
        sumzdy = sumzdy + zm(i)*uni%dydt(i)  !Sum of (charge)*(abundance flow).
      END DO
C..........CHANGES IN TEMPERATURE, hv, AND CHEMICAL POTENTIAL.
      dphdt9 = thm(12)*(-1.070e-4*uni%hv*sumzy/uni%t9 - thm(11))
      dphdln = -thm(12)*3.568e-5*uni%hv*sumzy
      dphdzy = thm(12)*3.568e-5*uni%hv
      bar    = 9.25e-5*uni%t9*sumy + 1.388e-4*uni%t9*sumdy/(3.*hubcst)
     |         + summdy/(3.*hubcst)
      dlndt9 = -(thm(2) + thm(5) + thm(6)*dphdt9 + thm(9)*1.388e-4*
     |         sumy)/(thm(1) + thm(3) + thm(4) + thm(7) + thm(9)*bar
     |         + thm(6)*(dphdln + dphdzy*sumzdy/(3.*hubcst)))   !(Ref 1).
      uni%dt9    = (3.*hubcst)/dlndt9
      dlt9dt = uni%dt9/uni%t9
      uni%dhv    = -uni%hv*((3.*hubcst) + 3.*dlt9dt)                    !(Ref 2).
      uni%dphie  = dphdt9*uni%dt9 + dphdln*(3.*hubcst) + dphdzy*sumzdy  !(Ref 3).

      RETURN

C-------REFERENCES--------------------------------------
C     1)  Kawano, L., 1992, Fermilab preprint FERMILAB-PUB-92/04-A,
C          Kellogg Radiation Lab preprint OAP-714,
C          equation D.35.
C     2)  Kawano, L. 1992, preprint, equation D.19.
C     3)  Kawano, L. 1992, preprint, equation D.20.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE accum

C-------LINKAGES.
C     CALLED BY - [subroutine] driver
C     CALLS     - none

C-------REMARKS.
C     Output accumulator.

      use computational_parameters
      use universe_module
      use variables

C==================PROCEDURE DIVISION======================

      it = it + 1                  !Set up accumulation counter.

!C10-----SET UP OUTPUT VARIABLES-------------------------------

C..........DIVIDE NUMBER FRACTION BY THAT OF PROTON.
      DO i = 1,isize
        xout(it,i) = uni%y(i)/uni%y(2)
      END DO
      xout(it,2) = uni%y(2)*am(2)      !Exception for proton.
      xout(it,6) = uni%y(6)*am(6)      !Exception for helium.
C..........SUM UP ABUNDANCES OF HEAVY NUCLIDES.
      xout(it,10) =  xout(it,10)+xout(it,11)+xout(it,12)+xout(it,13)
     |              +xout(it,14)+xout(it,15)+xout(it,16)+xout(it,17)
     |              +xout(it,18)+xout(it,19)+xout(it,20)+xout(it,21)
     |              +xout(it,22)+xout(it,23)+xout(it,24)+xout(it,25)
     |              +xout(it,26)   !Li8 to O16.
C..........RELABEL TEMPERATURE, TIME, THERMODYNAMIC VARIABLES, ETC.
      t9out(it)    = uni%t9            !Temperature.
      tout(it)     = t             !Time.
      thmout(it,1) = thm(1)        !rho photon.
      thmout(it,2) = thm(4)        !rho electron.
      thmout(it,3) = thm(8)        !rho neutrino.
      thmout(it,4) = thm(9)        !rho baryon.
      thmout(it,5) = uni%phie          !Chemical potential.
      thmout(it,6) = thm(10)       !rho total.
      dtout(it)    = dt            !Time step.
      etaout(it)   = uni%hv/(3.3683e+4)!Baryon to photon ratio.
      hubout(it)   = hubcst        !Expansion rate.

      itime = 5
      CALL check

!C20-----INDICATE TERMINATION OF ACCUMULATION IF APPROPRIATE------------

      IF ((it.eq.itmax).or.(ip.lt.inc)) ltime = 1
      RETURN        

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE therm

C-------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [subroutine] bessel, nudens
C               - [function] ex

C-------REMARKS.         
C     Computes various temperature dependent thermodynamic quantities.

      use scalarfield_module
      use bessel_module
      use computational_parameters
      use universe_module
      use variables

C-------PARAMETER.
      PARAMETER (q=2.531)          !(mass(neutron)-mass(proton))/m(electron)

C-------LOCAL VARIABLE.
      REAL    z                    !Defined by z = m(electron)*c**2/k*t9.
      REAL    t9
      REAL    phie
      
      t9 = uni%t9
      phie = uni%phie
      
C==================PROCEDURE DIVISION======================

!C10-----COMPUTE FACTORS------------------------------------

      z = 5.930/t9                 !z = m(electron)c**2/k(t9).
      tnu = ((rnb)**(1./3.))*t9i   !Neutrino temperature.
C..........FACTORS OF z.
      z1 = z
      z2 = z*z
      z3 = z*z*z
      z4 = z*z*z*z
      z5 = z*z*z*z*z
C..........TRIGNOMETRIC FUNCTION VALUES.
      IF (phie.le.17.) THEN        !No chance of overflow.
        cosh1 = cosh(phie)
        cosh2 = cosh(2.*phie)
        cosh3 = cosh(3.*phie)
        cosh4 = cosh(4.*phie)
        cosh5 = cosh(5.*phie)   
        sinh1 = sinh(phie)
        sinh2 = sinh(2.*phie)
        sinh3 = sinh(3.*phie)
        sinh4 = sinh(4.*phie)
        sinh5 = sinh(5.*phie)   
      ELSE
        cosh1 = 0.
        cosh2 = 0.
        cosh3 = 0.
        cosh4 = 0.
        cosh5 = 0.
        sinh1 = 0.
        sinh2 = 0.
        sinh3 = 0.
        sinh4 = 0.
        sinh5 = 0.
      END IF
      CALL bessel(z)

!C20-----COMPUTE THERMODYNAMIC VARIABLES--------------------------

      thm(1)  = 8.418*t9*t9*t9*t9                               !(Ref 1).
      thm(2)  = 4.*thm(1)/t9                                    !(Ref 2).
      thm(3)  = thm(1)/3.                                       !(Ref 3).
      thm(4)  = 3206.*(bm(1)*cosh1 - bm(2)*cosh2 + bm(3)*cosh3        !(Ref 4).
     |          - bm(4)*cosh4 + bm(5)*cosh5)
      thm(5)  = 3206.*(z/t9)*(bn(1)*cosh1 - 2.*bn(2)*cosh2          !(Ref 5).
     |          + 3.*bn(3)*cosh3 - 4.*bn(4)*cosh4 + 5.*bn(5)*cosh5)
      thm(6)  = 3206.*(bm(1)*sinh1 - 2.*bm(2)*sinh2 + 3.*bm(3)*sinh3  !(Ref 6).
     |          - 4.*bm(4)*sinh4 + 5.*bm(5)*sinh5)                  
      thm(7)  = 3206.*(bl(1)*cosh1/z - bl(2)*cosh2/(2.*z)           !(Ref 7).
     |          + bl(3)*cosh3/(3.*z) - bl(4)*cosh4/(4.*z)
     |          + bl(5)*cosh5/(5.*z))                             
      IF ((xi(1).eq.0.).and.(xi(2).eq.0.).and.(xi(3).eq.0)) THEN  !Nondegen.
        thm(8) = xnu*rhone0*(rnb**(4./3.))                      !(Ref 8).
      ELSE                         !Include effects of neutrino degeneracy.
        thm(8) = 0.
        DO nu = 1,int(xnu)              !For every neutrino family.
          CALL nudens              !Compute neutrino energy density.
          thm(8) = thm(8) + 12.79264*rhonu  !Have 12.79264 from units change.
        END DO
      END IF
      thm(9)  = rhob0*rnb                                       !(Ref 9).
      thm(10) = thm(1) + thm(4) + thm(8) + thm(9)               !(Ref 10).
      thm(11) = -(z**3/t9)*(sinh1*(3.*bl(1)-z*bm(1))-sinh2*(3.*bl(2)  !(Ref 11).
     |          -2.*z*bm(2)) + sinh3*(3.*bl(3)-3.*z*bm(3)) - sinh4
     |          *(3.*bl(4)-4.*z*bm(4)) + sinh5*(3.*bl(5)-5.*z*bm(5)))
      thm(12) = z**3*(cosh1*bl(1)- 2.*cosh2*bl(2)                   !(Ref 12).
     |          + 3.*cosh3*bl(3) - 4.*cosh4*bl(4) + 5.*cosh5*bl(5))
      IF (thm(12).ne.0.) thm(12) = 1./thm(12)
      thm(13) = 1.000 + 0.565/z1 - 6.382/z2 + 11.108/z3         !(Ref 13).
     |          + 36.492/z4 + 27.512/z5
      thm(14) = (5.252/z1 - 16.229/z2 + 18.059/z3 + 34.181/z4   !(Ref 14).
     |          + 27.617/z5)*ex(-q*z)

C     Add scalar field density to total density
      thm(10) = thm(10) + 0.5*uni%scalarfield%dscal*uni%scalarfield%dscal
     |          + ScalarPotential(uni%scalarfield%scal)

      RETURN       

C-------REFERENCES AND NOTES-------------------------------
C     1)  thm(1)  = rho photon
C         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967, Ap. J. 148,
C          page 43, equation A2.)
C     2)  thm(2)  = d(rho photon)/d(t9)
C     3)  thm(3)  = (p photon)/c**2
C         (Wagoner, R.V., Fowler, W.A., and Hoyle, F. 1967,
C          page 43, equation A3.)
C     4)  thm(4)  = rho electron+positron
C         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9, 
C          page 281, equation B44.)
C     5)  thm(5)  = d(rho electron+positron)/d(t9)
C     6)  thm(6)  = d(rho electron+positron)/d(phi e)
C     7)  thm(7)  = (p electron+positron)/c**2
C         (Fowler, W.A. and Hoyle, F., 1964, Ap. J. Suppl. No. 91, 9, 
C          page 279, equation B27.)
C     8)  thm(8)  = rho neutrino
C                 = # neutrino species x rho electron neutrino (nondegenerate)
C                 = rho nu(e) + rho nu(m) + rho nu(t)          (degenerate)
C     9)  thm(9)  = rho baryon
C     10) thm(10) = rho total 
C                 = rho photon + rho electron+positron + rho neutrino 
C                              + rho baryon
C     11) thm(11) = d     /pi**2(hbar*c)**3(ne- - ne+)*z**3\
C                   d(t9) \  2  (mc**2)**3                 /
C     12) thm(12) = d        /pi**2(hbar*c)**3(ne- - ne+)*z**3\
C                   d(phi e) \  2  (mc**2)**3                 /
C     13) thm(13) = rate for n->p
C     14) thm(14) = rate for p->n

      END


C===============IDENTIFICATION DIVISION====================

      SUBROUTINE nudens

C-------LINKAGES.
C     CALLED BY - [subroutine] therm
C     CALLS     - [function] xintd, eval

C-------REMARKS.
C     Computes energy density contribution from neutrinos.

      use physics_parameters
      use computational_parameters

C-------PARAMTER.
      PARAMETER (iter=50)          !Number of gaussian quads.

C-------EXTERNAL FUNCTIONS.
      EXTERNAL func5               !Integral for neutrinos.
      EXTERNAL func6               !Integral for antineutrinos.


C=================DECLARATION DIVISION=====================

C-------LOCAL VARIABLES.
      REAL    uplim1               !Upper limit for neutrino energy integral.
      REAL    uplim2               !Upper limit for antineu energy integral.


C==================PROCEDURE DIVISION======================

!C10-----COMPUTE NEUTRINO ENERGY DENSITIES------------------------

      IF (abs(xi(nu)).le.0.03) THEN
C..........SMALL xi APPROXIMATION.
        rhonu = 2.*(3.14159**2/30.)*(tnu)**4
     |          *(7./8.+(15./(4*3.14159**2))*xi(nu)**2
     |          +(15./(8.*3.14159**4))*xi(nu)**4)
      ELSE
        IF (abs(xi(nu)).ge.30.) THEN
C..........LARGE xi APPROXIMATION.
          rhonu = ((tnu)**4)/(8.*3.14159**2)*xi(nu)**4
     |            *(1+12.*1.645 /xi(nu)**2)
        ELSE
C..........DO INTEGRATION
          uplim1 = (88.029+xi(nu))*tnu
          uplim2 = (88.029-xi(nu))*tnu
          IF (uplim2.le.0.) THEN
            rhonu = xintd(0.,uplim1,func5,iter)
          ELSE
            rhonu= xintd(0.,uplim1,func5,iter)
     |             + xintd(0.,uplim2,func6,iter)
          END IF
        END IF !(abs(xi(nu)).ge.30.) 
      END IF !(abs(xi(nu)).le.0.03) 
      RETURN

C-------REFERENCES--------------------------------------
C     Forms of the integrals involved can be found in
C       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
C       Freese, K., Kolb, E.W., Turner, M.S., 1983, Phys. Rev. D, 27, 1689.

      END



C===============IDENTIFICATION DIVISION====================

	FUNCTION eval()

C-------LINKAGES.
C     CALLED BY - [subroutine] rate1, nudens
C     CALLS     - [function] ex

C-------REMARKS.
C     Contains integrands to be integrated.

      use computational_parameters
      use variables

C=================DECLARATION DIVISION=====================

C-------FUNCTIONS TO BE EVALUATED.
      REAL   func1                 !1st part n->p rate.
      REAL   func2                 !2nd part n->p rate.
      REAL   func3                 !1st part p->n rate.
      REAL   func4                 !2nd part p->n rate.
      REAL   func5                 !Energy density for neutrinos.
      REAL   func6                 !Energy density for antineutrinos.

C-------LOCAL VARIABLES.
      REAL    x                    !Value at which function is evaluated.
      REAL    part1                !Exponential expression with photon temp.
      REAL    part2                !Exponential expression with neutrino temp.


C==================PROCEDURE DIVISION======================

!C10-----1ST PART OF INTEGRAL FOR n->p RATE-----------------------

      ENTRY func1(x)
      IF (x.le.0.) THEN
        func1 = 0.
      ELSE
        part1 = 1./(1.+ex(-.511*x/t9mev))
        part2 = 1./(1.+ex(+(x-2.531)*(.511/tnmev)-xi(1)))
        func1 = cnorm*x*(x-2.531)**2*(x**2-1)**.5*part1*part2
      END IF
      RETURN

!C20-----2ND PART OF INTEGRAL FOR n->p RATE-----------------------

      ENTRY func2(x)
      IF (x.le.1.) THEN
        func2 = 0.
      ELSE
        part1 = 1./(1.+ex(+.511*x/t9mev))
        part2 = 1./(1.+ex(-(x+2.531)*(.511/tnmev)-xi(1)))
        func2 = cnorm*x*(x+2.531)**2*(x**2-1)**.5*part1*part2
      END IF
      RETURN

!C30-----1ST PART OF INTEGRAL FOR p->n RATE-----------------------

      ENTRY func3(x)
      IF (x.le.1.) THEN
        func3 = 0.
      ELSE
        part1 = 1./(1.+ex(-.511*x/t9mev))
        part2 = 1./(1.+ex(+(x+2.531)*(.511/tnmev)+xi(1)))
        func3 = cnorm*x*(x+2.531)**2*(x**2-1)**.5*part1*part2
      END IF
      RETURN

!C40-----2ND PART OF INTEGRAL FOR p->n RATE-----------------------

      ENTRY func4(x)
      IF (x.le.1.) THEN
        func4 = 0.
      ELSE
        part1 = 1./(1.+ex(+.511*x/t9mev))
        part2 = 1./(1.+ex(-(x-2.531)*(.511/tnmev)+xi(1)))
        func4 = cnorm*x*(x-2.531)**2*(x**2-1)**.5*part1*part2
      END IF
      RETURN

!C50-----INTEGRAL FOR ENERGY DENSITY OF NEUTRINO---------------------

      ENTRY func5(x)
      func5 = 1./(2*3.14159**2)*x**3/(1.+exp(x/tnu-xi(nu)))
      RETURN

!C60-----INTEGRAL FOR ENERGY DENSITY OF ANTINEUTRINO-----------------

      ENTRY func6(x)
      func6 = 1./(2*3.14159**2)*x**3/(1.+exp(x/tnu+xi(nu)))
      RETURN

C-------REFERENCES--------------------------------------
C     Forms of the integrals involved can be found in
C       Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683.
C       Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.

      END



C===============IDENTIFICATION DIVISION====================

      FUNCTION xintd (xlow,xhi,func,nq)

C-------LINKAGES.
C     CALLED BY - [subroutine] rate1, nudens
C     CALLS     - none

C-------REMARKS.
C     Computes the integral of the function "func".


C=================DECLARATION DIVISION=====================

C-------INPUT VARIABLES.
	external func
      REAL    xlow                 !Array of low limits.
      REAL    xhi                  !Array of high limits.
      INTEGER nq                   !Number of six point gaussian quads.

C-------COMPUTATION VARIABLES.
      REAL    dist                 !Size of quad interval.
      REAL    cent                 !Center of quad interval.
      REAL    x                    !Variables of integration.
      REAL    sum                  !Summation of terms.

C-------COUNTERS.
      INTEGER nint                 !Interval number.
      INTEGER npnt                 !Point number.
      INTEGER np                   !Total number of points in interval.

C-------ABSCISSAS AND WEIGHT FACTORS.
      REAL    u(6)                 !Abscissas.
      REAL    w(6)                 !Weight factor.


C=====================DATA DIVISION========================

C-------ABSCISSAS AND WEIGHT FACTORS.
      DATA u/-.93246951420315,-.66120938646627,-.23861918608320,
     |        .23861918608320, .66120938646627, .93246951420315/  
      DATA w/.17132449237917,.36076157304814,.46791393457269,
     |       .46791393457269,.36076157304814,.17132449237917/        
      DATA np/6/              !6 point Gaussian integration.


C==================PROCEDURE DIVISION======================

!C10-----DO INTEGRATION-------------------------------------

      sum   = 0.       
      dist  = (xhi-xlow)/float(nq) !Size of quad interval.
      DO nint = 1,nq
        cent = xlow+(float(nint)-0.5)*dist  !Center of interval.
        DO npnt = 1,np
          x   = cent+0.5*dist*u(npnt) !Integration point.
          f   = func(x)            !Evaluate function x(1).
          sum = sum+f*w(npnt)      !Add up sum.
        END DO
      END DO

!C20-----GET INTEGRAL VALUE---------------------------------

      xintd = sum*dist*0.5         !Do integral.
      RETURN        

      END 



C===============IDENTIFICATION DIVISION====================

      FUNCTION ex(x)

C-------LINKAGES.
C     CALLED BY - [subroutine] start, rate2, rate3, rate4, sol
C               - [function] eval
C     CALLS     - none

C-------REMARKS.
C     Exponential function with underflow precaution.


C==================PROCEDURE DIVISION======================

      IF (x.gt.88.029) THEN        !In danger of overflow.
        ex = exp(88.029)
      ELSE
        IF (x.lt.-88.722) THEN     !In danger of underflow.
          ex = 0.
        ELSE                       !Value of x in allowable range.
          ex = exp(x)
        END IF
      END IF
      RETURN       

C-------NOTE-----------------------------------------
C     The overflow limit for the VAX/VMS system is exp(88.029).
C     The underflow limit for the VAX/VMS system is exp(-88.722).

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE sol(loop)

C-------LINKAGES.
C     CALLED BY - [subroutine] derivs
C     CALLS     - [subroutine] eqslin
C               - [function] ex

C-------REMARKS.
C     Computes reverse strong and electromagnetic reaction rates.
C     Fills and solves matrix equation for dydt(i).

      use reactions_module
      use universe_module
      use variables

      INTEGER isize1               !Equals isize + 1.

C-------EVOLUTION EQUATION COEFFICIENTS.
      INTEGER i,j,k,l              !Equate to ii,jj,kk,ll.
      REAL    ri,rj,rk,rl          !Equate to si,sj,sk,sl.
      REAL    ci,cj,ck,cl          !Coefficients of rate equation.

C-------LOCAL VARIABLES.
      REAL    yy(nnuc)             !Abundances at end of iteration.
      REAL    bdln                 !(10**(-5))*volume expansion rate.
      INTEGER ind                  !Equate to iform.
      INTEGER ierror               !Element which does not converge.

      REAL, DIMENSION(nnuc) :: y, y0
      REAL    t9

      t9 = uni%t9
      y  = uni%y
      y0 = uni0%y

C==================PROCEDURE DIVISION======================

!C10-----TEMPERATURE FACTORS AND INITIAL VALUES--------------

C..........TEMPERATURE FACTORS.
      t932  = t9**1.5              !t9**(3/2).
      t9m32 = 1./t932              !t9**(-3/2).
C..........MATRIX SIZE.
      isize1 = isize + 1
C..........INITIALIZE A-MATRIX.
      DO i = 1,isize
        DO j = 1,isize
          a(j,i) = 0.d0            !Set a-matrix to zero.
        END DO
      END DO

!C20-----COMPUTE FACTORS FOR THE A-MATRIX-------------------------

      DO n = 1,jsize
C..........EQUATE VARIABLES TO ARRAYS.
        ind = iform(n)             !Type of reaction.
        i = ii(n)                  !ID # of incoming nuclide i.
        j = jj(n)                  !ID # of incoming nuclide j.
        k = kk(n)                  !ID # of outgoing nuclide k.
        l = ll(n)                  !ID # of outgoing nuclide l.
        IF ((ind.ne.0).and.(i.le.isize).and.(l.le.isize)) THEN !Reaction okay.
          ri = si(ind)             !# of incoming nuclide i.
          rj = sj(ind)             !# of incoming nuclide j.
          rk = sk(ind)             !# of outgoing nuclide k.
          rl = sl(ind)             !# of outgoing nuclide l.
C..........COMPUTE DIFFERENT REACTION RATES.
          GO TO (201,202,203,204,205,206,207,208,209,210,211),ind
 201      CONTINUE                 !1-0-0-1 configuration.
            ci = f(n)              !(Ref 1).
            cj = 0.
            ck = 0.
            cl = r(n)
            GO TO 212
 202      CONTINUE                 !1-1-0-1 configuration.
            r(n) = rev(n)*1.e+10*t932*ex(-q9(n)/t9)*f(n)  !(Ref 2).
            f(n) = rhob*f(n)
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = 0.
            cl = r(n)
            GO TO 212
 203      CONTINUE                 !1-1-1-1 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*ex(-q9(n)/t9)*f(n)  !(Ref 3).
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = y(l)*r(n)/2.
            cl = y(k)*r(n)/2.
            GO TO 212
 204      CONTINUE                 !1-0-0-2 configuration.
            ci = f(n)
            cj = 0.
            ck = 0.
            cl = y(l)*r(n)/2.
            GO TO 212
 205      CONTINUE                 !1-1-0-2 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*ex(-q9(n)/t9)*f(n)  !(Ref 3).
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = 0.
            cl = y(l)*r(n)/2.
            GO TO 212
 206      CONTINUE                 !2-0-1-1 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*ex(-q9(n)/t9)*f(n)  !(Ref 3).
            ci = y(i)*f(n)/2.
            cj = 0.
            ck = y(l)*r(n)/2.
            cl = y(k)*r(n)/2.
            GO TO 212
 207      CONTINUE                 !3-0-0-1 configuration.
            r(n) = rev(n)*1.e+20*t932*t932*ex(-q9(n)/t9)*f(n)  !(Ref 4).
            f(n) = rhob*rhob*f(n)
            ci = y(i)*y(i)*f(n)/6.
            cj = 0.
            ck = 0.
            cl = r(n)
            GO TO 212
 208      CONTINUE                 !2-1-0-1 configuration.
            r(n) = rev(n)*1.e+20*t932*t932*ex(-q9(n)/t9)*f(n)  !(Ref 4).
            f(n) = rhob*rhob*f(n)
            ci = y(j)*y(i)*f(n)/3.
            cj = y(i)*y(i)*f(n)/6.
            ck = 0.
            cl = r(n)
            GO TO 212
 209      CONTINUE                 !1-1-1-2 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*1.e-10*t9m32*rhob*ex(-q9(n)/t9)*f(n)  !(Ref 5).
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = y(l)*y(l)*r(n)/6.
            cl = y(k)*y(l)*r(n)/3.
            GO TO 212
 210      CONTINUE                 !1-1-0-3 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*1.e-10*t9m32*rhob*ex(-q9(n)/t9)*f(n)  !(Ref 5).
            ci = y(j)*f(n)/2.
            cj = y(i)*f(n)/2.
            ck = 0.
            cl = y(l)*y(l)*r(n)/6.
            GO TO 212
 211      CONTINUE                 !2-0-2-1 configuration.
            f(n) = rhob*f(n)
            r(n) = rev(n)*1.e-10*t9m32*rhob*ex(-q9(n)/t9)*f(n)  !(Ref 5).
            ci = y(i)*f(n)/2.
            cj = 0.
            ck = y(l)*y(k)*r(n)/3.
            cl = y(k)*y(k)*r(n)/6.
 212      CONTINUE

!C30-----CONSTRUCT THE A-MATRIX--------------------------------

          i = isize1 - i           !Invert i index.
          j = isize1 - j           !Invert j index.
          k = isize1 - k           !Invert k index.
          l = isize1 - l           !Invert l index.
C..........FILL I NUCLIDE COLUMN.
          IF (j.le.isize) a(j,i) = a(j,i) +  rj*ci
          IF (k.le.isize) a(k,i) = a(k,i) -  rk*ci
          a(i,i) = a(i,i) +  ri*ci
          a(l,i) = a(l,i) -  rl*ci
C..........FILL J NUCLIDE COLUMN.
          IF (j.le.isize) THEN
            a(j,j) = a(j,j) +  rj*cj
            IF (k.le.isize) a(k,j) = a(k,j) -  rk*cj
            a(i,j) = a(i,j) +  ri*cj
            a(l,j) = a(l,j) -  rl*cj
          END IF
C..........FILL K NUCLIDE COLUMN.
          IF (k.le.isize) THEN
            IF (j.le.isize) a(j,k) = a(j,k) -  rj*ck
            a(k,k) = a(k,k) +  rk*ck
            a(i,k) = a(i,k) -  ri*ck
            a(l,k) = a(l,k) +  rl*ck
          END IF
C..........FILL L NUCLIDE COLUMN.
          IF (j.le.isize) a(j,l) = a(j,l) -  rj*cl
          IF (k.le.isize) a(k,l) = a(k,l) +  rk*cl
          a(i,l) = a(i,l) -  ri*cl
          a(l,l) = a(l,l) +  rl*cl
        END IF !((ind.ne.0).and.(i.le.isize).and.(l.le.isize)) 
      END DO !n = 1,jsize

!C40-----PUT A-MATRIX AND B-VECTOR IN FINAL FORM OF MATRIX EQUATION--------

      bdln   = 1.e-5*(3.*hubcst)   !(10**(-5))*(Expansion rate).
      DO i = 1,isize
        i1 = isize1 - i            !Invert the rows.
        DO j = 1,isize
          j1 = isize1 - j          !Invert the columns.
          IF (dabs(a(j,i)).lt.bdln*y0(j1)/y0(i1)) THEN
            a(j,i) = 0.d0          !Set 0 if tiny.
          ELSE
            a(j,i) = a(j,i)*dt     !Bring dt over to other side.
          END IF
        END DO
        a(i,i) = 1.d0 + a(i,i)     !Add identity matrix to a-matrix.
        b(i1)  = y0(i)             !Initial abundances.
      END DO

!C50-----SOLVE EQUATIONS TO GET DERIVATIVE------------------------

C..........SET MONITOR FLAG AND SOLVE BY GAUSSIAN ELIMINATION.
      IF (loop.eq.1) THEN
        CALL eqslin(ip,ierror)
      ELSE
        CALL eqslin(0,ierror)
      END IF

C..........OBTAIN DERIVATIVE.
      DO i = 1,isize
        yy(i)   = yx(isize1-i)     !Abundance at t+dt.
        uni%dydt(i) = (yy(i) - y0(i))/dt         !Take derivative.
      END DO

!C60-----POSSIBLE ERROR MESSAGES AND EXIT-------------------------

      IF (mbad.ne.0) THEN          !Problem in gaussian elimination.
        IF (mbad.eq.-1) print 6000, ierror !Error message.
        IF (mbad.ge. 1) print 6002, mbad   !Error message.
 6000   FORMAT (' ','** y(', i2, ') fails to converge **')
 6002   FORMAT (' ','** ', i2, ' th diagonal term equals zero **')
      END IF
      RETURN

C-------REFERENCES--------------------------------------
C     1) The coefficients are given in general as:
C             ci = ri*(y(j)**rj)*(y(i)**(ri-1)*f(n)/
C                  ((ri+rj)*fac(ri)*fac(rj))
C             cj = rj*(y(i)**ri)*(y(j)**(rj-1)*f(n)/
C                  ((ri+rj)*fac(ri)*fac(rj))
C             ck = rk*(y(l)**rl)*(y(k)**(rk-1)*f(n)/
C                  ((rk+rl)*fac(rk)*fac(rl))
C             cl = rl*(y(k)**rk)*(y(l)**(rl-1)*f(n)/
C                  ((rk+rl)*fac(rk)*fac(rl))
C        in which fac(x) is the factorial of x.
C     2) Form of reverse rate given in 
C        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
C          tables 1B, 4B, 7B.
C     3) Form of reverse rate given in
C        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
C          tables 2B, 3B, 5B, 6B, 8B, 9B, 10B.
C     4) Form of reverse rate given in
C        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
C          table 11B.
C     5) Form of reverse rate given in
C        Wagoner, R.V.,1969, Ap. J. Suppl. No. 162, 18, 247,
C          tables 12B, 13B.


      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE eqslin(icnvm,ierror)

C-------LINKAGES.
C     CALLED BY - [subroutine] sol
C     CALLS     - none

C-------REMARKS.
C     Solves for new abundances using gaussian elimination
C     with back substitution, no pivoting.

      use computational_parameters
      use variables

C-------PARAMETERS.
      PARAMETER (mord=1)           !Higher order in correction.
      PARAMETER (eps=2.e-4)        !Tolerance for convergence (.ge. 1.e-7).

C=================DECLARATION DIVISION=====================

C-------LOCAL MATRICES AND VECTORS.
      DOUBLE PRECISION a0(nnuc,nnuc)!Coefficient array w/o manipulation.
      DOUBLE PRECISION x(nnuc)     !Right-hand vector.

C-------LOCAL COMPUTATION VARIABLES.
      DOUBLE PRECISION cx          !Scaling factor in triangularization.
      DOUBLE PRECISION sum         !Sum for backsubstitution.
      REAL   xdy                   !Relative error.
 
C-------LOCAL COUNTERS.
      INTEGER nord                 !Order of correction.
      INTEGER icnvm                !Convergence monitor.
      INTEGER ierror               !ith nuclide fails to converge.


C==================PROCEDURE DIVISION======================

!C10-----INITIALIZE VECTOR----------------------------------

C..........SET COUNTERS TO ZERO.
      nord = 0                     !No corrections yet.
      mbad = 0                     !No errors yet.
C..........SET RIGHT-HAND AND SOLUTION VECTORS TO INITIAL VALUES.
      DO i = 1,isize
        x(i) = b(i)                !Right-hand vector.
        yx(i) = 0.                  !Solution vector.
      END DO
C..........SAVE MATRIX.
      IF (icnvm.eq.inc) THEN       !Monitor convergence.
        DO i = 1,isize
          DO j = 1,isize
            a0(j,i) = a(j,i)       !Initial value of coefficient array.
          END DO
        END DO
      END IF

!C20-----TRIANGULARIZE MATRIX AND SAVE OPERATOR----------------------

C..........CHECK TO SEE THAT THERE ARE NO ZEROES AT PIVOT POINTS.
      DO i = 1,isize-1
        IF (a(i,i).eq.0.d0) THEN   !Don't want to divide by zero.
          mbad = i                 !Position of zero coefficient.
          RETURN                   !Terminate matrix evaluation.
        END IF
C..........TRIANGULARIZE MATRIX.
        DO j = i+1,isize
          IF (a(j,i).ne.0.d0) THEN !Progress diagonally down the column.
            cx = a(j,i)/a(i,i)     !Scaling factor down the column.
            DO k = i+1,isize       !Progress diagonally along row.
              a(j,k) = a(j,k) - cx*a(i,k)  !Subtract scaled coeff along row.
            END DO
            a(j,i) = cx            !Scaled coefficient.
C..........OPERATE ON RIGHT-HAND VECTOR.
            x(j) = x(j) - cx*x(i)  !Subtract off scaled coefficient.
          END IF
        END DO
      END DO

!C30-----DO BACK SUBSTITUTION-------------------------------

 300  CONTINUE
      x(isize) = x(isize)/a(isize,isize)   !Solution for ultimate position.
      yx(isize) = yx(isize) + x(isize)
      DO i = isize-1,1,-1          !From i = penultimate to i = 1.
        sum = 0.d0
        DO j = i+1,isize
          sum = sum + a(i,j)*x(j)  !Sum up all previous terms.
        END DO
        x(i) = (x(i) - sum)/a(i,i) 
        yx(i) = yx(i) + x(i)         !Add difference to initial value.
      END DO

!C40-----TESTS AND EXITS------------------------------------

      IF (icnvm.eq.inc) THEN
        DO i = 1,isize
          IF (yx(i).ne.0.) THEN
            xdy = dabs(x(i)/yx(i))  !Relative error.
            IF (xdy.gt.eps) THEN
              IF (nord.lt.mord) THEN !Continue to higher orders.
                nord = nord + 1
C..........FIND ERROR IN RIGHT-HAND VECTOR.
                DO j = 1,isize
                  r = 0.d0         !Initialize r.
                  DO k = 1,isize
                    r = r + a0(j,k)*yx(k) !Left side with approximate solution.
                  END DO
                  x(j) = b(j) - r  !Subtract difference from right side.
                END DO
C..........OPERATE ON RIGHT-HAND VECTOR.
                DO j = 1,isize-1
                  DO k = j+1,isize
                   x(k) = x(k) - a(k,j)*x(j) !Subtract off scaled coefficient.
                  END DO
                END DO
                GO TO 300       !Go for another iteratiion.
              ELSE
C..........NOT ENOUGH CONVERGENCE.
                mbad = -1          !Signal error problem.
                ierror = i         !ith nuclide for which x/yx checked.
                RETURN
              END IF !(nord.lt.mord)
            END IF !(xdy.gt.eps)
          END IF !(yx(i).ne.0)
        END DO !i = 1,isize
      END IF !(icnvm.eq.inc)
      RETURN                       !No more iterations & relative error small.

      END

C===============IDENTIFICATION DIVISION====================

      SUBROUTINE rate0

C-------LINKAGES.
C     CALLED BY - [subroutine] start
C     CALLS     - none

C-------REMARKS.
C     Generates weak decay rates.

      use reactions_module, only : f
      
C==================PROCEDURE DIVISION======================

!C10-----SET DECAY RATE COEFFICIENTS---------------------------

C.......H3 -> e- + v + He3.........(Tilly-Weller-Hasan 1987)
      f(2)  = 1.79e-9

C.......Li8 -> e- + v + 2He4.......(Ajzenberg-Selove 1988)
      f(3)  = 8.27e-1

C.......B12 -> e- + B + C12........(Ajzenberg-Selove 1990)
      f(4)  = 3.43e+1

C.......C14 -> e- + v + N14........(Ajzenberg-Selove 1986)
      f(5)  = 3.834e-12

C.......B8 -> e+ + v + 2He4........(Ajzenberg-Selove 1988)
      f(6)  = 9.00e-1

C.......C11 -> e+ + v + B11........(Ajzenberg-Selove 1990)
      f(7)  = 5.668e-4

C.......N12 -> e+ + v + C12........(Ajzenberg-Selove 1990)
      f(8)  = 6.301e+1

C.......N13 -> e+ + v + C13........(Ajzenberg-Selove 1986)
      f(9)  = 1.159e-3

C.......O14 -> e+ + v + N14........(Ajzenberg-Selove 1986)
      f(10) = 9.8171e-3

C.......O15 -> e+ + v + N15........(Ajzenberg-Selove 1986)
      f(11) = 5.6704e-3

      RETURN

C-------REFERENCES--------------------------------------
C     Ajzenberg-Selove, F., 1990, Nucl. Phys. A506, 1.
C     Ajzenberg-Selove, F., 1988, Nucl. Phys. A490, 1.
C     Ajzenberg-Selove, F., 1986, Nucl. Phys. A449, 1.
C     Tilley, D.R., Weller, H.R., and Hasan, H.H., 1987, Nucl. Phys. A474, 1.

      END



C===============IDENTIFICATION DIVISION====================

      SUBROUTINE rate1(tph)

C-------LINKAGES.
C     CALLED BY - [subroutine] start, derivs
C     CALLS     - [function] xintd, eval

C-------REMARKS.
C     Generates rate coefficients for weak n->p and p->n reactions.

      use physics_parameters
      use computational_parameters
      use reactions_module, only : f, r
      use variables, only : thm

C-------EXTERNAL FUNCTIONS.
      EXTERNAL func1               !Part 1 of n->p rate.
      EXTERNAL func2               !Part 2 of n->p rate.
      EXTERNAL func3               !Part 1 of p->n rate.
      EXTERNAL func4               !Part 2 of p->n rate.


C=================DECLARATION DIVISION=====================

C-------LOCAL VARIABLES.
      REAL    tph                  !Photon temperature.
      REAL    w(2),x(2),          !Upper limits for exponentials, forward rate.
     |        y(2),z(2)           !Upper limits for exponentials, reverse rate.
      REAL    uplim1,uplim2,      !Upper limits for integrals for forward rate.
     |        uplim3,uplim4       !Upper limits for integrals for reverse rate.
      REAL    part1,part2,         !Parts of integrals for forward rate.
     |        part3,part4          !Parts of integrals for reverse rate.


C==================PROCEDURE DIVISION======================

C10-----COMPUTE WEAK REACTION RATES (NONDEGENERATE)-----------------

      IF (xi(1).eq.0.) THEN
        f(1)  = thm(13)/tau        !Forward rate for weak np reaction.
        r(1)  = thm(14)/tau        !Reverse rate for weak np reaction.
      ELSE

!C20-----COMPUTE WEAK REACTION RATES (DEGENERATE)--------------------

        t9mev = tph*.086171        !Convert photon temp to units of MeV.
        tnmev = tnu*.086171        !Convert neutrino temp to units of MeV.
C..........COMPUTE OVERFLOW LIMITS FOR LIMITS OF INTEGRATION (Ref 1 & 2).
        w(1) = (-(t9mev/.511)*(-88.722))
        w(2) = ((tnmev/.511)*(88.029+xi(1))+2.531)
        x(1) = ((t9mev/.511)*(88.029))
        x(2) = (-(tnmev/.511)*(-88.722+xi(1))-2.531)
        y(1) = (-(t9mev/.511)*(-88.722))
        y(2) = ((tnmev/.511)*(88.029-xi(1))-2.531)
        z(1) = ((t9mev/.511)*(88.029))
        z(2) = (-(tnmev/.511)*(-88.722-xi(1))+2.531)
C..........COMPARE LIMITS AND TAKE LARGER OF THE TWO.
        uplim1 = abs(w(1))
        uplim2 = abs(x(1))
        uplim3 = abs(y(1))
        uplim4 = abs(z(1))
        IF (uplim1.lt.abs(w(2))) uplim1 = w(2)
        IF (uplim2.lt.abs(x(2))) uplim2 = x(2)
        IF (uplim3.lt.abs(y(2))) uplim3 = y(2)
        IF (uplim4.lt.abs(z(2))) uplim4 = z(2)
C..........EVALUATE THE INTEGRALS NUMERICALLY.
        part1 = xintd(1.,uplim1,func1,iter)
        part2 = xintd(1.,uplim2,func2,iter)
        part3 = xintd(1.,uplim3,func3,iter)
        part4 = xintd(1.,uplim4,func4,iter)
        f(1) = part1 + part2       !Add 2 integrals to get forward rate.
        r(1) = part3 + part4       !Add 2 integrals to get reverse rate.
      END IF !(xi(1).eq.0.)
      RETURN

C-------REFERENCES--------------------------------------
C     1) Forms of the integrals involved can be found in
C          Scherrer,R.J., 1983, Mon.Not.R.astr.Soc., 205, 683.
C          Beaudet,G. and Goret,P., 1976, Astron. & Astrophys., 49, 415.
C
C     2) The overflow limit for the VAX/VMS system is exp(88.029).
C        The underflow limit for the VAX/VMS system is exp(-88.722).

      END


C===============IDENTIFICATION DIVISION====================

      SUBROUTINE check

C-------REMARKS.
C     This is an interface subroutine,
C     a flexible module which allows user to manipulate physical quantities
C     of interest at certain key points during the computer run.
C     Included within this subroutine is a roster of all global variables 
C     and their respective COMMON areas.
C     Current instructions accumulate abundances of deuterium, helium-3,
C     helium-4, and lithium-7 for eventual plotting, taking into account
C     the contribution of beryllium-7 to lithium-7 and tritium to helium-3.

      use scalarfield_module
      use computational_parameters
      use reactions_module
      use universe_module
      use variables

      integer, save:: stage

C==================PROCEDURE DIVISION======================

!C10-----OPEN FILE---------------------------------------

      IF (itime.eq.1) THEN         !Beginning of program.
        OPEN (unit=3, file='nucint.dat',  status='old')
      END IF

      IF (.false. .and. itime.eq.2) THEN
        OPEN (unit=4, file='output.txt',  status='unknown')

        OPEN (unit=10, file='d_creation.txt',  status='unknown')
        OPEN (unit=11, file='d_destruction.txt',  status='unknown')
        OPEN (unit=12, file='t_creation.txt',  status='unknown')
        OPEN (unit=13, file='t_destruction.txt',  status='unknown')
        OPEN (unit=14, file='He3_creation.txt',  status='unknown')
        OPEN (unit=15, file='He3_destruction.txt',  status='unknown')
        OPEN (unit=16, file='He4_creation.txt',  status='unknown')
        OPEN (unit=17, file='He4_destruction.txt',  status='unknown')
        OPEN (unit=18, file='Li6_creation.txt',  status='unknown')
        OPEN (unit=19, file='Li6_destruction.txt',  status='unknown')
        OPEN (unit=20, file='Li7_creation.txt',  status='unknown')
        OPEN (unit=21, file='Li7_destruction.txt',  status='unknown')
        OPEN (unit=22, file='Be7_creation.txt',  status='unknown')
        OPEN (unit=23, file='Be7_destruction.txt',  status='unknown')
      END IF
      
      IF (itime.eq.3) THEN
        stage = 4
      END IF

C Start of RK loop - make dynamic changes here
      IF (.false. .and. itime.eq.4) THEN
        if((uni%t9 .lt. 0.833) .and. (stage.eq.0)) then
           print *, "T=", uni%t9, "; varying lnvar=", lnvar
           CALL InitialiseReactions
           CALL SetDeutronBindingVariation(lnvar)
           stage = 1
        else if((uni%t9 .lt. 0.01) .and. (stage.eq.1)) then
           print *, "T=", uni%t9, "; stop varying"
           CALL InitialiseReactions
           stage = 2
        end if        
      END IF

!C20-----PRINTINTO FILE------------------------------------

      IF (.false. .and. itime.eq.5) THEN
!         call rate2(uni%t9)

c      d creation: SKM Figure 4a
         write(10,"(5(e13.5,' '))") uni%t9, f(12)*uni%y(1)*uni%y(2)/hubcst,
     |      r(20)*uni%y(5)/hubcst, r(13)*uni%y(4)/hubcst, r(25)*uni%y(7)/hubcst
c      d destruction: SKM Figure 4b
         write(11,"(9(e13.5,' '))") uni%t9, r(12)*uni%y(3)/hubcst, f(30)*uni%y(3)*uni%y(4)/hubcst,
     |      f(28)*uni%y(3)*uni%y(3)/hubcst, f(29)*uni%y(3)*uni%y(3)/hubcst, f(20)*uni%y(3)*uni%y(2)/hubcst,
     |      f(31)*uni%y(3)*uni%y(5)/hubcst, f(13)*uni%y(3)*uni%y(1)/hubcst, f(33)*uni%y(3)*uni%y(8)/hubcst
c      3He creation: SKM Figure 6a
         write(14,"(5(e13.5,' '))") uni%t9, f(28)*uni%y(3)*uni%y(3)/hubcst, r(16)*uni%y(4)*uni%y(2)/hubcst,
     |      f(20)*uni%y(3)*uni%y(2)/hubcst, f(23)*uni%y(7)*uni%y(2)/hubcst
c      3He destruction: SKM Figure 6b
         write(15,"(7(e13.5,' '))") uni%t9, f(16)*uni%y(5)*uni%y(1)/hubcst, f(31)*uni%y(5)*uni%y(3)/hubcst,
     |      f(14)*uni%y(5)*uni%y(1)/hubcst, f(32)*uni%y(5)*uni%y(5)/hubcst, f(27)*uni%y(5)*uni%y(6)/hubcst,
     |      r(20)*uni%y(5)/hubcst
c     Li6 creation
         write(18,"(8(e13.5,' '))") uni%t9, r(15)*uni%y(8)/hubcst, r(18)*uni%y(4)*uni%y(6)/hubcst,
     |      r(22)*uni%y(9)/hubcst, r(23)*uni%y(5)*uni%y(6)/hubcst, f(25)*uni%y(3)*uni%y(6)/hubcst,
     |      f(46)*uni%y(12)*uni%y(2)/hubcst, r(49)*uni%y(13)/hubcst
c     Li6 destruction
         write(19,"(8(e13.5,' '))") uni%t9, f(15)*uni%y(7)*uni%y(1)/hubcst, f(18)*uni%y(7)*uni%y(1)/hubcst,
     |      f(22)*uni%y(7)*uni%y(2)/hubcst, f(23)*uni%y(7)*uni%y(2)/hubcst, r(25)*uni%y(7)*uni%y(2)/hubcst,
     |      r(46)*uni%y(7)*uni%y(6)/hubcst, f(49)*uni%y(7)*uni%y(6)/hubcst
c     Li7 creation
         write(20,"(5(e13.5,' '))") uni%t9, f(26)*uni%y(4)*uni%y(6)/hubcst, f(17)*uni%y(1)*uni%y(9)/hubcst,
     |      r(35)*uni%y(10)/hubcst, f(15)*uni%y(7)*uni%y(1)/hubcst
c     Li7 destruction
         write(21,"(6(e13.5,' '))") uni%t9, f(24)*uni%y(8)*uni%y(2)/hubcst, f(60)*uni%y(8)*uni%y(3)/hubcst,
     |      f(35)*uni%y(8)*uni%y(1)/hubcst, r(26)*uni%y(8)/hubcst, f(50)*uni%y(8)*uni%y(6)/hubcst
c     Be7 creation
         write(22,"(5(e13.5,' '))") uni%t9, f(27)*uni%y(5)*uni%y(6)/hubcst, r(17)*uni%y(8)*uni%y(2)/hubcst,
     |      f(22)*uni%y(7)*uni%y(2)/hubcst, r(40)*uni%y(11)/hubcst
c     Be7 destruction
         write(23,"(6(e13.5,' '))") uni%t9, f(17)*uni%y(9)*uni%y(1)/hubcst, f(19)*uni%y(9)*uni%y(1)/hubcst,
     |      r(27)*uni%y(9)/hubcst, f(34)*uni%y(9)*uni%y(3)/hubcst, f(40)*uni%y(9)*uni%y(2)/hubcst

!        write(4,200) t,hubcst,scalarfield%scal,scalarfield%dscal,ScalarPotential(scalarfield%scal)
      END IF

      IF (itime.eq.8) THEN         !Right after a run.
        xout(it,8) = xout(it,8) + xout(it,9)  !Add beryllium to lithium (electron-capture, half-life = 53 days)
        xout(it,5) = xout(it,5) + xout(it,4)  !Add tritium to helium-3 (beta-dacay, half-life = 12.33 y)
        !xout(it,6) = xout(it,6)-0.0025
                  !Radiative, coulomb, finite-temperature corrections (Ref 1).
        write(3,200) etaout(it),xout(it,3),
     |                xout(it,5),xout(it,6),xout(it,8)  
     |                             !Output eta, H2, He3, He4, and Li7.
 200    FORMAT (5(e13.5,' '))
      END IF

!C30-----close FILE--------------------------------------
      IF (.false. .and. itime.eq.9) THEN        !End of program.
        CLOSE (unit=4)
        CLOSE (unit=10)
        CLOSE (unit=11)
        CLOSE (unit=12)
        CLOSE (unit=13)
        CLOSE (unit=14)
        CLOSE (unit=15)
        CLOSE (unit=16)
        CLOSE (unit=17)
        CLOSE (unit=18)
        CLOSE (unit=19)
        CLOSE (unit=20)
        CLOSE (unit=21)
        CLOSE (unit=22)
        CLOSE (unit=23)
      END IF

      IF (itime.eq.10) THEN        !End of program.
        CLOSE (unit=3)
      END IF
      RETURN

C-------REFERENCES--------------------------------------
C     1) D.A. Dicus, E.W. Kolb, A.M. Gleeson, E.C.G. Sudarshan, V.L. Teplitz,
C        M.S. Turner, Phys. Rev. D., 26,2694 (1982).

      END
