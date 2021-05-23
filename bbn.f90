PROGRAM bbn

    use physics_parameters
    use computational_parameters
    use reactions_module
    use variables
    
    INTEGER :: i, j
    REAL :: ln_mq

!..........SET OUTPUT OPTION TO DEFAULT.
    nout    = 0                  !No output requests.

    ! Run options:
    ! 1 = full run, 2 = abridged network, 3 = small run
    CALL ResetComputationalParameters(1)
    CALL ResetPhysicsParameters
    CALL InitialiseReactions

    itime = 1       !Time = beginning of program.
    CALL check
    itime = 2       !Time = beginning of run section.
    CALL check

    cy = 0.1
    ct = 0.01
    inc = 10

!    CALL GetLinearResponse
!    CALL VaryLightQuarkMass

    if(.true.) then
        call InitialiseReactions
!        call SetDeutronBindingVariation(-0.538)
!        CALL SetHe5ResonancePositions(0.1)

        ln_mq = -0.04
        !CALL SetDeutronBindingVariation(ln_mq * (-1.39) * 25.8)   ! AV18
        !CALL SetHe5ResonancePositions(ln_mq * 17.59) ! delta E_r in MeV
        !CALL SetLi5ResonancePositions(ln_mq * 1.43) ! delta E_r in MeV
        !CALL SetHe5ResonancePositions(ln_mq * 0.81) ! delta E_r in MeV

        itime = 3       !Time = before computation
        CALL check

        CALL DriverRK2        !Do nucleosynthesis computation.

        itime = 8       !Time = after computation
        CALL check

        print *, it
        write (*, "(1p,8(e14.5))") t9out(it), (xout(it,3)), (xout(it,4)), (xout(it,5)), &
                                  (xout(it,6)), (xout(it,7)), (xout(it,8)), (xout(it,9))
    end if

    itime = 9       !Time = end of run section.
    CALL check
    
    ! Output file

    OPEN (unit=2, file='nuc123.dat', status='unknown')  !Output file.
    
    nout = nout + 1            !Keep track of number of output requests.
    IF (nout.eq.1) THEN
        write(2, "(54x,'NUCLIDE ABUNDANCE YIELDS',/,54x,'------- --------- ------',//)")
    END IF
    
    write (2, "(' Computational parameters:',/, &
              & '   cy = ',f5.3,'/  ct = ',f5.3, &
              & '/  initial temp = ',1pe8.2, &
              & '/  final temp = ',1pe8.2, &
              & '/  smallest abundances allowed = ',1pe8.2)") &
        cy, ct, t9i, t9f, ytmin
    
    write (2, "(' Model parameters:',/, &
              & '   g = ',f5.2,'/  tau = ',f6.2, &
              & '/  # nu = ',f5.2,'/  lambda = ',1pe10.3, &
              & '/  xi-e = ',e10.3,'/  xi-m = ',e10.3, &
              & '/  xi-t = ',e10.3,/)") &
        c(1), c(2), c(3), cosmo, xi(1), xi(2), xi(3)
    
    !PRINT HEADINGS, ABUNDANCES FOR NEUTRON TO LI8.         
    write (2, "(4x,'Temp',8x,'N/H',10x,'P',10x,'D/H',9x,'T/H',8x, &
              & 'He3/H',8x,'He4',8x,'Li6/H',7x,'Li7/H',7x, &
              & 'Be7/H',6x,'Li8/H&up',/,132('-'))")
    DO j = 1, it
        write(2,"(1pe10.3,1p10e12.3)") t9out(j), (xout(j,i),i=1,10)
    END DO

    !PRINT THERMODYNAMIC QUANTITIES.         
    write (2, "(' ',/,4x,'Temp',9x,'T',10x,'rhog',8x,'rhoe',7x, &
              & 'rhone',8x,'rhob',8x,'phie',9x,'dt',9x, &
              & 'eta',10x,'H',/,132('-'))")
    DO j = 1, it
        write (2, "(1pe10.3,9e12.3)") t9out(j), tout(j), (thmout(j,i),i=1,5), &
            dtout(j), etaout(j),hubout(j)
    END DO
    write (2, "(///)")
    
    close (unit=2,status='keep')   !Close output file.

    itime = 10      !Time = end of program.
    CALL check

END PROGRAM bbn

SUBROUTINE GetLinearResponse
    ! Get linear response to variation of a parameter X
    !   d(lnY_i)/d(lnX),
    ! where i is (D, 3He, 4He, 6Li, 7Li)
    use scalarfield_module
    use computational_parameters
    use reactions_module
    use variables

    INTEGER :: i, j
!    REAL :: lnvar, lnvar_step
    REAL:: lnvar_step
    INTEGER, PARAMETER :: numsteps = 5
    REAL, DIMENSION(numsteps, 6) :: lnY
    REAL :: dlnYdlnX, X0, binding_energy

    OPEN (unit=7, file='bbn.dat', status='unknown')
    write (7, "('eta  var  D  3He  4He  6Li  7Li  7Be  7Li+7Be')")

    lnvar = -0.02
    lnvar_step = 2.0 * abs(lnvar) / (numsteps - 1)
!    lnvar = -1.82
!    lnvar_step = 0.2

    X0 = 1.0

!    do j = 10, 80, 5
!        eta = j/10.0 * 1.e-10
!        print *, eta
!        lnvar = -0.18

    DO i = 1, numsteps
        ! Reset to SBBN values
        CALL InitialiseReactions

!        CALL SetDeutronBindingVariation(lnvar)
!        CALL SetDeutronBindingVariation(lnvar * (-1.39) * 25.8)   ! AV18
!        CALL SetBindingVariation(4, lnvar * (-1.44) * 98.39)  ! 3H = T
!        CALL SetBindingVariation(5, lnvar * (-1.55) * 89.53)  ! 3He
!        CALL SetBindingVariation(6, lnvar * (-1.08) * 328.2)  ! 4He
!        CALL SetBindingVariation(7, lnvar * (-1.36) * 371.1)  ! 6Li
!        CALL SetBindingVariation(8, lnvar * (-1.50) * 455.2)  ! 7Li
!        CALL SetBindingVariation(9, lnvar * (-1.57) * 436.2)  ! 7Be
!        do j = 10, 11
!            binding_energy = (zm(j)*km(2) + (am(j)-zm(j))*km(1) - km(j))/1000./kB
!            CALL SetBindingVariation(j, lnvar * (-1.45) * binding_energy)
!        end do
!        CALL SetBindingVariation(12, lnvar * (-1.59) * 674.7) ! 9Be
!        do j = 13, 26
!            binding_energy = (zm(j)*km(2) + (am(j)-zm(j))*km(1) - km(j))/1000./kB
!            CALL SetBindingVariation(j, lnvar * (-1.45) * binding_energy)
!        end do

!        CALL SetBindingVariation(4, lnvar)
!        c(2) = 889.0 * (1.0 + lnvar)

!        CALL SetLi5ResonancePositions(lnvar * 7.32) ! delta E_r in MeV
!        CALL SetHe5ResonancePositions(lnvar * 8.29) ! delta E_r in MeV
!        CALL SetLi5ResonanceWidth(lnvar) ! % change in Gamma
!        CALL SetHe5ResonanceWidth(lnvar) !

        itime = 3       !Time = before computation
        CALL check
        
        CALL DriverRK4        !Do nucleosynthesis computation.

        itime = 8       !Time = after computation
        CALL check

        ! dQ, T9, D, 3He, 4He, 7Li, 7Be
        write (*, "(1p,7(e14.5))") lnvar, t9out(it), xout(it,3), xout(it,5), xout(it,6), xout(it,8), xout(it,9)

        ! eta, dQ, D, 3He, 4He, 6Li, 7Li, 7Be, 7Li+7Be
        write (7, "(1p,9(e14.5))") eta, lnvar, (xout(it,3)), (xout(it,5)), &
                                  (xout(it,6)), (xout(it,7)), (xout(it,8)), (xout(it,9)), &
                                  (xout(it,8) + xout(it,9))
!        write (7, "(1p9(e14.5))") eta*1.e10, lnvar, xout(it,3)*1.e5, xout(it,5)*1.e6, xout(it,6), &
!                                  xout(it,7)*1.e14, xout(it,8)*1.e10, xout(it,9)*1.e10, &
!                                  (xout(it,8) + xout(it,9))*1.e10

        lnY(i, 1) = log(xout(it,3))
        lnY(i, 2) = log(xout(it,5))
        lnY(i, 3) = log(xout(it,6))
        lnY(i, 4) = log(xout(it,7))
        lnY(i, 5) = log(xout(it,8))
        lnY(i, 6) = log(xout(it,9))

        lnvar = lnvar + lnvar_step
    END DO
!    end do
    
    CLOSE (unit=7, status='keep')

    if(numsteps .eq. 5) then
    DO i = 1, 6
        dlnYdlnX = (2./3.*(lnY(4,i) - lnY(2,i)) - 1./12.*(lnY(5,i) - lnY(1,i)))/lnvar_step * X0
        write(*, "(i2, ' dlnY/dlnX = ', f8.3, ' (', f7.5,')'))") &
            i, dlnYdlnX, abs((lnY(5,i) - lnY(3, i))/2./lnvar_step*X0 - dlnYdlnX) 
    END DO
    end if
END SUBROUTINE

SUBROUTINE VaryLightQuarkMass
    ! Get response to variation of light quark mass, m_q
    !   d(ln Y_i)/d(ln m_q),
    ! where i is (D, 3He, 4He, 6Li, 7Li)
    use physics_parameters
    use scalarfield_module
    use computational_parameters
    use reactions_module
    use variables

    INTEGER :: i, j
    REAL :: ln_mq, ln_mq_step   ! d(m_q)/m_q and its step size
    INTEGER, PARAMETER :: numsteps = 2
    REAL :: binding_energy

    OPEN (unit=7, file='bbn.dat', status='unknown')
    !write (7, "('m_q, T9, D, 3He, 4He, 6Li, 7Li')")

    ln_mq = -0.01
    ln_mq_step = 0.04 / (numsteps - 1)

    DO i = 1, numsteps
        ! Reset to SBBN values
        CALL InitialiseReactions

        !CALL SetDeutronBindingVariation(ln_mq * (-0.75) * 25.8)   ! AV28
        CALL SetDeutronBindingVariation(ln_mq * (-1.39) * 25.8)   ! AV18
        CALL SetBindingVariation(4, ln_mq * (-1.44) * 98.39)  ! 3H = T
        CALL SetBindingVariation(5, ln_mq * (-1.55) * 89.53)  ! 3He
        CALL SetBindingVariation(6, ln_mq * (-1.08) * 328.2)  ! 4He
        CALL SetBindingVariation(7, ln_mq * (-1.36) * 371.1)  ! 6Li
        CALL SetBindingVariation(8, ln_mq * (-1.50) * 455.2)  ! 7Li
        CALL SetBindingVariation(9, ln_mq * (-1.57) * 436.2)  ! 7Be
        do j = 10, 11
            binding_energy = (zm(j)*km(2) + (am(j)-zm(j))*km(1) - km(j))/1000./kB
            CALL SetBindingVariation(j, ln_mq * (-1.45) * binding_energy)
        end do
        CALL SetBindingVariation(12, ln_mq * (-1.59) * 674.7) ! 9Be
        do j = 13, 26
            binding_energy = (zm(j)*km(2) + (am(j)-zm(j))*km(1) - km(j))/1000./kB
            CALL SetBindingVariation(j, ln_mq * (-1.45) * binding_energy)
        end do

!        CALL SetLi5ResonancePositions(ln_mq * 17.59) ! delta E_r in MeV
!        CALL SetHe5ResonancePositions(ln_mq * 18.68) ! delta E_r in MeV
        CALL SetLi5ResonancePositions(ln_mq * 7.32) ! delta E_r in MeV
        CALL SetHe5ResonancePositions(ln_mq * 8.29) ! delta E_r in MeV

!        CALL PrintQValues

        itime = 3       !Time = before computation
        CALL check
        
        CALL DriverRK4        !Do nucleosynthesis computation.

        itime = 8       !Time = after computation
        CALL check

        ! d(m_q)/m_q, T9, D, 3He, 4He, 6Li, 7Li
        write (7, "(1p,7(e14.5))") ln_mq, (xout(it,3)), (xout(it,5)), &
                                  (xout(it,6)), (xout(it,7)), (xout(it,8))

        write (*, "(1p,9(e14.5))") ln_mq, t9out(it), (xout(it,3)), (xout(it,4)), (xout(it,5)), &
                                  (xout(it,6)), (xout(it,7)), (xout(it,8)), (xout(it,9))

        ln_mq = ln_mq + ln_mq_step
    END DO
    write (7, "('m_q, D, 3He+T, 4He, 6Li, 7Li+7Be')")

    CLOSE (unit=7, status='keep')

END SUBROUTINE
