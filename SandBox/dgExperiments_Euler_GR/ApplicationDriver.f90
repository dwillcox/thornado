PROGRAM ApplicationDriver

  USE KindModule, ONLY: &
    DP, Zero, One, Two, Pi, Four, TwoPi
  USE ProgramInitializationModule, ONLY: &
    InitializeProgram, &
    FinalizeProgram
  USE ReferenceElementModuleX, ONLY: &
    InitializeReferenceElementX, &
    FinalizeReferenceElementX
  USE ReferenceElementModuleX_Lagrange, ONLY: &
    InitializeReferenceElementX_Lagrange, &
    FinalizeReferenceElementX_Lagrange
  USE EquationOfStateModule, ONLY: &
    InitializeEquationOfState, &
    FinalizeEquationOfState
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1
  USE GeometryComputationModule, ONLY: &
    ComputeGeometryX
  USE InitializationModule_GR, ONLY: &
    InitializeFields_GR, ReadParameters
  USE Euler_GR_SlopeLimiterModule, ONLY: &
    Euler_GR_InitializeSlopeLimiter, &
    Euler_GR_FinalizeSlopeLimiter, &
    Euler_GR_ApplySlopeLimiter
  USE Euler_GR_PositivityLimiterModule, ONLY: &
    Euler_GR_InitializePositivityLimiter, &
    Euler_GR_FinalizePositivityLimiter, &
    Euler_GR_ApplyPositivityLimiter
  USE Euler_GR_UtilitiesModule, ONLY: &
    ComputeFromConserved_GR, &
    ComputeTimeStep_GR
  USE InputOutputModuleHDF, ONLY: &
    WriteFieldsHDF
  USE FluidFieldsModule, ONLY: &
    uCF, uPF, uAF
  USE GeometryFieldsModule, ONLY: &
    uGF
  USE Euler_GR_dgDiscretizationModule, ONLY: &
    Euler_GR_ComputeIncrement_DG_Explicit
  USE TimeSteppingModule_SSPRK, ONLY: &
    InitializeFluid_SSPRK, &
    FinalizeFluid_SSPRK, &
    UpdateFluid_SSPRK
  USE UnitsModule, ONLY: &
    Millisecond
  USE Euler_GR_TallyModule, ONLY: &
    Euler_GR_InitializeTally, &
    Euler_GR_FinalizeTally, &
    Euler_GR_ComputeTally

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  CHARACTER(32) :: ProgramName
  CHARACTER(32) :: RiemannProblemName, RiemannProblem2dName
  CHARACTER(32) :: SphericalRiemannProblemName
  CHARACTER(32) :: CoordinateSystem
  LOGICAL       :: wrt
  LOGICAL       :: UseFixed_dt, UseSourceTerm
  LOGICAL       :: UseSlopeLimiter
  LOGICAL       :: UseCharacteristicLimiting
  LOGICAL       :: UseTroubledCellIndicator
  LOGICAL       :: UsePositivityLimiter
  INTEGER       :: iCycle, iCycleD, iCycleW = 0
  INTEGER       :: nX(3), bcX(3), swX(3), nNodes
  INTEGER       :: nStagesSSPRK
  REAL(DP)      :: SlopeTolerance
  REAL(DP)      :: Min_1, Min_2
  REAL(DP)      :: xL(3), xR(3), Gamma
  REAL(DP)      :: t, dt, t_end, dt_wrt, t_wrt, wTime, CFL
  REAL(DP)      :: BetaTVD, BetaTVB
  REAL(DP)      :: LimiterThresholdParameter

  ! --- Sedov blast wave ---
  REAL(DP) :: Eblast
  INTEGER  :: nDetCells
  REAL(DP) :: Vmax, LorentzFactor

  ! --- Standing accretion shock ---
  REAL(DP), ALLOCATABLE :: FluidFieldParameters(:)
  REAL(DP)              :: M_PNS = Zero, Ri, R_PNS, R_shock, Rf

  LOGICAL :: DEBUG = .FALSE.
  LOGICAL :: WriteGF = .FALSE., WriteFF = .TRUE.
  REAL(DP) :: wTime_UF, wTime_CTS, Timer_Evolution
  
!  ProgramName = 'RiemannProblem'
  ProgramName = 'RiemannProblem2d'
!  ProgramName = 'SphericalRiemannProblem'
!  ProgramName = 'SphericalSedov'
!  ProgramName = 'KelvinHelmholtz_Relativistic'
!  ProgramName = 'KelvinHelmholtz'
!  ProgramName = 'StandingAccretionShock'

  SELECT CASE ( TRIM( ProgramName ) )

    CASE( 'RiemannProblem' )

      RiemannProblemName = 'Sod'

      SELECT CASE ( TRIM( RiemannProblemName ) )

        CASE( 'Sod' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'MBProblem1' )
          Gamma = 4.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'MBProblem4' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'PerturbedShockTube' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.35d0
          bcX   = [ 2, 0, 0 ]

        CASE( 'CartesianSedov' )
          Gamma = 4.0_DP / 3.0_DP
          t_end = 0.2d0
          bcX   = [ 32, 0, 0 ]
          nDetCells = 1
          Eblast    = 1.0d-3

        CASE( 'ShockReflection' )
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.75d0
          bcX   = [ 23, 0, 0 ]

        CASE DEFAULT ! Sod
          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 0, 0 ]

      END SELECT

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 256, 1, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'RiemannProblem2d' )

      RiemannProblem2dName = 'DzB2002'

      SELECT CASE ( TRIM( RiemannProblem2dName ) )

        CASE( 'DzB2002' )

          Gamma = 5.0_DP / 3.0_DP
          t_end = 0.4d0
          bcX   = [ 2, 2, 0 ]

      END SELECT

      CoordinateSystem = 'CARTESIAN'

      nX  = [ 64, 64, 1 ]
      xL  = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR  = [ 1.0_DP, 1.0_DP, 1.0_DP ]

    CASE( 'SphericalRiemannProblem' )

      SphericalRiemannProblemName = 'SphericalSod'

      CoordinateSystem = 'SPHERICAL'

      Gamma = 5.0_DP / 3.0_DP

      nX = [ 128, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 2.0_DP, Pi, TwoPi ]

      bcX = [ 2, 0, 0 ]

      t_end = 5.0d-1

      WriteGF = .TRUE.

    CASE( 'SphericalSedov' )

      nDetCells = 1
      Eblast    = 1.0d-3

      CoordinateSystem = 'SPHERICAL'

      Gamma = 1.4_DP!4.0_DP / 3.0_DP

      nX = [ 256, 1, 1 ]
      xL = [ 0.0_DP, 0.0_DP, 0.0_DP ]
      xR = [ 1.2_DP, Pi, TwoPi ]

      bcX = [ 3, 0, 0 ]

      t_end = 1.0d0

      WriteGF = .TRUE.

    CASE( 'KelvinHelmholtz_Relativistic' )

       CoordinateSystem = 'CARTESIAN'

       Gamma = 4.0d0 / 3.0d0

       nX = [ 16, 32, 1 ]
       xL = [ -0.5d0, -1.0d0, 0.0d0 ]
       xR = [  0.5d0,  1.0d0, 1.0d0 ]

       bcX = [ 1, 1, 0 ]

       t_end = 3.0d0

    CASE( 'KelvinHelmholtz' )

       CoordinateSystem = 'CARTESIAN'

       Gamma = 5.0d0 / 3.0d0

       nX = [ 32, 32, 1 ]
       xL = [ 0.0d0, 0.0d0, 0.0d0 ]
       xR = [ 1.0d0, 1.0d0, 1.0d0 ]

       bcX = [ 1, 1, 0 ]

       t_end = 1.5d0

    CASE( 'StandingAccretionShock' )

      CoordinateSystem = 'SPHERICAL'

      CALL ReadParameters &
             ( '../StandingAccretionShock_Parameters.dat', &
                 FluidFieldParameters )

      M_PNS   = FluidFieldParameters(1)
      Gamma   = FluidFieldParameters(2)
      Ri      = FluidFieldParameters(3)
      R_PNS   = FluidFieldParameters(4)
      R_shock = FluidFieldParameters(5)

      nX = [ 256, 1, 1 ]
      xL = [ R_PNS, 0.0_DP, 0.0_DP ]
      xR = [ Two * R_shock, Pi, Four ]

      bcX = [ 11, 0, 0 ]

      t_end = 3.0d2 * Millisecond

      WriteGF = .TRUE.

    CASE DEFAULT

      WRITE(*,*)
      WRITE(*,'(A21,A)') 'Invalid ProgramName: ', ProgramName
      WRITE(*,'(A)')     'Valid choices:'
      WRITE(*,'(A)')     '  RiemannProblem'
      WRITE(*,'(A)')     '  RiemannProblem2d'
      WRITE(*,'(A)')     '  SphericalRiemannProblem'
      WRITE(*,'(A)')     '  SphericalSedov'
      WRITE(*,'(A)')     '  KelvinHelmholtz2D_Relativistic'
      WRITE(*,'(A)')     '  StandingAccretionShock'
      WRITE(*,'(A)')     'Stopping...'
      STOP

  END SELECT

  nNodes = 3
  IF( .NOT. nNodes .LE. 4 ) &
    STOP 'nNodes must be less than or equal to four.'

  BetaTVD = 1.75d0
  BetaTVB = 0.0d0

  UseSlopeLimiter           = .TRUE.
  SlopeTolerance            = 1.0d-6
  UseCharacteristicLimiting = .TRUE.

  UseTroubledCellIndicator  = .TRUE.
  LimiterThresholdParameter = 0.015_DP

  UsePositivityLimiter = .TRUE.
  Min_1 = 1.0d-13
  Min_2 = 1.0d-13

  UseFixed_dt   = .FALSE.
  UseSourceTerm = .TRUE.

  iCycleD = 10
!!$  iCycleW = 1000; dt_wrt = -1.0d0
  dt_wrt = 1.0d-2 * t_end; iCycleW = -1

  IF( dt_wrt .GT. Zero .AND. iCycleW .GT. 0 ) &
    STOP 'dt_wrt and iCycleW cannot both be present'

  nStagesSSPRK = 3
  IF( .NOT. nStagesSSPRK .LE. 3 ) &
    STOP 'nStagesSSPRK must be less than or equal to three.'

  ! --- Cockburn & Shu, (2001), JSC, 16, 173 ---
  CFL = 0.5d0 / ( 2.0d0 * DBLE( nNodes - 1 ) + 1.0d0 )

  CALL WriteProgramHeader

  CALL InitializeProgram &
         ( ProgramName_Option &
             = TRIM( ProgramName ), &
           nX_Option &
             = nX, &
           swX_Option &
             = [ 1, 1, 1 ], &
           bcX_Option &
             = bcX, &
           xL_Option &
             = xL, &
           xR_Option &
             = xR, &
           nNodes_Option &
             = nNodes, &
           CoordinateSystem_Option &
             = TRIM( CoordinateSystem ), &
           BasicInitialization_Option &
             = .TRUE. )

  CALL InitializeReferenceElementX

  CALL InitializeReferenceElementX_Lagrange

  CALL ComputeGeometryX &
         ( iX_B0, iX_E0, iX_B1, iX_E1, uGF, Mass_Option = M_PNS )

  CALL InitializeEquationOfState &
         ( EquationOfState_Option = 'IDEAL', &
           Gamma_IDEAL_Option = Gamma )

  CALL Euler_GR_InitializeSlopeLimiter &
         ( BetaTVD_Option = BetaTVD, &
           BetaTVB_Option = BetaTVB, &
           SlopeTolerance_Option &
             = SlopeTolerance, &
           UseSlopeLimiter_Option &
             = UseSlopeLimiter, &
           UseCharacteristicLimiting_Option &
             = UseCharacteristicLimiting, &
           UseTroubledCellIndicator_Option &
             = UseTroubledCellIndicator, &
           LimiterThresholdParameter_Option &
             = LimiterThresholdParameter )

  CALL Euler_GR_InitializePositivityLimiter &
         ( Min_1_Option = Min_1, &
           Min_2_Option = Min_2, &
           UsePositivityLimiter_Option = UsePositivityLimiter )

  CALL InitializeFields_GR &
         ( RiemannProblemName_Option = TRIM( RiemannProblemName ), &
           RiemannProblem2dName_Option = TRIM( RiemannProblem2dName ), &
           SphericalRiemannProblemName_Option &
             = TRIM( SphericalRiemannProblemName ), &
           nDetCells_Option = nDetCells, Eblast_Option = Eblast )

  CALL Euler_GR_ApplySlopeLimiter &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:),&
           uCF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:) )

  CALL Euler_GR_ApplyPositivityLimiter &
         ( iX_B0, iX_E0, iX_B1, iX_E1, &
           uGF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:),&
           uCF(:,iX_B1(1):iX_E1(1),iX_B1(2):iX_E1(2),iX_B1(3):iX_E1(3),:) )

  CALL WriteFieldsHDF &
       ( 0.0_DP, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  CALL InitializeFluid_SSPRK( nStages = nStagesSSPRK )

  WRITE(*,*)
  WRITE(*,'(A2,A,ES11.3E3)') '', 'CFL: ', CFL
  WRITE(*,*)
  WRITE(*,'(A2,A)') '', 'Evolving Fields...'
  WRITE(*,*)

  t     = 0.0_DP
  t_wrt = dt_wrt
  wrt   = .FALSE.

  IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL Euler_GR_InitializeTally'
  CALL Euler_GR_InitializeTally &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  iCycle = 0
  Timer_Evolution = MPI_WTIME()
  DO WHILE( t .LT. t_end )

    iCycle = iCycle + 1

    IF( .NOT. UseFixed_dt )THEN
      IF( DEBUG ) wTime_CTS = MPI_WTIME( )
      IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL ComputeTimeStep_GR'
      CALL ComputeTimeStep_GR &
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               CFL = CFL, TimeStep = dt, UseSourceTerm_Option = UseSourceTerm )
      IF( DEBUG ) wTime_CTS = MPI_WTIME( ) - wTime_CTS
    ELSE
      dt = CFL * ( xR(1) - xL(1) ) / DBLE( nX(1) )
    END IF

    IF( t + dt .LT. t_end )THEN
      t = t + dt
    ELSE
      dt = t_end - t
      t  = t_end
    END IF

    IF( MOD( iCycle, iCycleD ) .EQ. 0 )THEN

      IF( ProgramName .EQ. 'StandingAccretionShock' )THEN
        WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A4,A5,ES13.6E3,A3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t / Millisecond, ' ms ', &
          'dt = ', dt / Millisecond, ' ms'
      ELSE
        WRITE(*,'(A8,A8,I8.8,A2,A4,ES13.6E3,A1,A5,ES13.6E3)') &
          '', 'Cycle = ', iCycle, '', 't = ',  t, '', 'dt = ', dt
      END IF

      IF( DEBUG )THEN
        IF( .NOT. UseFixed_dt ) &
          WRITE(*,'(A8,A,ES10.3E3,A)') &
            '', 'Time to call ComputeTimeStep:   ', wTime_CTS, ' s'
        WRITE(*,'(A8,A,ES10.3E3,A)') &
          '', 'Time to call UpdateFluid_SSPRK: ', wTime_UF,  ' s'
      END IF

    END IF

    IF( DEBUG ) wTime_UF = MPI_WTIME()
    IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL UpdateFluid_SSPRK'
    CALL UpdateFluid_SSPRK &
           ( t, dt, uGF, uCF, Euler_GR_ComputeIncrement_DG_Explicit )
    IF( DEBUG ) wTime_UF = MPI_WTIME() - wTime_UF

    IF( DEBUG )THEN
      WRITE(*,'(A)') 'AD: CALL ComputeFromConserved_GR'
      CALL ComputeFromConserved_GR &
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

      Vmax = MAXVAL( ABS( &
               uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),2)))
      LorentzFactor = One / SQRT( One - Vmax**2 )

      IF( MOD( iCycle, iCycleD ) .EQ. 0 )THEN
        WRITE(*,'(A8,A4,F12.10)') '', 'V = ', Vmax
        WRITE(*,'(A8,A4,F10.5)')  '', 'W = ', LorentzFactor
        WRITE(*,*)
      END IF
   END IF

    IF( iCycleW .GT. 0 )THEN
      IF( MOD( iCycle, iCycleW ) .EQ. 0 ) &
        wrt = .TRUE.
    ELSE
      IF( t + dt .GT. t_wrt )THEN
        t_wrt = t_wrt + dt_wrt
        wrt   = .TRUE.
      END IF
    END IF

    IF( wrt )THEN

      IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL ComputeFromConserved_GR'
      CALL ComputeFromConserved_GR &
             ( iX_B0, iX_E0, &
               uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
               uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

      IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL WriteFieldsHDF'
      CALL WriteFieldsHDF &
             ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

      IF( DEBUG ) WRITE(*,'(A)') 'AD: CALL Euler_GR_ComputeTally'
      CALL Euler_GR_ComputeTally &
           ( iX_B0, iX_E0, &
             uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
             Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

      wrt = .FALSE.

    END IF

  END DO
  Timer_Evolution = MPI_WTIME() - Timer_Evolution
  WRITE(*,*)
  WRITE(*,'(A,ES13.6E3)') 'Total evolution time: ', Timer_Evolution

  CALL ComputeFromConserved_GR &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uPF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uAF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:) )

  CALL WriteFieldsHDF &
         ( t, WriteGF_Option = WriteGF, WriteFF_Option = WriteFF )

  CALL Euler_GR_ComputeTally &
         ( iX_B0, iX_E0, &
           uGF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           uCF(:,iX_B0(1):iX_E0(1),iX_B0(2):iX_E0(2),iX_B0(3):iX_E0(3),:), &
           Time = t, iState_Option = 1, DisplayTally_Option = .TRUE. )

  CALL Euler_GR_FinalizePositivityLimiter

  CALL Euler_GR_FinalizeSlopeLimiter

  CALL FinalizeFluid_SSPRK

  CALL FinalizeReferenceElementX 

  CALL FinalizeReferenceElementX_Lagrange

  CALL FinalizeEquationOfState

  CALL FinalizeProgram

CONTAINS

  SUBROUTINE WriteProgramHeader

    IF( ( .NOT. nNodes .GT. 1 ) .AND. UseSlopeLimiter )THEN
      WRITE(*,*)
      WRITE(*,'(A)') 'Slope limiter requires nNodes > 1'
      WRITE(*,'(A)') 'Setting UseSlopeLimiter to .FALSE.'
      WRITE(*,'(A)') 'Setting UseCharacteristicLimiting to .FALSE.'
      WRITE(*,'(A)') 'Setting UseTroubledCellIndicator to .FALSE.'
      WRITE(*,*)
      UseSlopeLimiter           = .FALSE.
      UseCharacteristicLimiting = .FALSE.
      UseTroubledCellIndicator  = .FALSE.
    END IF
    IF( ( .NOT. nNodes .GT. 1 ) .AND. UsePositivityLimiter )THEN
      WRITE(*,*)
      WRITE(*,'(A)') 'Positivity limiter requires nNodes > 1'
      WRITE(*,'(A)') 'Setting UsePositivityLimiter to .FALSE.'
      UsePositivityLimiter = .FALSE.
      WRITE(*,*)
    END IF
    IF( .NOT. UseSlopeLimiter .AND. UseCharacteristicLimiting )THEN
      WRITE(*,*)
      WRITE(*,'(A)') 'Characteristic limiting requires use of slope limiter'
      WRITE(*,'(A)') 'Setting UseCharacteristicLimiting to .FALSE.'
      WRITE(*,*)
      UseCharacteristicLimiting = .FALSE.
    END IF
    IF( .NOT. UseSlopeLimiter .AND. UseTroubledCellIndicator )THEN
      WRITE(*,*)
      WRITE(*,'(A)') 'Troubled cell indicator requires use of slope limiter'
      WRITE(*,'(A)') 'Setting UseTroubledCellIndicator to .FALSE.'
      WRITE(*,*)
      UseTroubledCellIndicator = .FALSE.
    END IF
    IF( TRIM(CoordinateSystem) .EQ. 'CARTESIAN' .AND. UseSourceTerm )THEN
      WRITE(*,*)
      WRITE(*,'(A)') 'Cartesian coordinates demand that UseSourceTerm = .FALSE.'
      WRITE(*,'(A)') 'Setting UseSourceTerm to .FALSE.'
      WRITE(*,*)
      UseSourceTerm = .FALSE.
    END IF

    ! --- Write program parameters to header file ---
    OPEN( 100, FILE = '../Output/.ProgramHeader' )
      WRITE(100,'(A,A)')         'Program Name: ', TRIM(ProgramName)
      WRITE(100,*)
      IF( TRIM( ProgramName ) .EQ. 'RiemannProblem' ) &
        WRITE(100,'(A,A)') &
          'Riemann Problem Name: ', RiemannProblemName
      IF( TRIM( ProgramName ) .EQ. 'RiemannProblem2d' ) &
        WRITE(100,'(A,A)') &
          '2D Riemann Problem Name: ', RiemannProblem2dName
      IF( TRIM( ProgramName ) .EQ. 'SphericalSedov' )THEN
        WRITE(100,'(A,I4.4)')     'nDetCells: ', nDetCells
        WRITE(100,'(A,ES10.3E3)') 'Eblast:    ', Eblast
      END IF
      IF( TRIM( ProgramName ) .EQ. 'StandingAccretionShock' )THEN
        WRITE(100,'(A,ES10.3E3)') 'PNS Mass:     ', M_PNS
        WRITE(100,'(A,ES10.3E3)') 'Inner radius: ', Ri
        WRITE(100,'(A,ES10.3E3)') 'PNS Radius:   ', R_PNS
        WRITE(100,'(A,ES10.3E3)') 'Shock Radius: ', R_shock
      END IF
      WRITE(100,*)
      WRITE(100,'(A,F5.3)')      'Gamma_IDEAL: ', Gamma
      WRITE(100,*)
      WRITE(100,'(A)')           'Mesh'
      WRITE(100,'(A)')           '----'
      WRITE(100,'(A,A)')         'Coordinate System: ', TRIM( CoordinateSystem)
      WRITE(100,'(A,3I5.4)')     'nX:     ', nX
      WRITE(100,'(A,3I3.2)')     'bcX:    ', bcX
      WRITE(100,'(A,3ES12.3E3)') 'xL:     ', xL
      WRITE(100,'(A,3ES12.3E3)') 'xR:     ', xR
      WRITE(100,'(A,I2.2)')      'nNodes: ', nNodes
      WRITE(100,*)
      WRITE(100,'(A)')           'Time-Stepping'
      WRITE(100,'(A)')           '-------------'
      WRITE(100,'(A,ES10.3E3)')  't_end:         ', t_end
      WRITE(100,'(A,F4.2)')      'CFL:           ', CFL
      WRITE(100,'(A,I1.1)')      'nStagesSSPRK:  ', nStagesSSPRK
      WRITE(100,'(A,L)')         'UseFixed_dt:   ', UseFixed_dt
      WRITE(100,'(A,L)')         'UseSourceTerm: ', UseSourceTerm
      WRITE(100,*)
      WRITE(100,'(A)')           'Slope Limiter'
      WRITE(100,'(A)')           '------------------'
      WRITE(100,'(A,L)')         'UseSlopeLimiter:           ', UseSlopeLimiter
      WRITE(100,'(A,L)')         'UseTroubledCellIndicator:  ', &
                                   UseTroubledCellIndicator
      WRITE(100,'(A,L)')         'UseCharacteristicLimiting: ', &
                                   UseCharacteristicLimiting
      WRITE(100,'(A,ES10.3E3)')  'BetaTVD:                   ', BetaTVD
      WRITE(100,'(A,ES10.3E3)')  'BetaTVB:                   ', BetaTVB
      WRITE(100,'(A,ES10.3E3)')  'SlopeTolerance:            ', SlopeTolerance
      WRITE(100,'(A,F5.3)')      'LimiterThresholdParameter: ', &
                                   LimiterThresholdParameter
      WRITE(100,*)
      WRITE(100,'(A)')           'Positivity Limiter'
      WRITE(100,'(A)')           '------------------'
      WRITE(100,'(A,L)')         'UsePositivityLimiter: ', UsePositivityLimiter
      WRITE(100,'(A,ES11.3E3)')  'Min_1: ', Min_1
      WRITE(100,'(A,ES11.3E3)')  'Min_2: ', Min_2
    CLOSE(100)

  END SUBROUTINE WriteProgramHeader

END PROGRAM ApplicationDriver
