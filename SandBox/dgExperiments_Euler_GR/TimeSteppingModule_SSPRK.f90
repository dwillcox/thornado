MODULE TimeSteppingModule_SSPRK

  USE KindModule, ONLY: &
    DP, Zero, One
  USE ProgramHeaderModule, ONLY: &
    iX_B0, iX_B1, iX_E0, iX_E1, &
    nDOFX
  USE FluidFieldsModule, ONLY: &
    nCF
  USE Euler_GR_SlopeLimiterModule, ONLY: &
    Euler_GR_ApplySlopeLimiter
  USE Euler_GR_PositivityLimiterModule, ONLY: &
    Euler_GR_ApplyPositivityLimiter
  
  IMPLICIT NONE
  PRIVATE

  INTEGER :: nStages_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: c_SSPRK
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: w_SSPRK
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: a_SSPRK

  REAL(DP), DIMENSION(:,:,:,:,:),   ALLOCATABLE :: U_SSPRK
  REAL(DP), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: D_SSPRK

  PUBLIC :: InitializeFluid_SSPRK
  PUBLIC :: UpdateFluid_SSPRK
  PUBLIC :: FinalizeFluid_SSPRK

  INTERFACE
    SUBROUTINE FluidIncrement( iX_B0, iX_E0, iX_B1, iX_E1, G, U, dU )
      USE KindModule, ONLY: DP
      INTEGER, INTENT(in)     :: &
        iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
      REAL(DP), INTENT(in)    :: &
        G (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      REAL(DP), INTENT(inout) :: &
        U (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
      REAL(DP), INTENT(out)   :: &
        dU(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)
    END SUBROUTINE FluidIncrement
  END INTERFACE

CONTAINS


  SUBROUTINE InitializeFluid_SSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: i

    nStages_SSPRK = nStages

    CALL InitializeSSPRK( nStages )

    WRITE(*,*)
    WRITE(*,'(A5,A,I1)') '', 'SSP RK Scheme: ', nStages

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'Butcher Table:'
    WRITE(*,'(A5,A)') '', '--------------'
    DO i = 1, nStages
      WRITE(*,'(A5,4ES14.4E3)') '', c_SSPRK(i), a_SSPRK(i,1:nStages)
    END DO
    WRITE(*,'(A5,A14,3ES14.4E3)') '', '', w_SSPRK(1:nStages)

    ALLOCATE( U_SSPRK &
                (1:nDOFX, &
                 iX_B1(1):iX_E1(1), &
                 iX_B1(2):iX_E1(2), &
                 iX_B1(3):iX_E1(3), &
                 1:nCF) )

    ALLOCATE( D_SSPRK &
                (1:nDOFX, &
                 iX_B0(1):iX_E0(1), &
                 iX_B0(2):iX_E0(2), &
                 iX_B0(3):iX_E0(3), &
                 1:nCF,1:nStages) )

  END SUBROUTINE InitializeFluid_SSPRK


  SUBROUTINE FinalizeFluid_SSPRK

    DEALLOCATE( a_SSPRK, c_SSPRK, w_SSPRK )

    DEALLOCATE( U_SSPRK, D_SSPRK )

  END SUBROUTINE FinalizeFluid_SSPRK


  SUBROUTINE InitializeSSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    INTEGER :: iS

    CALL AllocateButcherTables_SSPRK( nStages )

    SELECT CASE ( nStages )
      CASE ( 1 )

        a_SSPRK(1,1) = 0.0_DP
        w_SSPRK(1)   = 1.0_DP

      CASE ( 2 )

        a_SSPRK(1,1:2) = [ 0.0_DP, 0.0_DP ]
        a_SSPRK(2,1:2) = [ 1.0_DP, 0.0_DP ]
        w_SSPRK(1:2)   = [ 0.5_DP, 0.5_DP ]

      CASE ( 3 )

        a_SSPRK(1,1:3) = [ 0.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(2,1:3) = [ 1.00_DP, 0.00_DP, 0.00_DP ]
        a_SSPRK(3,1:3) = [ 0.25_DP, 0.25_DP, 0.00_DP ]
        w_SSPRK(1:3)   = [ 1.0_DP / 6.0_DP, &
                           1.0_DP / 6.0_DP, &
                           2.0_DP / 3.0_DP ]

    END SELECT

    DO iS = 1, nStages
      c_SSPRK(iS) = SUM( a_SSPRK(iS,1:iS-1) )
    END DO

  END SUBROUTINE InitializeSSPRK


  SUBROUTINE AllocateButcherTables_SSPRK( nStages )

    INTEGER, INTENT(in) :: nStages

    ALLOCATE( a_SSPRK(nStages,nStages) )
    ALLOCATE( c_SSPRK(nStages) )
    ALLOCATE( w_SSPRK(nStages) )

    a_SSPRK = Zero
    c_SSPRK = Zero
    w_SSPRK = Zero

  END SUBROUTINE AllocateButcherTables_SSPRK


  SUBROUTINE UpdateFluid_SSPRK( t, dt, G, U, ComputeIncrement_Fluid )

    REAL(DP), INTENT(in) :: &
      t, dt
    REAL(DP), INTENT(in) :: &
      G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    PROCEDURE (FluidIncrement) :: &
      ComputeIncrement_Fluid
    LOGICAL :: DEBUG = .FALSE.

    INTEGER :: iS, jS

    U_SSPRK = Zero ! --- State
    D_SSPRK = Zero ! --- Increment

    DO iS = 1, nStages_SSPRK

      U_SSPRK = U

      DO jS = 1, iS - 1

        IF( a_SSPRK(iS,jS) .NE. Zero )THEN

          IF( DEBUG ) WRITE(*,'(A)') 'CALL AddIncrement_Fluid (1)'
          CALL AddIncrement_Fluid &
                 ( One, &
                   U_SSPRK(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                   dt * a_SSPRK(iS,jS), &
                   D_SSPRK(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,jS) )

        END IF

      END DO

      IF( ANY( a_SSPRK(:,iS) .NE. Zero ) &
          .OR. ( w_SSPRK(iS) .NE. Zero ) )THEN

        IF( DEBUG ) WRITE(*,'(A)') 'CALL Euler_GR_ApplySlopeLimiter'
        CALL Euler_GR_ApplySlopeLimiter &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                 G      (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                 U_SSPRK(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )

        IF( DEBUG ) WRITE(*,'(A)') 'CALL Euler_GR_ApplyPositivityLimiter'
        CALL Euler_GR_ApplyPositivityLimiter &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                 G      (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                 U_SSPRK(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )

        IF( DEBUG ) WRITE(*,'(A)') 'CALL ComputeIncrement_Fluid'
        CALL ComputeIncrement_Fluid &
               ( iX_B0, iX_E0, iX_B1, iX_E1, &
                 G      (1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                 U_SSPRK(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                 D_SSPRK(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,iS) )

      END IF

    END DO

    DO iS = 1, nStages_SSPRK

      IF( w_SSPRK(iS) .NE. Zero )THEN
          
        IF( DEBUG ) WRITE(*,'(A)') 'CALL AddIncrement_Fluid (2)'
        CALL AddIncrement_Fluid &
               ( One, &
                 U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
                 dt * w_SSPRK(iS), &
                 D_SSPRK(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:,iS) )

      END IF

    END DO

    IF( DEBUG ) WRITE(*,'(A)') 'CALL Euler_GR_ApplySlopeLimiter (2)'
    CALL Euler_GR_ApplySlopeLimiter &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
             U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )

    IF( DEBUG ) WRITE(*,'(A)') 'CALL Euler_GR_ApplyPositivityLimiter (2)'
    CALL Euler_GR_ApplyPositivityLimiter &
           ( iX_B0, iX_E0, iX_B1, iX_E1, &
             G(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:), &
             U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:) )


  END SUBROUTINE UpdateFluid_SSPRK
 

  SUBROUTINE AddIncrement_Fluid( alpha, U, beta, D )

    REAL(DP), INTENT(in)    :: &
      alpha, beta
    REAL(DP), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    REAL(DP), INTENT(in)    :: &
      D(1:,iX_B0(1):,iX_B0(2):,iX_B0(3):,1:)

    INTEGER :: iCF, iX1, iX2, iX3

    DO iCF = 1, nCF
      DO iX3 = iX_B0(3), iX_E0(3)
        DO iX2 = iX_B0(2), iX_E0(2)
          DO iX1 = iX_B0(1), iX_E0(1)

            U(:,iX1,iX2,iX3,iCF) &
              = alpha * U(:,iX1,iX2,iX3,iCF) &
                  + beta * D(:,iX1,iX2,iX3,iCF)

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE AddIncrement_Fluid


END MODULE TimeSteppingModule_SSPRK
