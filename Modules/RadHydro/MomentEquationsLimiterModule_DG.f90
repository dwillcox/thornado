MODULE MomentEquationsLimiterModule_DG

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOF, nDimsE, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumber, &
    MinModB, &
    GetRoots_Quadratic
  USE PolynomialBasisModule_Lagrange, ONLY: &
    evalL, L_E, L_X1, L_X2, L_X3
  USE PolynomialBasisModule_Legendre, ONLY: &
    evalP, P_E, P_X1, P_X2, P_X3, &
    IndP_Q, MassP
  USE PolynomialBasisMappingModule, ONLY: &
    MapNodalToModal_Radiation, &
    MapModalToNodal_Radiation
  USE MeshModule, ONLY: &
    MeshE, &
    MeshX
  USE GeometryFieldsModule, ONLY: &
    Vol, VolJac
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    WeightsR, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3, nCR
  USE RadiationFieldsUtilitiesModule, ONLY: &
    CellAverage

  IMPLICIT NONE
  PRIVATE

  LOGICAL,  PARAMETER                   :: Debug = .FALSE.
  LOGICAL                               :: ApplySlopeLimiter
  LOGICAL                               :: ApplyPositivityLimiter
  INTEGER                               :: nPoints
  REAL(DP)                              :: BetaTVB
  REAL(DP)                              :: BetaTVD
  REAL(DP), PARAMETER                   :: Tol_TVD = 1.0d-2
  REAL(DP), PARAMETER                   :: Tol_N = 1.0d-100
  REAL(DP), PARAMETER                   :: Tol_G = 1.0d-12
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: Points_E
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: Points_X1
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: Points_X2
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: Points_X3
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: uCR_P
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: uCR_M
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Lagrange
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: Legendre

  PUBLIC :: InitializeLimiters_M1_DG
  PUBLIC :: ApplySlopeLimiter_M1_DG
  PUBLIC :: ApplyPositivityLimiter_M1_DG

CONTAINS


  SUBROUTINE InitializeLimiters_M1_DG &
               ( ApplySlopeLimiter_Option, BetaTVB_Option, BetaTVD_Option, &
                 ApplyPositivityLimiter_Option )

    LOGICAL,  INTENT(in), OPTIONAL :: ApplySlopeLimiter_Option
    REAL(DP), INTENT(in), OPTIONAL :: BetaTVB_Option
    REAL(DP), INTENT(in), OPTIONAL :: BetaTVD_Option
    LOGICAL,  INTENT(in), OPTIONAL :: ApplyPositivityLimiter_Option

    INTEGER :: iPol, iNodeE, iNodeX1, iNodeX2, iNodeX3, iPoint, iNode
    REAL(DP), DIMENSION(:), ALLOCATABLE :: NodesX1

    ! --- Limiter Parameters ---

    ApplySlopeLimiter = .TRUE.
    IF( PRESENT( ApplySlopeLimiter_Option ) ) &
      ApplySlopeLimiter = ApplySlopeLimiter_Option

    BetaTVB = 1.0d0
    IF( PRESENT( BetaTVB_Option ) ) &
      BetaTVB = BetaTVB_Option

    BetaTVD = 2.0d0
    IF( PRESENT( BetaTVD_Option ) ) &
      BetaTVD = BetaTVD_Option

    ApplyPositivityLimiter = .TRUE.
    IF( PRESENT( ApplyPositivityLimiter_Option ) ) &
      ApplyPositivityLimiter = ApplyPositivityLimiter_Option

    WRITE(*,*)
    WRITE(*,'(A5,A)') '', 'InitializeLimiters_M1_DG'
    WRITE(*,*)
    WRITE(*,'(A7,A20,L1)') &
      '', 'ApplySlopeLimiter = ', ApplySlopeLimiter
    WRITE(*,'(A7,A10,ES8.2E2)') '', 'BetaTVB = ', BetaTVB
    WRITE(*,'(A7,A10,ES8.2E2)') '', 'BetaTVD = ', BetaTVD
    WRITE(*,'(A7,A25,L1)') &
      '', 'ApplyPositivityLimiter = ', ApplyPositivityLimiter
    WRITE(*,*)

    ! --- Legendre Polynomials in Gaussian Quadrature Points ---

    ALLOCATE( Legendre(nDOF,nDOF) )

    DO iPol = 1, nDOF

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX2 = 1, nNodesX(2)
          DO iNodeX1 = 1, nNodesX(1)
            DO iNodeE = 1, nNodesE

              iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

              Legendre(iNode,iPol) &
                = P_E   (IndP_Q(0,iPol)) % P( MeshE    % Nodes(iNodeE)  ) &
                  * P_X1(IndP_Q(1,iPol)) % P( MeshX(1) % Nodes(iNodeX1) ) &
                  * P_X2(IndP_Q(2,iPol)) % P( MeshX(2) % Nodes(iNodeX2) ) &
                  * P_X3(IndP_Q(3,iPol)) % P( MeshX(3) % Nodes(iNodeX3) )

            END DO
          END DO
        END DO
      END DO

    END DO

    ! --- ---

    ALLOCATE( NodesX1(nNodesX(1)+2) )

    NodesX1 = [ - 0.5_DP, MeshX(1) % Nodes, + 0.5_DP ]

    nPoints = nNodesE * ( nNodesX(1) + 2 ) * nNodesX(2) * nNodesX(3)

    IF( nNodesX(2) > 1 ) &
      nPoints = nPoints + 2 * nNodesE * nNodesX(1) * nNodesX(3)

    IF( nNodesX(3) > 1 ) &
      nPoints = nPoints + 2 * nNodesE * nNodesX(2) * nNodesX(3)

    ALLOCATE( uCR_P(nPoints,nCR) )

    ALLOCATE( uCR_M(nDOF,nCR) )

    ! --- Coordinates of Points Where Positivity is Required:

    ALLOCATE &
      ( Points_E (nPoints), Points_X1(nPoints), &
        Points_X2(nPoints), Points_X3(nPoints) )

    iPoint = 0

    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, SIZE( NodesX1 )
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = NodesX1(iNodeX1)
            Points_X2(iPoint) = MeshX(2) % Nodes(iNodeX2)
            Points_X3(iPoint) = MeshX(3) % Nodes(iNodeX3)

          END DO
        END DO
      END DO
    END DO

    IF( nNodesX(2) > 1 )THEN

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X2(iPoint) = - 0.5_DP
            Points_X3(iPoint) = MeshX(3) % Nodes(iNodeX3)

          END DO
        END DO
      END DO

      DO iNodeX3 = 1, nNodesX(3)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X2(iPoint) = + 0.5_DP
            Points_X3(iPoint) = MeshX(3) % Nodes(iNodeX3)

          END DO
        END DO
      END DO

    END IF

    IF( nNodesX(3) > 1 )THEN

      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X2(iPoint) = MeshX(2) % Nodes(iNodeX2)
            Points_X3(iPoint) = - 0.5_DP

          END DO
        END DO
      END DO

      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iPoint = iPoint + 1

            Points_E (iPoint) = MeshE    % Nodes(iNodeE)
            Points_X1(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X2(iPoint) = MeshX(1) % Nodes(iNodeX1)
            Points_X3(iPoint) = + 0.5_DP

          END DO
        END DO
      END DO

    END IF

    DEALLOCATE( NodesX1 )

    ALLOCATE( Lagrange(nDOF,nPoints) )

    DO iNodeX3 = 1, nNodesX(3)
      DO iNodeX2 = 1, nNodesX(2)
        DO iNodeX1 = 1, nNodesX(1)
          DO iNodeE = 1, nNodesE

            iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

            DO iPoint = 1, nPoints

              Lagrange(iNode,iPoint) &
                = L_E(iNodeE) % P( Points_E(iPoint) ) &
                    * L_X1(iNodeX1) % P( Points_X1(iPoint) ) &
                        * L_X2(iNodeX2) % P( Points_X2(iPoint) ) &
                            * L_X3(iNodeX3) % P( Points_X3(iPoint) )

            END DO

          END DO
        END DO
      END DO
    END DO

  END SUBROUTINE InitializeLimiters_M1_DG


  SUBROUTINE ApplySlopeLimiter_M1_DG

    LOGICAL :: LimitPolynomial
    INTEGER :: iS, iE, iX1, iX2, iX3, iCR, iOS
    REAL(DP), DIMENSION(nCR) :: uCR_A_P_X1, uCR_A_N_X1
    REAL(DP), DIMENSION(nDOF,nCR) :: uCR_M_P_X1, uCR_M_N_X1
    REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE :: uCR_A, uCR_X1, uCR_X1_T

    IF( nDOF == 1 ) RETURN

    IF( .NOT. ApplySlopeLimiter ) RETURN

    IF( Debug )THEN
      WRITE(*,*)
      WRITE(*,'(A6,A)') '', 'ApplySlopeLimiter_M1_DG'
      WRITE(*,*)
    END IF

    ALLOCATE &
      ( uCR_A(1:nCR,1:nE,1:nX(1),1:nX(2),1:nX(3)), &
        uCR_X1  (1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR), &
        uCR_X1_T(1:nE,1:nX(1),1:nX(2),1:nX(3),1:nCR) )

    ASSOCIATE( dX1 => MeshX(1) % Width(1:nX(1)) )

    iOS = 1 ! --- Offset to Access Position Space Slopes
    IF( nDimsE == 0 ) iOS = 0

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              ! --- Limiting Using Modal Representation ---

              DO iCR = 1, nCR

                ! --- Map To Modal Representation ---

                CALL MapNodalToModal_Radiation &
                       ( VolJac(:,iE,iX1-1,iX2,iX3) &
                           * uCR(:,iE,iX1-1,iX2,iX3,iCR,iS), uCR_M_P_X1(:,iCR) )
                CALL MapNodalToModal_Radiation &
                       ( VolJac(:,iE,iX1,  iX2,iX3) &
                           * uCR(:,iE,iX1,  iX2,iX3,iCR,iS), uCR_M     (:,iCR) )
                CALL MapNodalToModal_Radiation &
                       ( VolJac(:,iE,iX1+1,iX2,iX3) &
                           * uCR(:,iE,iX1+1,iX2,iX3,iCR,iS), uCR_M_N_X1(:,iCR) )

              END DO

              uCR_A(1:nCR,iE,iX1,iX2,iX3) = uCR_M     (1,1:nCR)
              uCR_A_P_X1(1:nCR)           = uCR_M_P_X1(1,1:nCR)
              uCR_A_N_X1(1:nCR)           = uCR_M_N_X1(1,1:nCR)

              ! --- Slope From Modal Representation ---

              uCR_X1(iE,iX1,iX2,iX3,1:nCR) &
                = uCR_M(iOS+2,1:nCR) ! X1-Dimension

              ! --- Compute Limited Slopes ---

              uCR_X1_T(iE,iX1,iX2,iX3,1:nCR) &
                = MinModB &
                    ( uCR_X1(iE,iX1,iX2,iX3,1:nCR), &
                      BetaTVD * ( uCR_A(1:nCR,iE,iX1,iX2,iX3) - uCR_A_P_X1 ), &
                      BetaTVD * ( uCR_A_N_X1 - uCR_A(1:nCR,iE,iX1,iX2,iX3) ), &
                      dX1(iX1), BetaTVB )

            END DO
          END DO
        END DO
      END DO

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              ! --- Component-wise limiting ---

              DO iCR = 1, nCR

                LimitPolynomial = .FALSE.
                LimitPolynomial &
                  = ( ABS( uCR_X1(iE,iX1,iX2,iX3,iCR) &
                           - uCR_X1_T(iE,iX1,iX2,iX3,iCR) ) &
                        / Vol(iE,iX1,iX2,iX3) > Tol_TVD )

                IF( LimitPolynomial )THEN

                  IF( Debug )THEN

                    WRITE(*,*)
                    WRITE(*,'(A4,A,1I2.2)') &
                      '', 'Limiting Radiation Field: ', iCR
                    WRITE(*,'(A6,A,4I5.4)') '', &
                      'iE, iX1, iX2, iX3 = ', iE, iX1, iX2, iX3
                    WRITE(*,*)

                    print*, " ux       = ", uCR_X1  (iE,iX1,iX2,iX3,iCR)
                    print*, " ux_tilde = ", uCR_X1_T(iE,iX1,iX2,iX3,iCR)

                  END IF

                  CALL MapNodalToModal_Radiation &
                         ( VolJac(:,iE,iX1,iX2,iX3) &
                             * uCR(:,iE,iX1,iX2,iX3,iCR,iS), uCR_M(:,iCR) )

                  ! --- Cell-Integrated Moments ---

                  uCR_M(:,iCR) = 0.0_DP
                  uCR_M(1,iCR) &     ! -- Cell-Average
                    = uCR_A(iCR,iE,iX1,iX2,iX3)
                  uCR_M(iOS+2,iCR) & ! -- Slope X1-Direction
                    = uCR_X1_T(iE,iX1,iX2,iX3,iCR)

                  ! --- Back to Nodal Representation ---

                  CALL MapModalToNodal_Radiation &
                         ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), uCR_M(:,iCR) )

                  uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                    = uCR(:,iE,iX1,iX2,iX3,iCR,iS) &
                        / VolJac(:,iE,iX1,iX2,iX3)

                END IF

              END DO

            END DO
          END DO
        END DO
      END DO

    END DO

    END ASSOCIATE ! dX1, etc.

    DEALLOCATE( uCR_A, uCR_X1, uCR_X1_T )

  END SUBROUTINE ApplySlopeLimiter_M1_DG


  SUBROUTINE ApplyPositivityLimiter_M1_DG

    INTEGER  :: iE, iX1, iX2, iX3, iS
    INTEGER  :: iCR, iPoint
    REAL(DP) :: Theta_1, Theta_2
    REAL(DP) :: dN, dG1, dG2, dG3
    REAL(DP) :: a, b, c, r1, r2
    REAL(DP), DIMENSION(nPoints) :: absG
    REAL(DP), DIMENSION(1:nCR)   :: uCR_K

    IF( nDOF == 1 ) RETURN

    IF( .NOT. ApplyPositivityLimiter ) RETURN

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iX1 = 1, nX(1)
            DO iE = 1, nE

              DO iCR = 1, nCR
                DO iPoint = 1, nPoints

                  uCR_P(iPoint,iCR) &
                    = DOT_PRODUCT &
                        ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), Lagrange(:,iPoint) )

                END DO
              END DO

              absG &
                = SQRT( uCR_P(:,iCR_G1)**2 + uCR_P(:,iCR_G2)**2 &
                        + uCR_P(:,iCR_G3)**2 )

              IF( ANY( uCR_P(:,iCR_N) < Tol_N ) .OR. &
                  ANY( uCR_P(:,iCR_N) - absG(:) < Tol_G ) )THEN

                ! --- Limiting Using Modal Representation ---

                DO iCR = 1, nCR

                  CALL MapNodalToModal_Radiation &
                         ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), uCR_M(:,iCR) )

                  uCR_K(iCR) &
                    = CellAverage &
                        ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), &
                          VolJac(:,iE,iX1,iX2,iX3) )

                END DO

                ! --- Ensure Positive Number Density ---

                Theta_1 = 1.0_DP
                DO iPoint = 1, nPoints
                  IF( uCR_P(iPoint,iCR_N) < Tol_N ) &
                    Theta_1 &
                      = MIN( Theta_1, &
                             ( Tol_N - uCR_K(iCR_N) ) &
                               / ( uCR_P(iPoint,iCR_N) - uCR_K(iCR_N) ) )
                END DO

                IF( Theta_1 < 1.0_DP )THEN

                  Theta_1 = ( 1.0_DP - SQRT( EPSILON( 1.0_DP ) ) ) * Theta_1

                  uCR_M(1,iCR_N) &
                    = ( 1.0_DP - Theta_1 ) * uCR_K(iCR_N) &
                        + Theta_1 * uCR_M(1,iCR_N)

                  uCR_M(2:nDOF,iCR_N) &
                    = Theta_1 * uCR_M(2:nDOF,iCR_N)

                  DO iPoint = 1, nPoints

                    uCR_P(iPoint,iCR_N) &
                      = evalP( uCR_M(:,iCR_N), &
                               Points_E (iPoint), Points_X1(iPoint), &
                               Points_X2(iPoint), Points_X3(iPoint) )

                  END DO

                END IF

                ! --- Ensure Limited Number Flux ---

                absG &
                  = SQRT( uCR_P(:,iCR_G1)**2 + uCR_P(:,iCR_G2)**2 &
                          + uCR_P(:,iCR_G3)**2 )

                Theta_2 = 1.0_DP
                DO iPoint = 1, nPoints

                  IF( uCR_P(iPoint,iCR_N) - absG(iPoint) < Tol_G )THEN

                    dN  = uCR_K(iCR_N)  - uCR_P(iPoint,iCR_N)
                    dG1 = uCR_K(iCR_G1) - uCR_P(iPoint,iCR_G1)
                    dG2 = uCR_K(iCR_G2) - uCR_P(iPoint,iCR_G2)
                    dG3 = uCR_K(iCR_G3) - uCR_P(iPoint,iCR_G3)

                    a = dN**2 - ( dG1**2 + dG2**2 + dG3**2 )

                    b = - 2.0_DP * ( uCR_K(iCR_N) * dN &
                                       * ( 1.0_DP - Tol_G / uCR_K(iCR_N) ) &
                                     - ( uCR_K(iCR_G1) * dG1 &
                                         + uCR_K(iCR_G2) * dG2 &
                                         + uCR_K(iCR_G3) * dG3 ) )

                    c = uCR_K(iCR_N)**2 &
                          * ( 1.0_DP - 2.0_DP * Tol_G / uCR_K(iCR_N) ) &
                        + Tol_G**2 &
                        - ( uCR_K(iCR_G1)**2 + uCR_K(iCR_G2)**2 &
                            + uCR_K(iCR_G3)**2 )

                    CALL GetRoots_Quadratic &
                           ( a, b, c, r1, r2, Tol_Option = 0.0d-14 )

                    IF( r1 < 0.0_DP .AND. r2 < 0.0_DP )THEN
                      Theta_2 = 0.0_DP
                    ELSE
                      IF( r1 > 0.0_DP ) Theta_2 = MIN( r1, Theta_2 )
                      IF( r2 > 0.0_DP ) Theta_2 = MIN( r2, Theta_2 )
                    END IF

                    Theta_2 = 0.0_DP ! --- DEBUG ---

                  END IF

                END DO

                IF( Theta_2 < 1.0_DP )THEN

                  DO iCR = 1, nCR

                    uCR_M(1,iCR) &
                      = ( 1.0_DP - Theta_2 ) * uCR_K(iCR) &
                          + Theta_2 * uCR_M(1,iCR)

                    uCR_M(2:nDOF,iCR) &
                      = Theta_2 * uCR_M(2:nDOF,iCR)

                  END DO

                END IF

                ! --- Back to Nodal Representation ---

                DO iCR = 1, nCR

                  CALL MapModalToNodal_Radiation &
                         ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), uCR_M(:,iCR) )

                END DO

                IF( Debug )THEN

                  DO iCR = 1, nCR
                    DO iPoint = 1, nPoints

                      uCR_P(iPoint,iCR) &
                        = evalL( uCR(:,iE,iX1,iX2,iX3,iCR,iS), &
                                 Points_E (iPoint), Points_X1(iPoint), &
                                 Points_X2(iPoint), Points_X3(iPoint) )

                    END DO

                    uCR_K(iCR) &
                      = CellAverage &
                          ( uCR(:,iE,iX1,iX2,iX3,iCR,iS), &
                            VolJac(:,iE,iX1,iX2,iX3) )

                  END DO

                  absG &
                    = SQRT( uCR_P(:,iCR_G1)**2 + uCR_P(:,iCR_G2)**2 &
                            + uCR_P(:,iCR_G3)**2 )

                  IF( ANY( uCR_P(:,iCR_N) < 0.0_DP ) .OR. &
                      ANY( uCR_P(:,iCR_N) - absG(:) < 0.0_DP ) )THEN

                    PRINT*
                    PRINT*, "ApplyPositivityLimiter_M1_DG"
                    PRINT*
                    PRINT*, "  Problem with Positivity Limiter!"
                    PRINT*
                    PRINT*, "    iE, iX1, iX2, iX3 = ", iE, iX1, iX2, iX3
                    PRINT*, "    Theta_1 = ", Theta_1
                    PRINT*, "    Theta_2 = ", Theta_2
                    PRINT*
                    PRINT*, "  Conserved Radiation Fields (Nodal):"
                    PRINT*, "  N     = ", uCR_P(:,iCR_N)
                    PRINT*, "  G1    = ", uCR_P(:,iCR_G1)
                    PRINT*, "  G2    = ", uCR_P(:,iCR_G2)
                    PRINT*, "  G3    = ", uCR_P(:,iCR_G3)
                    PRINT*, "  N-|G| = ", &
                      uCR_P(:,iCR_N) &
                      - SQRT( uCR_P(:,iCR_G1)**2 &
                              + uCR_P(:,iCR_G2)**2 &
                              + uCR_P(:,iCR_G3)**2 )
                    PRINT*
                    PRINT*, "  Cell-Averages:"
                    PRINT*, "  N_K       = ", uCR_K(iCR_N)
                    PRINT*, "  G1_K      = ", uCR_K(iCR_G1)
                    PRINT*, "  G2_K      = ", uCR_K(iCR_G2)
                    PRINT*, "  G3_K      = ", uCR_K(iCR_G3)
                    PRINT*, "  N_K-|G_K| = ", &
                      uCR_K(iCR_N) &
                      - SQRT( uCR_K(iCR_G1)**2 &
                              + uCR_K(iCR_G2)**2 &
                              + uCR_K(iCR_G3)**2 )
                    PRINT*

                    STOP

                  END IF

                END IF ! Debug

              END IF

            END DO
          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ApplyPositivityLimiter_M1_DG


END MODULE MomentEquationsLimiterModule_DG
