MODULE ApplicationBoundaryConditionsModule

  USE KindModule, ONLY: &
    DP, Pi
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    xL, xR, &
    nX, nNodesX, &
    nE, nNodesE
  USE UtilitiesModule, ONLY: &
    NodeNumberX, &
    NodeNumber
  USE MeshModule, ONLY: &
    MeshX, &
    NodeCoordinate
  USE FluidFieldsModule, ONLY: &
    uCF, iCF_S1, nCF, &
    uPF, iPF_D, iPF_V1, iPF_V2, iPF_V3, iPF_E, iPF_Ne, nPF, &
    uAF, iAF_P, iAF_T, iAF_Ye, iAF_E, iAF_Gm, iAF_Cs, nAF
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, iCR_N, iCR_G1, iCR_G2, iCR_G3
  USE EquationOfStateModule, ONLY: &
    ComputeInternalEnergyDensityFromPressure, &
    ComputeAuxiliary_Fluid
  USE EulerEquationsUtilitiesModule, ONLY: &
    ComputeConserved

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ApplyApplicationBoundaryConditions_Fluid_X1
  PUBLIC :: ApplyApplicationBoundaryConditions_Radiation_X1

CONTAINS


  SUBROUTINE ApplyApplicationBoundaryConditions_Fluid_X1

    SELECT CASE ( TRIM( ProgramName ) )

      CASE ( 'HydrostaticPolytrope1D' )

        CALL ApplyBC_Fluid_X1_HydrostaticPolytrope1D

      CASE ( 'EvrardsCollapse1D' )

        CALL ApplyBC_Fluid_X1_EvrardsCollapse1D

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A38)') &
          '', 'Application Boundary Condition Missing'
        WRITE(*,*)
        WRITE(*,'(A7,A)') &
          '', TRIM( ProgramName )
        WRITE(*,*)
        STOP

    END SELECT

  END SUBROUTINE ApplyApplicationBoundaryConditions_Fluid_X1


  SUBROUTINE ApplyBC_Fluid_X1_HydrostaticPolytrope1D

    INTEGER  :: iX2, iX3
    INTEGER  :: iNodeX1, jNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iNodeX, jNodeX
    INTEGER  :: iCF, iPF, iAF
    REAL(DP) :: D_M, Kappa, Gamma

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)

        ! --- Inner Boundary (Reflecting) ---

        DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
            DO iNodeX1 = 1, nNodesX(1)

              jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

              iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
              jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

              ! -- Conserved --

              DO iCF = 1, nCF

                uCF(iNodeX,0,iX2,iX3,iCF) &
                  = uCF(jNodeX,1,iX2,iX3,iCF)

              END DO

              uCF(iNodeX,0,iX2,iX3,iCF_S1) &
                = - uCF(jNodeX,1,iX2,iX3,iCF_S1)

              ! -- Primitive --

              DO iPF = 1, nPF

                uPF(iNodeX,0,iX2,iX3,iPF) &
                  = uPF(jNodeX,1,iX2,iX3,iPF)

              END DO

              uPF(iNodeX,0,iX2,iX3,iPF_V1) &
                = - uPF(jNodeX,1,iX2,iX3,iPF_V1)

              ! -- Auxiliary --

              DO iAF = 1, nAF

                uAF(iNodeX,0,iX2,iX3,iAF) &
                  = uAF(jNodeX,1,iX2,iX3,iAF)

              END DO

            END DO
          END DO
        END DO

        ! --- Outer Boundary (Dirichlet) ---

        D_M   = 1.0d-6
        Kappa = 2.0_DP / Pi
        Gamma = 2.0_DP

        uPF(:,nX(1)+1,iX2,iX3,iPF_D)  = D_M
        uPF(:,nX(1)+1,iX2,iX3,iPF_V1) = 0.0_DP
        uPF(:,nX(1)+1,iX2,iX3,iPF_V2) = 0.0_DP
        uPF(:,nX(1)+1,iX2,iX3,iPF_V3) = 0.0_DP
        uAF(:,nX(1)+1,iX2,iX3,iAF_P)  = Kappa * D_M**Gamma

        CALL ComputeInternalEnergyDensityFromPressure &
               ( uPF(:,nX(1)+1,iX2,iX3,iPF_D),  uAF(:,nX(1)+1,iX2,iX3,iAF_P), &
                 uAF(:,nX(1)+1,iX2,iX3,iAF_Ye), uPF(:,nX(1)+1,iX2,iX3,iPF_E) )

        CALL ComputeAuxiliary_Fluid &
               ( uPF(:,nX(1)+1,iX2,iX3,iPF_D),  uPF(:,nX(1)+1,iX2,iX3,iPF_E),  &
                 uPF(:,nX(1)+1,iX2,iX3,iPF_Ne), uAF(:,nX(1)+1,iX2,iX3,iAF_P),  &
                 uAF(:,nX(1)+1,iX2,iX3,iAF_T),  uAF(:,nX(1)+1,iX2,iX3,iAF_Ye), &
                 uAF(:,nX(1)+1,iX2,iX3,iAF_E),  uAF(:,nX(1)+1,iX2,iX3,iAF_Gm), &
                 uAF(:,nX(1)+1,iX2,iX3,iAF_Cs) )

        CALL ComputeConserved &
               ( uPF(:,nX(1)+1,iX2,iX3,1:nPF), uCF(:,nX(1)+1,iX2,iX3,1:nCF) )

      END DO
    END DO

  END SUBROUTINE ApplyBC_Fluid_X1_HydrostaticPolytrope1D


  SUBROUTINE ApplyBC_Fluid_X1_EvrardsCollapse1D

    INTEGER  :: iX2, iX3
    INTEGER  :: iNodeX1, jNodeX1, iNodeX2, iNodeX3
    INTEGER  :: iNodeX, jNodeX
    INTEGER  :: iCF, iPF, iAF

    DO iX3 = 1, nX(3)
      DO iX2 = 1, nX(2)

        ! --- Inner Boundary (Reflecting) ---

        DO iNodeX3 = 1, nNodesX(3)
          DO iNodeX2 = 1, nNodesX(2)
            DO iNodeX1 = 1, nNodesX(1)

              jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1

              iNodeX = NodeNumberX( iNodeX1, iNodeX2, iNodeX3 )
              jNodeX = NodeNumberX( jNodeX1, iNodeX2, iNodeX3 )

              ! -- Conserved --

              DO iCF = 1, nCF

                uCF(iNodeX,0,iX2,iX3,iCF) &
                  = uCF(jNodeX,1,iX2,iX3,iCF)

              END DO

              uCF(iNodeX,0,iX2,iX3,iCF_S1) &
                = - uCF(jNodeX,1,iX2,iX3,iCF_S1)

              ! -- Primitive --

              DO iPF = 1, nPF

                uPF(iNodeX,0,iX2,iX3,iPF) &
                  = uPF(jNodeX,1,iX2,iX3,iPF)

              END DO

              uPF(iNodeX,0,iX2,iX3,iPF_V1) &
                = - uPF(jNodeX,1,iX2,iX3,iPF_V1)

              ! -- Auxiliary --

              DO iAF = 1, nAF

                uAF(iNodeX,0,iX2,iX3,iAF) &
                  = uAF(jNodeX,1,iX2,iX3,iAF)

              END DO

            END DO
          END DO
        END DO

        ! --- Outer Boundary (Dirichlet) ---

        uPF(:,nX(1)+1,iX2,iX3,iPF_D)  &
          = 1.0_DP / ( 4.0_DP * Pi )
        uPF(:,nX(1)+1,iX2,iX3,iPF_V1) = 0.0_DP
        uPF(:,nX(1)+1,iX2,iX3,iPF_V2) = 0.0_DP
        uPF(:,nX(1)+1,iX2,iX3,iPF_V3) = 0.0_DP
        uAF(:,nX(1)+1,iX2,iX3,iAF_P) &
          = 0.1_DP * uPF(:,nX(1)+1,iX2,iX3,iPF_D) / 3.0_DP

        CALL ComputeInternalEnergyDensityFromPressure &
               ( uPF(:,nX(1)+1,iX2,iX3,iPF_D),  uAF(:,nX(1)+1,iX2,iX3,iAF_P), &
                 uAF(:,nX(1)+1,iX2,iX3,iAF_Ye), uPF(:,nX(1)+1,iX2,iX3,iPF_E) )

        CALL ComputeAuxiliary_Fluid &
               ( uPF(:,nX(1)+1,iX2,iX3,iPF_D),  uPF(:,nX(1)+1,iX2,iX3,iPF_E),  &
                 uPF(:,nX(1)+1,iX2,iX3,iPF_Ne), uAF(:,nX(1)+1,iX2,iX3,iAF_P),  &
                 uAF(:,nX(1)+1,iX2,iX3,iAF_T),  uAF(:,nX(1)+1,iX2,iX3,iAF_Ye), &
                 uAF(:,nX(1)+1,iX2,iX3,iAF_E),  uAF(:,nX(1)+1,iX2,iX3,iAF_Gm), &
                 uAF(:,nX(1)+1,iX2,iX3,iAF_Cs) )

        CALL ComputeConserved &
               ( uPF(:,nX(1)+1,iX2,iX3,1:nPF), uCF(:,nX(1)+1,iX2,iX3,1:nCF) )

      END DO
    END DO

  END SUBROUTINE ApplyBC_Fluid_X1_EvrardsCollapse1D


  SUBROUTINE ApplyApplicationBoundaryConditions_Radiation_X1 &
               ( Time, LimiterBC_Option )

    REAL(DP), INTENT(in) :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: LimiterBC_Option

    SELECT CASE ( TRIM( ProgramName ) )

      CASE ( 'GaussianSphericalWave1D' )

        CALL ApplyBC_Radiation_X1_GaussianSphericalWave1D &
               ( Time, LimiterBC_Option )

      CASE ( 'HomogeneousSphere1D' )

        CALL ApplyBC_Radiation_X1_OutflowSphericalSymmetry

      CASE ( 'GaussianSphericalDiffusion1D' )

        CALL ApplyBC_Radiation_X1_OutflowSphericalSymmetry

      CASE ( 'CoolingProblem1D' )

        CALL ApplyBC_Radiation_X1_OutflowSphericalSymmetry

      CASE DEFAULT

        WRITE(*,*)
        WRITE(*,'(A5,A38)') &
          '', 'Application Boundary Condition Missing'
        WRITE(*,*)
        WRITE(*,'(A7,A)') &
          '', TRIM( ProgramName )
        WRITE(*,*)
        STOP

    END SELECT

  END SUBROUTINE ApplyApplicationBoundaryConditions_Radiation_X1


  SUBROUTINE ApplyBC_Radiation_X1_GaussianSphericalWave1D &
               ( Time, LimiterBC_Option )

    REAL(DP), INTENT(in) :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: LimiterBC_Option

    LOGICAL  :: LimiterBC
    INTEGER  :: iS, iX1, iX2, iX3, iE
    INTEGER  :: iNodeX1, iNodeX2, iNodeX3, iNodeE, iNode
    REAL(DP) :: X1

    LimiterBC = .FALSE.
    IF( PRESENT( LimiterBC_Option ) ) &
      LimiterBC = LimiterBC_Option

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)
                  DO iNodeE = 1, nNodesE

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    ! -- Inner Boundary:

                    IF( LimiterBC )THEN

                      X1 = NodeCoordinate( MeshX(1), 0, iNodeX1 )

                    ELSE

                      X1 = xL(1)

                    END IF

                    uCR(iNode,iE,0,iX2,iX3,iCR_N,iS) &
                      = EXP( - ( X1 - Time )**2 ) / X1**2

                    uCR(iNode,iE,0,iX2,iX3,iCR_G1,iS) &
                      = uCR(iNode,iE,0,iX2,iX3,iCR_N,iS)

                    uCR(iNode,iE,0,iX2,iX3,iCR_G2,iS) &
                      = 0.0_DP

                    uCR(iNode,iE,0,iX2,iX3,iCR_G3,iS) &
                      = 0.0_DP

                    ! -- Outer Boundary:

                    IF( LimiterBC )THEN

                      X1 = NodeCoordinate( MeshX(1), nX(1)+1, iNodeX1 )

                    ELSE

                      X1 = xR(1)

                    END IF

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_N,iS) &
                      = EXP( - ( X1 - Time )**2 ) / X1**2

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G1,iS) &
                      = uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_N,iS)

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G2,iS) &
                      = 0.0_DP

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G3,iS) &
                      = 0.0_DP

                  END DO
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ApplyBC_Radiation_X1_GaussianSphericalWave1D


  SUBROUTINE ApplyBC_Radiation_X1_OutflowSphericalSymmetry

    INTEGER :: iS, iX1, iX2, iX3, iE
    INTEGER :: iNodeX1, iNodeX2, iNodeX3, iNodeE, iNode
    INTEGER :: jNodeX1, jNode

    DO iS = 1, nSpecies

      DO iX3 = 1, nX(3)
        DO iX2 = 1, nX(2)
          DO iE = 1, nE

            DO iNodeX3 = 1, nNodesX(3)
              DO iNodeX2 = 1, nNodesX(2)
                DO iNodeX1 = 1, nNodesX(1)
                  DO iNodeE = 1, nNodesE

                    iNode = NodeNumber( iNodeE, iNodeX1, iNodeX2, iNodeX3 )

                    ! -- Inner Boundary: Reflecting

                    jNodeX1 = ( nNodesX(1) - iNodeX1 ) + 1
                    jNode   = NodeNumber( iNodeE, jNodeX1, iNodeX2, iNodeX3 )

                    uCR(iNode,iE,0,iX2,iX3,iCR_N,iS) &
                      = uCR(jNode,iE,1,iX2,iX3,iCR_N,iS)

                    uCR(iNode,iE,0,iX2,iX3,iCR_G1,iS) &
                      = - uCR(jNode,iE,1,iX2,iX3,iCR_G1,iS)

                    uCR(iNode,iE,0,iX2,iX3,iCR_G2,iS) &
                      = + uCR(jNode,iE,1,iX2,iX3,iCR_G2,iS)

                    uCR(iNode,iE,0,iX2,iX3,iCR_G3,iS) &
                      = + uCR(jNode,iE,1,iX2,iX3,iCR_G3,iS)

                    ! -- Outer Boundary: Homogeneous

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_N,iS) &
                      = uCR(iNode,iE,nX(1),iX2,iX3,iCR_N,iS)

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G1,iS) &
                      = uCR(iNode,iE,nX(1),iX2,iX3,iCR_G1,iS)

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G2,iS) &
                      = uCR(iNode,iE,nX(1),iX2,iX3,iCR_G2,iS)

                    uCR(iNode,iE,nX(1)+1,iX2,iX3,iCR_G3,iS) &
                      = + uCR(iNode,iE,nX(1),iX2,iX3,iCR_G3,iS)

                  END DO
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO

    END DO

  END SUBROUTINE ApplyBC_Radiation_X1_OutflowSphericalSymmetry


END MODULE ApplicationBoundaryConditionsModule
