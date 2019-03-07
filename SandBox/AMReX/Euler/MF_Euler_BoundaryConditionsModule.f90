MODULE MF_Euler_BoundaryConditionsModule

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ConstructEdgeMap
  PUBLIC :: ApplyNonPeriodicBoundaryConditions_Fluid

  TYPE, PUBLIC :: EdgeMap
     LOGICAL :: is_lo_edge(3)
     LOGICAL :: is_hi_edge(3)
   CONTAINS
     PROCEDURE :: get_euler_bc => edgemap_get_euler_bc
  END type EdgeMap

CONTAINS

  SUBROUTINE edgemap_get_euler_bc(this, dim, Apply_Euler_BC)
    USE Euler_BoundaryConditionsModule, ONLY: &
         APPLY_EULER_BC_INNER, APPLY_EULER_BC_OUTER, APPLY_EULER_BC_BOTH, APPLY_EULER_BC_NONE

    CLASS(EdgeMap), INTENT(in) :: this
    INTEGER, INTENT(in)  :: dim
    INTEGER, INTENT(out) :: Apply_Euler_BC

    IF (this % is_lo_edge(dim) .AND. this % is_hi_edge(dim)) THEN
       Apply_Euler_BC = APPLY_EULER_BC_BOTH
    ELSE IF (this % is_lo_edge(dim)) THEN
       Apply_Euler_BC = APPLY_EULER_BC_INNER
    ELSE IF (this % is_hi_edge(dim)) THEN
       Apply_Euler_BC = APPLY_EULER_BC_OUTER
    ELSE
       Apply_Euler_BC = APPLY_EULER_BC_NONE
    END IF
  END SUBROUTINE edgemap_get_euler_bc


  SUBROUTINE ConstructEdgeMap( Geom, Box, Edge_Map )

    USE amrex_base_module, ONLY: &
         amrex_box,      &
         amrex_geometry
    USE amrex_fort_module, ONLY: &
         amrex_spacedim

    TYPE(amrex_geometry), INTENT(in)    :: Geom
    TYPE(amrex_box),      INTENT(in)    :: Box
    TYPE(EdgeMap),        INTENT(inout) :: Edge_Map

    INTEGER :: i

    DO i = 1, 3
       Edge_Map % is_lo_edge(i) = .FALSE.
       Edge_Map % is_hi_edge(i) = .FALSE.
       IF (i <= amrex_spacedim) THEN
          IF ( Box % lo(i) <= Geom % domain % lo(i) ) THEN
             Edge_Map % is_lo_edge(i) = .TRUE.
          END IF
          IF ( Box % hi(i) >= Geom % domain % hi(i) ) THEN
             Edge_Map % is_hi_edge(i) = .TRUE.
          END IF
       END IF
    END DO

  END SUBROUTINE ConstructEdgeMap


  SUBROUTINE ApplyNonperiodicBoundaryConditions_Fluid &
       ( iX_B0, iX_E0, iX_B1, iX_E1, U, Edge_Map )

    USE amrex_fort_module, ONLY: &
         amrex_real
    USE Euler_BoundaryConditionsModule, ONLY: &
         ApplyBC_Fluid_X1, ApplyBC_Fluid_X2, ApplyBC_Fluid_X3

    INTEGER,  INTENT(in)    :: &
      iX_B0(3), iX_E0(3), iX_B1(3), iX_E1(3)
    REAL(amrex_real), INTENT(inout) :: &
      U(1:,iX_B1(1):,iX_B1(2):,iX_B1(3):,1:)
    TYPE(EdgeMap), INTENT(in) :: Edge_Map

    INTEGER :: Apply_Euler_BC

    CALL Edge_Map % get_euler_bc(1, Apply_Euler_BC)
    CALL ApplyBC_Fluid_X1 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, Apply_Euler_BC )

    CALL Edge_Map % get_euler_bc(2, Apply_Euler_BC)
    CALL ApplyBC_Fluid_X2 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, Apply_Euler_BC )

    CALL Edge_Map % get_euler_bc(3, Apply_Euler_BC)
    CALL ApplyBC_Fluid_X3 &
           ( iX_B0, iX_E0, iX_B1, iX_E1, U, Apply_Euler_BC )

  END SUBROUTINE ApplyNonperiodicBoundaryConditions_Fluid

END MODULE MF_Euler_BoundaryConditionsModule
