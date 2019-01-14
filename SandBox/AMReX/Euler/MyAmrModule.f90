MODULE MyAmrModule

  USE iso_c_binding
  USE amrex_amr_module
  USE amrex_fort_module, ONLY: &
    amrex_real

  USE MyAmrDataModule

  IMPLICIT NONE

  REAL(amrex_real)                    :: t_end, dt_wrt, Gamma_IDEAL
  INTEGER                             :: nNodes, nStages, CoordSys, nLevels
  INTEGER,          ALLOCATABLE       :: MaxGridSize(:), nX(:), swX(:)
  REAL(amrex_real), ALLOCATABLE       :: xL(:), xR(:), dt(:)
  CHARACTER(LEN=:), ALLOCATABLE       :: ProgramName
  INTEGER,          ALLOCATABLE, SAVE :: StepNo(:), nSubSteps(:)

CONTAINS

  SUBROUTINE MyAmrInit

    TYPE(amrex_parmparse) :: PP
    INTEGER               :: iLevel

    IF( .NOT. amrex_amrcore_initialized() ) &
      CALL amrex_amrcore_init()

!!$    CALL amrex_init_virtual_functions ( MyMakeNewLevelFromScratch,     &
!!$                                        my_make_new_level_from_coarse, &
!!$                                        my_remake_level,               &
!!$                                        MyClearLevel,                  &
!!$                                        my_error_estimate )
 
    CALL amrex_parmparse_build( PP )
      CALL PP % get( 'dt_wrt',      dt_wrt )
      CALL PP % get( 't_end',       t_end )
      CALL PP % get( 'nNodes',      nNodes )
      CALL PP % get( 'nStages',     nStages )
      CALL PP % get( 'ProgramName', ProgramName )
      CALL PP % get( 'Gamma',       Gamma_IDEAL )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters geometry.* ---
    CALL amrex_parmparse_build( PP, 'geometry' )
      CALL PP % get   ( 'CoordSys', CoordSys )
      CALL PP % getarr( 'prob_lo',  xL )
      CALL PP % getarr( 'prob_hi',  xR )
    CALL amrex_parmparse_destroy( PP )

    ! --- Parameters amr.*
    CALL amrex_parmparse_build( PP, 'amr' )
      CALL PP % getarr( 'n_cell',      nX )
      CALL PP % getarr( 'swX',         swX )
      CALL PP % getarr( 'MaxGridSize', MaxGridSize )
      CALL PP % get   ( 'max_level',   nLevels )
    CALL amrex_parmparse_destroy( PP )

!!$    if (.not. amrex_is_all_periodic()) then
!!$       lo_bc = amrex_bc_foextrap
!!$       hi_bc = amrex_bc_foextrap
!!$    end if

    ALLOCATE( StepNo(0:nLevels) )
    StepNo = 0

    ALLOCATE( nSubSteps(0:nLevels) )
    nSubSteps(0) = 1
    DO iLevel = 1, nLevels
      nSubSteps(iLevel) = amrex_ref_ratio(iLevel-1)
    END DO

    ALLOCATE( dt(0:nLevels) )
    dt = 1.0e-4_amrex_real

    CALL InitializeDataAMReX

  END SUBROUTINE MyAmrInit

  SUBROUTINE MyAmrFinalize

    WRITE(*,*)
    WRITE(*,'(A)') 'Calling MyAmrFinalize...'
    WRITE(*,*)

    CALL FinalizeDataAMReX
    DEALLOCATE( dt )
    DEALLOCATE( nSubSteps )
    DEALLOCATE( StepNo )

  END SUBROUTINE MyAmrFinalize


!!$  SUBROUTINE MyMakeNewLevelFromScratch( iLevel, Time, pBA, pDM ) BIND(c)
!!$
!!$    USE MF_GeometryModule,        ONLY: &
!!$      MF_ComputeGeometryX
!!$    USE MF_InitializationModule,  ONLY: &
!!$      MF_InitializeFields
!!$    USE MF_Euler_UtilitiesModule, ONLY: &
!!$      MF_ComputeFromConserved
!!$
!!$    INTEGER,          INTENT(in), VALUE :: iLevel
!!$    REAL(amrex_real), INTENT(in), VALUE :: Time
!!$    TYPE(c_ptr),      INTENT(in), VALUE :: pBA, pDM
!!$
!!$    TYPE(amrex_boxarray)  :: BA
!!$    TYPE(amrex_distromap) :: DM
!!$    TYPE(amrex_mfiter)    :: MFI_GF, MFI_CF, MFI_PF, MFI_AF
!!$    TYPE(amrex_box)       :: BX
!!$    REAL(amrex_real), CONTIGUOUS, POINTER :: &
!!$      MF_uGF(:,:,:,:), MF_uCF(:,:,:,:), MF_uPF(:,:,:,:), MF_uAF(:,:,:,:)
!!$
!!$    BA = pBA
!!$    DM = pDM
!!$
!!$    t_new(iLevel) = Time
!!$    t_old(iLevel) = Time - 1.0e200_amrex_real
!!$
!!$    CALL MyClearLevel( iLevel )
!!$
!!$    CALL amrex_multifab_build( MF_uGF_new(iLevel), BA, DM, nComp, nGhost )
!!$    CALL amrex_multifab_build( MF_uCF_new(iLevel), BA, DM, nComp, nGhost )
!!$    CALL amrex_multifab_build( MF_uPF_new(iLevel), BA, DM, nComp, nGhost )
!!$    CALL amrex_multifab_build( MF_uAF_new(iLevel), BA, DM, nComp, nGhost )
!!$
!!$    CALL amrex_multifab_build( MF_uGF_old(iLevel), BA, DM, nComp, nGhost )
!!$    CALL amrex_multifab_build( MF_uCF_old(iLevel), BA, DM, nComp, nGhost )
!!$    CALL amrex_multifab_build( MF_uPF_old(iLevel), BA, DM, nComp, nGhost )
!!$    CALL amrex_multifab_build( MF_uAF_old(iLevel), BA, DM, nComp, nGhost )
!!$
!!$    IF( iLevel .GT. 0 .AND. do_reflux ) &
!!$      CALL amrex_fluxregister_build &
!!$             ( flux_reg(iLevel), BA, DM, amrex_ref_ratio(iLevel-1), &
!!$                iLevel, nComp )
!!$
!!$    CALL amrex_mfiter_build( MFI_GF, MF_uGF_new(iLevel) )
!!$    CALL amrex_mfiter_build( MFI_CF, MF_uCF_new(iLevel) )
!!$    CALL amrex_mfiter_build( MFI_PF, MF_uPF_new(iLevel) )
!!$    CALL amrex_mfiter_build( MFI_AF, MF_uAF_new(iLevel) )
!!$
!!$    DO WHILE( MFI_GF % next() )
!!$
!!$      BX = MFI_GF % tilebox()
!!$      MF_uGF => MF_uGF_new( iLevel ) % DataPtr( MFI_GF )
!!$
!!$      BX = MFI_CF % tilebox()
!!$      MF_uCF => MF_uCF_new( iLevel ) % DataPtr( MFI_CF )
!!$
!!$      BX = MFI_PF % tilebox()
!!$      MF_uPF => MF_uPF_new( iLevel ) % DataPtr( MFI_PF )
!!$
!!$      BX = MFI_AF % tilebox()
!!$      MF_uAF => MF_uAF_new( iLevel ) % DataPtr( MFI_AF )
!!$
!!$      CALL MF_ComputeGeometryX( MF_uGF )
!!$      CALL MF_InitializeFields &
!!$             ( TRIM( ProgramName ), MF_uGF, MF_uCF )
!!$      CALL MF_ComputeFromConserved &
!!$             ( MF_uGF, MF_uCF, MF_uPF, MF_uAF )
!!$
!!$    END DO
!!$
!!$    CALL amrex_mfiter_destroy( MFI_GF )
!!$    CALL amrex_mfiter_destroy( MFI_CF )
!!$    CALL amrex_mfiter_destroy( MFI_PF )
!!$    CALL amrex_mfiter_destroy( MFI_AF )
!!$
!!$  END SUBROUTINE MyMakeNewLevelFromScratch
!!$
!!$
!!$  SUBROUTINE MyClearLevel( iLevel ) BIND(c)
!!$
!!$    INTEGER, INTENT(in), VALUE :: iLevel
!!$
!!$    CALL amrex_multifab_destroy( MF_uGF_new(iLevel) )
!!$    CALL amrex_multifab_destroy( MF_uCF_new(iLevel) )
!!$    CALL amrex_multifab_destroy( MF_uPF_new(iLevel) )
!!$    CALL amrex_multifab_destroy( MF_uAF_new(iLevel) )
!!$
!!$    CALL amrex_multifab_destroy( MF_uGF_old(iLevel) )
!!$    CALL amrex_multifab_destroy( MF_uCF_old(iLevel) )
!!$    CALL amrex_multifab_destroy( MF_uPF_old(iLevel) )
!!$    CALL amrex_multifab_destroy( MF_uAF_old(iLevel) )
!!$
!!$    CALL amrex_fluxregister_destroy( flux_reg(iLevel) )
!!$
!!$  END SUBROUTINE MyClearLevel

END MODULE MyAmrModule
