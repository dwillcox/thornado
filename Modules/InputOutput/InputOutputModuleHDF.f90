MODULE InputOutputModuleHDF

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    ProgramName, &
    nE, nNodesE, &
    nX, nNodesX, &
    nDOF
  USE ReferenceElementModule_Beta, ONLY: &
    NodeNumberTable
  USE MeshModule, ONLY: &
    MeshE, MeshX
  USE InputOutputUtilitiesModule, ONLY: &
    NodeCoordinates, &
    Field4D
  USE GeometryFieldsModule, ONLY: &
    uGF, nGF, namesGF
  USE FluidFieldsModule, ONLY: &
    uCF, nCF, namesCF, &
    uPF, nPF, namesPF, &
    uAF, nAF, namesAF
  USE RadiationFieldsModule, ONLY: &
    nSpecies, &
    uCR, nCR, namesCR, &
    uPR, nPR, namesPR, &
    uAR, nAR, namesAR

  USE HDF5

  IMPLICIT NONE
  PRIVATE

  CHARACTER(9),  PARAMETER :: &
    OutputDirectory = '../Output'
  CHARACTER(14), PARAMETER :: &
    GeometrySuffix  = 'GeometryFields'
  CHARACTER(11), PARAMETER :: &
    FluidSuffix     = 'FluidFields'
  CHARACTER(15), PARAMETER :: &
    RadiationSuffix = 'RadiationFields'
  INTEGER :: FileNumber = 0

  INTEGER :: HDFERR

  PUBLIC :: WriteFieldsHDF

CONTAINS


  SUBROUTINE WriteFieldsHDF &
              ( Time, WriteGF_Option, WriteFF_Option, WriteRF_Option )

    REAL(DP), INTENT(in) :: Time
    LOGICAL,  INTENT(in), OPTIONAL :: WriteGF_Option
    LOGICAL,  INTENT(in), OPTIONAL :: WriteFF_Option
    LOGICAL,  INTENT(in), OPTIONAL :: WriteRF_Option

    LOGICAL :: WriteGF
    LOGICAL :: WriteFF
    LOGICAL :: WriteRF

    WriteGF = .FALSE.
    IF( PRESENT( WriteGF_Option ) ) &
      WriteGF = WriteGF_Option

    WriteFF = .FALSE.
    IF( PRESENT( WriteFF_Option ) ) &
      WriteFF = WriteFF_Option

    WriteRF = .FALSE.
    IF( PRESENT( WriteRF_Option ) ) &
      WriteRF = WriteRF_Option

    IF( WriteGF )THEN

      CALL WriteGeometryFieldsHDF( Time )

    END IF

    IF( WriteFF )THEN

      CALL WriteFluidFieldsHDF( Time )

    END IF

    IF( WriteRF )THEN

      CALL WriteRadiationFieldsHDF( Time )

    END IF

    FileNumber = FileNumber + 1

  END SUBROUTINE WriteFieldsHDF


  SUBROUTINE WriteGeometryFieldsHDF( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: GroupName
    CHARACTER(256) :: DatasetName
    INTEGER        :: iGF
    INTEGER(HID_T) :: FILE_ID
    REAL(DP)       :: Dummy3D(2,2,2) = 0.0_DP

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        GeometrySuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    ! --- Write Time ---

    DatasetName = '/Time'

    CALL WriteDataset1DHDF &
           ( [ Time ], DatasetName, FILE_ID )

    ! --- Write Spatial Grid ---

    GroupName = 'Spatial Grid'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ) , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X1'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(1),nX(1),nNodesX(1)), &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X2'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(2),nX(2),nNodesX(2)), &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X3'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(3),nX(3),nNodesX(3)), &
             DatasetName, FILE_ID )

    ! --- Write Geometry Variables ---

    GroupName = 'Geometry Fields'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

    DO iGF = 1, nGF

      DatasetName = TRIM( GroupName ) // '/' // TRIM( namesGF(iGF) )

      CALL WriteDataset3DHDF &
             ( Dummy3D, DatasetName, FILE_ID )

    END DO

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

  END SUBROUTINE WriteGeometryFieldsHDF


  SUBROUTINE WriteFluidFieldsHDF( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: GroupName
    CHARACTER(256) :: GroupName2
    CHARACTER(256) :: DatasetName
    INTEGER        :: iFF
    INTEGER(HID_T) :: FILE_ID
    REAL(DP)       :: Dummy3D(2,2,2) = 0.0_DP

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        FluidSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    ! --- Write Time ---

    DatasetName = '/Time'

    CALL WriteDataset1DHDF &
           ( [ Time ], DatasetName, FILE_ID )

    ! --- Write Spatial Grid ---

    GroupName = 'Spatial Grid'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ) , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X1'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(1),nX(1),nNodesX(1)), &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X2'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(2),nX(2),nNodesX(2)), &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X3'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(3),nX(3),nNodesX(3)), &
             DatasetName, FILE_ID )

    ! --- Write Fluid Variables ---

    GroupName = 'Fluid Fields'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

    ! --- Conserved ---

    GroupName2 = TRIM( GroupName ) // '/Conserved'

    CALL CreateGroupHDF( FileName, TRIM( GroupName2 ), FILE_ID )

    DO iFF = 1, nCF

      DatasetName = TRIM( GroupName2 ) // '/' // TRIM( namesCF(iFF) )

      CALL WriteDataset3DHDF &
               ( Dummy3D, DatasetName, FILE_ID )

    END DO

    ! --- Primitive ---

    GroupName2 = TRIM( GroupName ) // '/Primitive'

    CALL CreateGroupHDF( FileName, TRIM( GroupName2 ), FILE_ID )

    DO iFF = 1, nPF

      DatasetName = TRIM( GroupName2 ) // '/' // TRIM( namesPF(iFF) )

      CALL WriteDataset3DHDF &
               ( Dummy3D, DatasetName, FILE_ID )

    END DO

    ! --- Auxiliary ---

    GroupName2 = TRIM( GroupName ) // '/Auxiliary'

    CALL CreateGroupHDF( FileName, TRIM( GroupName2 ), FILE_ID )

    DO iFF = 1, nAF

      DatasetName = TRIM( GroupName2 ) // '/' // TRIM( namesAF(iFF) )

      CALL WriteDataset3DHDF &
               ( Dummy3D, DatasetName, FILE_ID )

    END DO

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

  END SUBROUTINE WriteFluidFieldsHDF


  SUBROUTINE WriteRadiationFieldsHDF( Time )

    REAL(DP), INTENT(in) :: Time

    CHARACTER(2)   :: String2
    CHARACTER(6)   :: FileNumberString
    CHARACTER(256) :: FileName
    CHARACTER(256) :: GroupName
    CHARACTER(256) :: GroupNameSpecies
    CHARACTER(256) :: DatasetName
    INTEGER        :: iS, iRF
    INTEGER(HID_T) :: FILE_ID

    WRITE( FileNumberString, FMT='(i6.6)') FileNumber

    FileName &
      = OutputDirectory // '/' // &
        TRIM( ProgramName ) // '_' // &
        RadiationSuffix // '_' // &
        FileNumberString // '.h5'

    CALL H5OPEN_F( HDFERR )

    CALL H5FCREATE_F( TRIM( FileName ), H5F_ACC_TRUNC_F, FILE_ID, HDFERR )

    ! --- Write Time ---

    DatasetName = '/Time'

    CALL WriteDataset1DHDF &
           ( [ Time ], DatasetName, FILE_ID )

    ! --- Write Spatial Grid ---

    GroupName = 'Spatial Grid'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ) , FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X1'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(1),nX(1),nNodesX(1)), &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X2'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(2),nX(2),nNodesX(2)), &
             DatasetName, FILE_ID )

    DatasetName = TRIM( GroupName ) // '/X3'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshX(3),nX(3),nNodesX(3)), &
             DatasetName, FILE_ID )

    ! --- Write Energy Grid ---

    GroupName = 'Energy Grid'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

    DatasetName = TRIM( GroupName ) // '/E'

    CALL WriteDataset1DHDF &
           ( NodeCoordinates(MeshE,nE,nNodesE), DatasetName, FILE_ID )

    ! --- Write Radiation Variables ---

    GroupName = 'Radiation Fields'

    CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

    DO iS = 1, nSpecies

      WRITE( String2, FMT='(i2.2)') iS

      GroupNameSpecies = 'Radiation Fields/' // 'Species_' // String2

      CALL CreateGroupHDF( FileName, TRIM( GroupNameSpecies ), FILE_ID )

      ! --- Conserved ---

      GroupName = TRIM( GroupNameSpecies ) // '/Conserved'

      CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

      DO iRF = 1, nCR

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesCR(iRF) )

        CALL WriteDataset4DHDF &
               ( Field4D &
                   ( uCR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iRF,iS), &
                     [ nE, nX(1), nX(2), nX(3) ],                     &
                     [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                     nDOF, NodeNumberTable ), DatasetName, FILE_ID )

      END DO

      ! --- Primitive ---

      GroupName = TRIM( GroupNameSpecies  ) // '/Primitive'

      CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

      DO iRF = 1, nPR

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesPR(iRF) )

        CALL WriteDataset4DHDF &
               ( Field4D &
                   ( uPR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iRF,iS), &
                     [ nE, nX(1), nX(2), nX(3) ],                     &
                     [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                     nDOF, NodeNumberTable ), DatasetName, FILE_ID )

      END DO

      ! --- Auxiliary ---

      GroupName = TRIM( GroupNameSpecies  ) // '/Auxiliary'

      CALL CreateGroupHDF( FileName, TRIM( GroupName ), FILE_ID )

      DO iRF = 1, nAR

        DatasetName = TRIM( GroupName ) // '/' // TRIM( namesAR(iRF) )

        CALL WriteDataset4DHDF &
               ( Field4D &
                   ( uAR(1:nDOF,1:nE,1:nX(1),1:nX(2),1:nX(3),iRF,iS), &
                     [ nE, nX(1), nX(2), nX(3) ],                     &
                     [ nNodesE, nNodesX(1), nNodesX(2), nNodesX(3) ], &
                     nDOF, NodeNumberTable ), DatasetName, FILE_ID )

      END DO

    END DO

    CALL H5FCLOSE_F( FILE_ID, HDFERR )

    CALL H5CLOSE_F( HDFERR )

  END SUBROUTINE WriteRadiationFieldsHDF


  SUBROUTINE CreateGroupHDF( FileName, GroupName, FILE_ID )

    CHARACTER(len=*), INTENT(in) :: FileName
    CHARACTER(len=*), INTENT(in) :: GroupName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HID_T) :: GROUP_ID

    CALL H5GCREATE_F( FILE_ID, TRIM( GroupName ), GROUP_ID, HDFERR )

    CALL H5GCLOSE_F( GROUP_ID, HDFERR )

  END SUBROUTINE CreateGroupHDF


  SUBROUTINE WriteDataset1DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(1)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SIZE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 1, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    call H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    call H5SCLOSE_F( DATASPACE_ID, HDFERR )

    call H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset1DHDF


  SUBROUTINE WriteDataset3DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:,:,:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(3)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SHAPE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 3, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    CALL H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    CALL H5SCLOSE_F( DATASPACE_ID, HDFERR )

    CALL H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset3DHDF


  SUBROUTINE WriteDataset4DHDF( Dataset, DatasetName, FILE_ID )

    REAL(DP),         INTENT(in) :: Dataset(:,:,:,:)
    CHARACTER(LEN=*), INTENT(in) :: DatasetName
    INTEGER(HID_T),   INTENT(in) :: FILE_ID

    INTEGER(HSIZE_T) :: DATASIZE(4)
    INTEGER(HID_T)   :: DATASPACE_ID
    INTEGER(HID_T)   :: DATASET_ID

    DATASIZE = SHAPE( Dataset )

    CALL H5SCREATE_F( H5S_SIMPLE_F, DATASPACE_ID, HDFERR )

    CALL H5SSET_EXTENT_SIMPLE_F &
           ( DATASPACE_ID, 4, DATASIZE, DATASIZE, HDFERR )

    CALL H5DCREATE_F &
           ( FILE_ID, TRIM( DatasetName ), H5T_NATIVE_DOUBLE, &
             DATASPACE_ID, DATASET_ID, HDFERR )

    call H5DWRITE_F &
           ( DATASET_ID, H5T_NATIVE_DOUBLE, Dataset, DATASIZE, HDFERR )

    call H5SCLOSE_F( DATASPACE_ID, HDFERR )

    call H5DCLOSE_F( DATASET_ID, HDFERR )

  END SUBROUTINE WriteDataset4DHDF


END MODULE InputOutputModuleHDF
