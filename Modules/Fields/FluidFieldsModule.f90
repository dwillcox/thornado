MODULE FluidFieldsModule

  USE KindModule, ONLY: &
    DP
  USE ProgramHeaderModule, ONLY: &
    nDOFX

  IMPLICIT NONE
  PRIVATE

  ! --- Conserved Fluid Fields ---

  INTEGER, PUBLIC, PARAMETER :: iCF_D  = 1  ! Conserved Baryon Density
  INTEGER, PUBLIC, PARAMETER :: iCF_S1 = 2  ! Conserved Momentum Density 1
  INTEGER, PUBLIC, PARAMETER :: iCF_S2 = 3  ! Conserved Momentum Density 2
  INTEGER, PUBLIC, PARAMETER :: iCF_S3 = 4  ! Conserved Momentum Density 3
  INTEGER, PUBLIC, PARAMETER :: iCF_E  = 5  ! Conserved Energy Density
  INTEGER, PUBLIC, PARAMETER :: nCF    = 5  ! n Conserved Fluid Fields

  CHARACTER(32), DIMENSION(nCF), PUBLIC, PARAMETER :: &
    namesCF = [ 'Conserved Baryon Density        ', &
                'Conserved Momentum Density (1)  ', &
                'Conserved Momentum Density (2)  ', &
                'Conserved Momentum Density (3)  ', &
                'Conserved Energy Density        ' ]

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: uCF, rhsCF

  ! --- Primitive Fluid Fields ---

  INTEGER, PUBLIC, PARAMETER :: iPF_D  = 1  ! Comoving Baryon Density
  INTEGER, PUBLIC, PARAMETER :: iPF_V1 = 2  ! Three-Velocity 1
  INTEGER, PUBLIC, PARAMETER :: iPF_V2 = 3  ! Three-Velocity 2
  INTEGER, PUBLIC, PARAMETER :: iPF_V3 = 4  ! Three-Velocity 3
  INTEGER, PUBLIC, PARAMETER :: iPF_E  = 5  ! Internal Energy Density
  INTEGER, PUBLIC, PARAMETER :: nPF    = 5  ! n Primitive Fluid Fields

  CHARACTER(32), DIMENSION(nPF), PUBLIC, PARAMETER :: &
    namesPF = [ 'Comoving Baryon Density         ', &
                'Three-Velocity (1)              ', &
                'Three-Velocity (2)              ', &
                'Three-Velocity (3)              ', &
                'Internal Energy Density         ' ]

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: uPF

  ! --- Auxiliary Fluid Fields ---

  INTEGER, PUBLIC, PARAMETER :: iAF_P  = 1  ! Pressure
  INTEGER, PUBLIC, PARAMETER :: iAF_T  = 2  ! Temperature
  INTEGER, PUBLIC, PARAMETER :: iAF_Ye = 3  ! Electron Fraction
  INTEGER, PUBLIC, PARAMETER :: iAF_S  = 4  ! Entropy Per Baryon
  INTEGER, PUBLIC, PARAMETER :: iAF_E  = 5  ! Specific Internal Energy
  INTEGER, PUBLIC, PARAMETER :: iAF_Me = 6  ! Electron Chemical Potential
  INTEGER, PUBLIC, PARAMETER :: iAF_Mp = 7  ! Proton Chemical Potential
  INTEGER, PUBLIC, PARAMETER :: iAF_Mn = 8  ! Neutron Chemical Potential
  INTEGER, PUBLIC, PARAMETER :: iAF_Gm = 9  ! Ratio of Specific Heats
  INTEGER, PUBLIC, PARAMETER :: iAF_Cs = 10 ! Sound Speed
  INTEGER, PUBLIC, PARAMETER :: nAF    = 10 ! n Auxiliary Fluid Fields

  CHARACTER(32), DIMENSION(nAF), PUBLIC, PARAMETER :: &
    namesAF = [ 'Pressure                        ', &
                'Temperature                     ', &
                'Electron Fraction               ', &
                'Entropy Per Baryon              ', &
                'Specific Internal Energy        ', &
                'Electron Chemical Potential     ', &
                'Proton Chemical Potential       ', &
                'Neutron Chemical Potential      ', &
                'Ratio of Specific Heats (Gamma) ', &
                'Sound Speed                     ' ]

  REAL(DP), DIMENSION(:,:,:,:,:), ALLOCATABLE, PUBLIC :: uAF

  PUBLIC :: CreateFluidFields
  PUBLIC :: DestroyFluidFields

CONTAINS


  SUBROUTINE CreateFluidFields( nX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX

    CALL CreateFluidFields_Conserved( nX, swX )
    CALL CreateFluidFields_Primitive( nX, swX )
    CALL CreateFluidFields_Auxiliary( nX, swX )

  END SUBROUTINE CreateFluidFields


  SUBROUTINE CreateFluidFields_Conserved( nX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX

    INTEGER :: iCF

    WRITE(*,*)
    WRITE(*,'(A5,A24)') '', 'Fluid Fields (Conserved)'
    WRITE(*,*)
    DO iCF = 1, nCF
      WRITE(*,'(A5,A32)') '', TRIM( namesCF(iCF) )
    END DO

    ALLOCATE( uCF &
                (1:nDOFX, &
                 1-swX(1):nX(1)+swX(1), &
                 1-swX(2):nX(2)+swX(2), &
                 1-swX(3):nX(3)+swX(3), &
                 1:nCF) )

    ALLOCATE( rhsCF &
                (1:nDOFX, &
                 1-swX(1):nX(1)+swX(1), &
                 1-swX(2):nX(2)+swX(2), &
                 1-swX(3):nX(3)+swX(3), &
                 1:nCF) )

  END SUBROUTINE CreateFluidFields_Conserved


  SUBROUTINE CreateFluidFields_Primitive( nX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX

    INTEGER :: iPF

    WRITE(*,*)
    WRITE(*,'(A5,A24)') '', 'Fluid Fields (Primitive)'
    WRITE(*,*)
    DO iPF = 1, nPF
      WRITE(*,'(A5,A32)') '', TRIM( namesPF(iPF) )
    END DO

    ALLOCATE( uPF(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nPF) )

  END SUBROUTINE CreateFluidFields_Primitive


  SUBROUTINE CreateFluidFields_Auxiliary( nX, swX )

    INTEGER, DIMENSION(3), INTENT(in) :: nX, swX

    INTEGER :: iAF

    WRITE(*,*)
    WRITE(*,'(A5,A24)') '', 'Fluid Fields (Auxiliary)'
    WRITE(*,*)
    DO iAF = 1, nAF
      WRITE(*,'(A5,A32)') '', TRIM( namesAF(iAF) )
    END DO

    ALLOCATE( uAF(1:nDOFX, &
                  1-swX(1):nX(1)+swX(1), &
                  1-swX(2):nX(2)+swX(2), &
                  1-swX(3):nX(3)+swX(3), &
                  1:nAF) )

  END SUBROUTINE CreateFluidFields_Auxiliary


  SUBROUTINE DestroyFluidFields

    DEALLOCATE( uCF, rhsCF, uPF, uAF )

  END SUBROUTINE DestroyFluidFields


END MODULE FluidFieldsModule
