!=================================================================
! MODULES for FFTP v3
! Parallel Fast Fourier Transform in 2D and 3D
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!=================================================================

!=================================================================

  MODULE fftplans
!
! Set the variable ikind to:  4 in 32 bits machines
!                             8 in 64 bits machines
! Set the variable csize to:  8 if L1 cache is <= 64 kb
!                            16 if L1 cache is 128 kb
! The variable nstrip controls strip mining during the 
! transposition. Often set to 1.
!
      USE fprecision
      USE mpi_f08, only   :  MPI_Datatype
      USE, INTRINSIC :: iso_c_binding

      IMPLICIT NONE

      INCLUDE 'fftw3.f'
#include "gfft_defs.h"
! Interfaces de SHTns, para las transformadas angulares
      INCLUDE 'shtns.f03'
 
      INTEGER, PARAMETER  :: ikind = IKIND_
      INTEGER, PARAMETER  :: csize = CSIZE_
      INTEGER, PARAMETER  :: nstrip = NSTRIP_
      INTEGER, PARAMETER  :: FFTW_REAL_TO_COMPLEX = FFTW_FORWARD
      INTEGER, PARAMETER  :: FFTW_COMPLEX_TO_REAL = FFTW_BACKWARD
      INTEGER             :: hcom,hfft,hmem,htra,hcont,htot
      DOUBLE PRECISION    :: comtime = 0.0
      DOUBLE PRECISION    :: ffttime = 0.0
      DOUBLE PRECISION    :: tratime = 0.0
      DOUBLE PRECISION    :: conttime = 0.0
      DOUBLE PRECISION    :: tottime = 0.0
      TYPE FFTPLAN
         COMPLEX(KIND=GP), DIMENSION (:,:), POINTER  :: ccarr
         COMPLEX(KIND=GP), DIMENSION (:,:), POINTER  :: carr1,carr2,carr3
         REAL(KIND=GP),    DIMENSION (:,:,:), POINTER  :: rarr
         TYPE(MPI_Datatype), DIMENSION (:), POINTER :: itype1, itype2
         ! planrxy: el plan de SHTns para la transformada 2D en
         ! (theta,phi), en lugar del de la FFT 2D de FFTW
         TYPE(C_PTR)         :: planrxy
         INTEGER(kind=ikind) :: plancz      ! host (FFTW) plan
         GFFT_HANDLE_T       :: dplancz     ! device plan
         INTEGER :: nlat,nphi,nz,nlm,lmax
      END TYPE FFTPLAN
!
! Arreglos del camino en la GPU, compartidos por todos los planes (la
! directa y la inversa nunca se usan a la vez) y residentes en la
! GPU: la salida de las transformadas 2D (gcarr1, gcarr2 y gcarr3,
! uno por componente en las vectoriales; las escalares usan gcarr1)
! y los mensajes empaquetados (gsbuf). No hay gc1: lo que llega de
! cada tarea es contiguo en la salida (lm,z) y va directo ahí. Las
! tablas tienen los rangos de lm (glmoff/glmlen) y de z
! (gkoff/gklen) de cada tarea, y los desplazamientos y tamaños de
! los bloques que se mandan empaquetados (gsoff/gscnt) y de los
! bloques contiguos que se reciben en la salida (groff/grcnt).
! gexchange elige el intercambio entre tareas:
!   1: MPI recibe los arreglos de la GPU (MPI GPU-aware, el default)
!   2: pasando por las copias en la CPU (para MPI que no es
!      GPU-aware; variable de entorno GHOST_GPU_EXCHANGE=staged)
      COMPLEX(KIND=GP), DIMENSION (:,:), ALLOCATABLE, TARGET :: gcarr1
      COMPLEX(KIND=GP), DIMENSION (:,:), ALLOCATABLE, TARGET :: gcarr2
      COMPLEX(KIND=GP), DIMENSION (:,:), ALLOCATABLE, TARGET :: gcarr3
      COMPLEX(KIND=GP), DIMENSION (:)  , ALLOCATABLE, TARGET :: gsbuf
      INTEGER, DIMENSION (:), ALLOCATABLE :: glmoff,glmlen,gkoff,gklen
      INTEGER, DIMENSION (:), ALLOCATABLE :: gsoff,gscnt,groff,grcnt
      INTEGER :: gexchange = 1
      LOGICAL :: gbuffers_ready = .FALSE.
      SAVE

  END MODULE fftplans
!=================================================================

  MODULE threads
!
! nth: number of threads used in OpenMP loops and FFTs
!$    USE omp_lib
      INTEGER :: nth
      SAVE

  END MODULE threads
!=================================================================

  MODULE mpivars
      ! pkend is last k-index corresponding to the
      ! actual physical (non-continuated) domain.
      INTEGER, SAVE :: ista,iend
      INTEGER, SAVE :: jsta,jend
      INTEGER, SAVE :: ksta,kend,pkend
      INTEGER, SAVE :: lmsta,lmend
      INTEGER, SAVE :: nprocs,myrank
      INTEGER, SAVE :: provided
      INTEGER, SAVE :: ierr

  END MODULE mpivars
!=================================================================
