!=================================================================
! MODULES for 3D codes
!
! 2003 Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar
!=================================================================

!=================================================================

  MODULE grid
      USE fprecision
!
! lmax: grado máximo de los armónicos esféricos, y nlm: número de
! coeficientes (l,m) con m >= 0, (lmax+1)(lmax+2)/2
      INTEGER :: lmax
      INTEGER :: nlm
!
! nlat, nphi: número de puntos de la grilla angular, de Gauss en
! theta y uniforme en phi, que salen de lmax con la regla de 3/2
! (ver grid_init)
      INTEGER :: nlat
      INTEGER :: nphi
!
! nz: number of points in the radial grid (z is the radius),
! including the continuation points
      INTEGER :: nz
!
! Cz: number of continuation points, and oz: order at the
! boundaries, in z
      INTEGER :: Cz
      INTEGER :: oz
!
! Unidimensional grids
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:) :: theta,phi,z
      SAVE
    CONTAINS
!-----------------------------------------------------------------
      SUBROUTINE grid_init(infile_)
!
! Initializes the grid spatial resolution and the FC-Gram
! continuation. Cz and oz are zero if they are not in the
! namelist. La grilla angular sale de lmax con la regla de 3/2
! (sin aliasing en el producto de dos campos): nlat es el menor
! múltiplo de 4 que no baja de 3(lmax+1)/2 (SHTns en la GPU pide
! múltiplos de 4) y nphi = 2 nlat.
      USE commtypes
      USE mpivars
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: infile_
      NAMELIST/ grid / lmax,nz,Cz,oz
      Cz = 0
      oz = 0
      IF ( myrank .EQ. 0 ) THEN
         OPEN(1,file=infile_,status='unknown',form="formatted")
         READ(1,NML=grid)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(lmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(nz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(Cz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(oz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      nlm  = (lmax+1)*(lmax+2)/2
      nlat = 4*((3*(lmax+1)+7)/8)
      nphi = 2*nlat
      END SUBROUTINE grid_init

  END MODULE grid

!=================================================================

  MODULE filefmt
!
! Change the length of the string 'ext' to change the number
! of characters used to number the binary files and files with
! the spectra. The format fmtext should be consistent with the
! length of the string, e.g. for a 5 character extent
! len=7 and fmtext = '(i5.5)'.
      CHARACTER(len=4)      :: ext
      CHARACTER(len=6),SAVE :: fmtext = '(i4.4)'

  END MODULE filefmt
!=================================================================

  MODULE fft
      USE fcgram
      TYPE(FCPLAN), SAVE, TARGET :: planfc

  END MODULE fft
!=================================================================

  MODULE var
      USE fprecision
      REAL(KIND=GP), PARAMETER    :: pi = 3.14159265358979323846_GP
      COMPLEX(KIND=GP), PARAMETER :: im = (0.0_GP,1.0_GP)
      SAVE

  END MODULE var
!=================================================================

  MODULE kes
!
! kz: wavenumbers in z; ell, emm: grado l y orden m de cada índice
! lm, en el orden de SHTns (ver box_init)
      USE fprecision
      REAL(KIND=GP), ALLOCATABLE, DIMENSION (:) :: kz
      INTEGER, ALLOCATABLE, DIMENSION (:)       :: ell,emm
      SAVE
  END MODULE kes

!=================================================================

  MODULE boxsize
!
! rint, rext: radios interno y externo de la cáscara, en unidades
! de 1. Dkz is the spacing between wavenumbers in z.
      USE fprecision
      REAL(KIND=GP)    :: rint,rext
      REAL(KIND=GP)    :: Dkz
    CONTAINS

!-----------------------------------------------------------------
      SUBROUTINE box_init(infile_)
!
! Reads the radii from the namelist 'boxparams' and builds the
! grids in physical and spectral space. La grilla radial z va de
! rint (z(1)) a rext (z(nz-Cz)), uniforme, seguida por los Cz puntos
! de la continuación de FC-Gram.
!     rint : radio interno
!     rext : radio externo
      USE commtypes
      USE mpivars
      USE fftplans
      USE grid
      USE kes
      USE var
      USE gmem
!$    USE threads
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: infile_
      TYPE(C_PTR)                  :: sht
      TYPE(shtns_info), POINTER    :: info
      REAL(C_DOUBLE), POINTER      :: ct(:)
      INTEGER(C_SHORT), POINTER    :: li(:),mi(:)
      REAL(KIND=GP)                :: rmp
      INTEGER                      :: i,j,k
      NAMELIST / boxparams / rint,rext

      IF (myrank.eq.0) THEN
         OPEN(1,file=infile_,status='unknown',form="formatted")
         READ(1,NML=boxparams)
         CLOSE(1)
      ENDIF
      CALL MPI_BCAST(rint,1,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(rext,1,GC_REAL,0,MPI_COMM_WORLD,ierr)

      ! Grids and wavenumbers (with device copies in GPU builds)
      CALL galloc( theta, nlat )
      CALL galloc( phi, nphi )
      CALL galloc( z, nz )
      CALL galloc( kz, nz )
      CALL galloc( ell, nlm )
      CALL galloc( emm, nlm )

! Construct grids
! Los nodos de Gauss en theta y el (l,m) de cada índice lm salen de
! una configuración de SHTns que se usa solo acá, con la grilla del
! plan de sphp (fftp3d_create_plan), así box_init se puede llamar
! antes de crear el plan, como en SPECTER
      sht = shtns_create(lmax,lmax,1,SHT_ORTHONORMAL+SHT_NO_CS_PHASE)
      CALL shtns_set_grid(sht,SHT_GAUSS+SHT_THETA_CONTIGUOUS,1.0D-14, &
                          nlat,nphi)
      CALL C_F_POINTER(sht,info)
      CALL C_F_POINTER(info%ct,ct,(/nlat/))
      CALL C_F_POINTER(info%li,li,(/nlm/))
      CALL C_F_POINTER(info%mi,mi,(/nlm/))
      DO i = 1,nlat
         theta(i) = ACOS(ct(i))
      END DO
      DO i = 1,nlm
         ell(i) = li(i)
         emm(i) = mi(i)
      END DO
      CALL shtns_destroy(sht)

      rmp = 2.0_GP*pi/nphi
      DO j = 1,nphi
        phi(j) = rmp*(j-1)
      ENDDO

      ! Non-periodic in z
      rmp = (rext-rint)/(nz-Cz-1)
      Dkz = 2.0_GP*pi/(rmp*nz)
      DO k=1,nz
        z(k) = rint+rmp*(k-1)
      ENDDO
      pkend = min(nz-Cz,kend)  ! kend of physical domain

      ! Save grid
      IF ( myrank .eq. 0 ) THEN
         OPEN(1,file='theta.txt',action='write')
         OPEN(2,file='phi.txt',action='write')
         OPEN(3,file='z.txt',action='write')
         DO i = 1,nlat
            WRITE(1,FMT='(1P E23.15)') theta(i)
         ENDDO
         DO i = 1,nphi
            WRITE(2,FMT='(1P E23.15)') phi(i)
         ENDDO
         DO i = 1,nz-Cz
            WRITE(3,FMT='(1P E23.15)') z(i)
         ENDDO
         CLOSE(1)
         CLOSE(2)
         CLOSE(3)
      ENDIF

! Construct grids for wavenumber domain
      DO k = 1,nz/2
         kz(k) = real(k-1,kind=GP)
         kz(k+nz/2) = real(k-nz/2-1,kind=GP)
      END DO
      kz = kz*Dkz

      ! Device copies (no-op in CPU builds)
      CALL gupdate_to( theta )
      CALL gupdate_to( phi )
      CALL gupdate_to( z )
      CALL gupdate_to( kz )
      CALL gupdate_to( ell )
      CALL gupdate_to( emm )
      END SUBROUTINE box_init

  END MODULE boxsize

!=================================================================

  MODULE random
      USE var
      USE fprecision
      CONTAINS
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       REAL(KIND=GP) FUNCTION randu(idum)
!
! Uniform distributed random numbers between -1 and
! 1. The seed idum must be between 0 and the value
! of mask

       INTEGER, PARAMETER :: iq=127773,ir=2836,mask=123459876
       INTEGER, PARAMETER :: ia=16807,im=2147483647
       INTEGER            :: k,idum
       REAL(KIND=GP), PARAMETER :: am=1./im

       idum = ieor(idum,mask)
       k = idum/iq
       idum = ia*(idum-k*iq)-ir*k
       IF (idum.lt.0) idum = idum+im
       randu = am*idum
       randu = (randu-.5)*2
       idum = ieor(idum,mask)
       END FUNCTION randu
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       REAL(KIND=GP) FUNCTION randn(idum)
!
! Normally distributed random numbers with zero mean
! and unit variance. The seed idum must be between 0
! and the value of mask in randu.

       REAL(KIND=GP)      :: v1,v2,ran1
       REAL(KIND=GP)      :: fac,rsq
       REAL(KIND=GP),SAVE :: gset
       INTEGER       :: idum
       INTEGER, SAVE :: iset

       IF ((iset.ne.0).or.(iset.ne.1)) iset=0
       IF (idum.lt.0) iset=0
       IF (iset.eq.0) THEN
          rsq = 2.
          DO WHILE ((rsq.ge.1.).or.(rsq.eq.0.))
             v1 = randu(idum)
             v2 = randu(idum)
             rsq = v1**2+v2**2
          END DO
          fac = sqrt(-2.*log(rsq)/rsq)
          gset = v1*fac
          randn = v2*fac
          iset = 1
       ELSE
          randn = gset
          iset = 0
       ENDIF
       END FUNCTION randn
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       SUBROUTINE randn_cmplx(x,y,idum)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Normally distributed complex random numbers with
! unit variance. The seed idum must be between 0
! and the value of mask in randu.
!-----------------------------------------------------------------

       REAL(KIND=GP), INTENT (OUT) :: x,y
       REAL(KIND=GP)       :: fac,phi
       INTEGER, INTENT(IN) :: idum

       x = randu(idum)
       y = randu(idum)
       do while (y.eq.0)
          y = randu(idum)
       end do
       phi = 2.*pi*abs(x)
       fac = sqrt(-2.*log(abs(y)))
       x = fac*cos(phi)
       y = fac*sin(phi)

       END SUBROUTINE randn_cmplx
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       SUBROUTINE prandom_seed(iseed)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Computes a unique seed for each MPI task. Note: there is a limit
! of 2^(size(integer)) for the number of tasks that can have a
! unique seed. iseed must be an integer seed, and iseed+myrank
! must not be the value of mask in randu (randu stays at that
! value, and all the words of the seed would be equal).
! (Fixed as in GHOST, rama new_paradigm: the old array constructor
! had size(seed)+1 elements, an out-of-bounds write, and seeds with
! many small or correlated words give biased first numbers.)
!
! Parameters
!     iseed : integer seed
!-----------------------------------------------------------------
       USE fprecision
       USE commtypes
       USE mpivars

       INTEGER, INTENT(INOUT)   :: iseed
       INTEGER, ALLOCATABLE     :: iseed1(:)
       INTEGER                  :: j,k,idum
       REAL(KIND=GP)            :: r

       iseed     = mod(iseed+myrank,abs(huge(0)-iseed)-1)
       CALL random_seed(size=k)
       ALLOCATE (iseed1(k))
! Fills the seed with successive states of randu, since some
! compilers (nvfortran) do not scramble the seed, and seeds with
! many small or correlated words give biased first numbers. The
! first states are discarded because the seeds of adjacent ranks
! start close.
       idum = iseed
       DO j = 1,3
          r = randu(idum)
       END DO
       DO j = 1,k
          r = randu(idum)
          iseed1(j) = idum
       END DO
       CALL random_seed(put=iseed1)
       DEALLOCATE(iseed1)

       END SUBROUTINE prandom_seed
!-----------------------------------------------------------------
!-----------------------------------------------------------------

       SUBROUTINE prandom_number(r)
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!
! Generates a random number for this MPI task on the interval [0,1].
! Seed for this must be generated by using prandom_seed call.
!
! Parameters
!     iseed : integer seed
!-----------------------------------------------------------------
       USE fprecision
       USE commtypes
       USE mpivars

       REAL(KIND=GP), INTENT(OUT)   :: r

       CALL random_number(r)

       END SUBROUTINE prandom_number
!-----------------------------------------------------------------
!-----------------------------------------------------------------

  END MODULE random
!=================================================================
