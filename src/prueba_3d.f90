!=================================================================
! Prueba de las transformadas 3D de sphp (fftp3d_*): la 2D en
! (theta,phi) con SHTns, la transposición entre tareas y la FFT en
! z. Corre con cualquier cantidad de tareas hasta nz.
! (a) z periódica: f = sin(th) cos(ph) cos(2 pi q z/L)
!     + cos(th) sin(2 pi p z/L) + 1 tiene cinco coeficientes
!     exactos en (lm,kz); la ida y vuelta devuelve el campo (la 2D
!     de SHTns y la FFT en z vienen normalizadas).
! (b) FC-Gram en z (C = 25, d = 5, tablas de SPECTER): la derivada
!     en z de Y(th,ph) g(z), con Y = sin(th) cos(ph) + cos(th) y g
!     suave y no periódica en [0,1], contra la exacta en los planos
!     físicos, para nz = 64, 128 y 256.
!=================================================================

      PROGRAM prueba_3d
      USE fprecision
      USE commtypes
      USE mpivars
      USE fcgram
      IMPLICIT NONE

      INTEGER, PARAMETER :: lmax = 15, nlat = 3*(lmax+1)/2
      INTEGER, PARAMETER :: nphi = 2*nlat, nlm = (lmax+1)*(lmax+2)/2
      INTEGER, PARAMETER :: q = 1, p = 3, C = 25, od = 5
      REAL(KIND=GP), PARAMETER :: pi = 3.14159265358979323846_GP
      COMPLEX(KIND=GP), PARAMETER :: im = (0.0_GP,1.0_GP)
      CHARACTER(len=*), PARAMETER :: tdir = &
               '/home/mramdan/Doctorado/REDONDO/SPECTER/tables'

      TYPE(FCPLAN) :: plan
      TYPE(shtns_info), POINTER :: info
      REAL(C_DOUBLE), POINTER   :: ct(:),st(:)
! cth, sth y ph: cos(theta), sin(theta) y phi en la grilla
      REAL(KIND=GP), DIMENSION(nlat,nphi) :: cth,sth,ph
      REAL(KIND=GP), ALLOCATABLE    :: f(:,:,:),f0(:,:,:)
      COMPLEX(KIND=GP), ALLOCATABLE :: b(:,:),e(:,:)
      REAL(KIND=GP)    :: err(2),emax(2),dz,zk,kn,eant
      INTEGER :: nz,j,k,ic

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,nlm,nprocs,myrank,lmsta,lmend)

! (a) z periódica
      nz = 8
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL fcgram_create_plan(plan,(/nlat,nphi,nz/),(/0,0,0/), &
                              (/0,0,0/),(/1.0_GP,1.0_GP,1.0_GP/), &
                              '.',FFTW_ESTIMATE,lmax)
      CALL grilla()
      ALLOCATE( f(nlat,nphi,ksta:kend),f0(nlat,nphi,ksta:kend) )
      ALLOCATE( b(lmsta:lmend,nz),e(nlm,nz) )
      DO k = ksta,kend
         f0(:,:,k) = sth*COS(ph)*COS(2*pi*q*(k-1)/nz)+ &
                     cth*SIN(2*pi*p*(k-1)/nz)+1.0_GP
      END DO
      e = 0.0_GP
      e(lm(1,1),q+1)    = SQRT(2*pi/3)/2
      e(lm(1,1),nz-q+1) = SQRT(2*pi/3)/2
      e(lm(1,0),p+1)    = -im*SQRT(4*pi/3)/2
      e(lm(1,0),nz-p+1) = im*SQRT(4*pi/3)/2
      e(lm(0,0),1)      = SQRT(4*pi)
! La directa destruye su entrada: se transforma una copia
      f = f0
      CALL fftp3d_real_to_complex(plan,f,b,MPI_COMM_WORLD)
      err(1) = MAXVAL(ABS(b-e(lmsta:lmend,:)))
      CALL fftp3d_complex_to_real(plan,b,f,MPI_COMM_WORLD)
      err(2) = MAXVAL(ABS(f-f0))
      CALL MPI_REDUCE(err,emax,2,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         PRINT*, 'tareas:', nprocs
         PRINT*, '(a) z periódica, nz =', nz
         PRINT '(3X,A,ES10.2)', 'coeficientes, max|F-F_exacto|   :', emax(1)
         PRINT '(3X,A,ES10.2)', 'ida y vuelta, max|f-f_original|:', emax(2)
      ENDIF
      DEALLOCATE( f,f0,b,e )
      CALL fcgram_destroy_plan(plan)

! (b) FC-Gram en z: derivada en z de Y g(z)
      IF (myrank.eq.0) THEN
         PRINT*, '(b) FC-Gram en z (C = 25, d = 5): derivada en z, ', &
                 'max|df - df_exacta|/max|df_exacta|'
      ENDIF
      eant = 0.0_GP
      DO ic = 1,3
         nz = 32*2**ic
         dz = 1.0_GP/(nz-C-1)
         CALL range(1,nz,nprocs,myrank,ksta,kend)
         CALL fcgram_create_plan(plan,(/nlat,nphi,nz/),(/0,0,C/), &
                                 (/0,0,od/),(/1.0_GP,1.0_GP,dz/), &
                                 tdir,FFTW_ESTIMATE,lmax)
         CALL grilla()
         ALLOCATE( f(nlat,nphi,ksta:kend),f0(nlat,nphi,ksta:kend) )
         ALLOCATE( b(lmsta:lmend,nz) )
! Los planos físicos son los primeros nz-C; los demás los pone la
! continuación
         f = 0.0_GP
         f0 = 0.0_GP
         DO k = ksta,MIN(kend,nz-C)
            zk = (k-1)*dz
            f(:,:,k)  = (sth*COS(ph)+cth)*g(zk)
            f0(:,:,k) = (sth*COS(ph)+cth)*dg(zk)
         END DO
         CALL fftp3d_real_to_complex(plan,f,b,MPI_COMM_WORLD)
! Derivada: i kz, con kz = 2 pi k/(nz dz); el de Nyquist se anula
         DO j = 1,nz
            kn = j-1
            IF (j-1.gt.nz/2) kn = j-1-nz
            IF (j-1.eq.nz/2) kn = 0.0_GP
            b(:,j) = im*2*pi*kn/(nz*dz)*b(:,j)
         END DO
         CALL fftp3d_complex_to_real(plan,b,f,MPI_COMM_WORLD)
         err = 0.0_GP
         DO k = ksta,MIN(kend,nz-C)
            err(1) = MAX(err(1),MAXVAL(ABS(f(:,:,k)-f0(:,:,k))))
            err(2) = MAX(err(2),MAXVAL(ABS(f0(:,:,k))))
         END DO
         CALL MPI_REDUCE(err,emax,2,GC_REAL,MPI_MAX,0, &
                         MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            emax(1) = emax(1)/emax(2)
            IF (eant.gt.0) THEN
               PRINT '(3X,A,I4,A,ES10.2,A,F5.2)', 'nz =',nz,':', &
                     emax(1),'   orden',LOG(eant/emax(1))/LOG(2.0_GP)
            ELSE
               PRINT '(3X,A,I4,A,ES10.2)', 'nz =',nz,':',emax(1)
            ENDIF
            eant = emax(1)
         ENDIF
         DEALLOCATE( f,f0,b )
         CALL fcgram_destroy_plan(plan)
      END DO

      CALL MPI_FINALIZE(ierr)

      CONTAINS

!-----------------------------------------------------------------
! La grilla angular del plan: los nodos de Gauss los da SHTns
      SUBROUTINE grilla()
      CALL C_F_POINTER(plan%planrc%planrxy,info)
      CALL C_F_POINTER(info%ct,ct,(/nlat/))
      CALL C_F_POINTER(info%st,st,(/nlat/))
      DO j = 1,nphi
         cth(:,j) = ct
         sth(:,j) = st
         ph(:,j)  = 2*pi*(j-1)/nphi
      END DO
      END SUBROUTINE grilla

!-----------------------------------------------------------------
! Lugar del modo (l,m) en los arreglos de coeficientes
      INTEGER FUNCTION lm(l,m)
      INTEGER, INTENT(IN) :: l,m
      lm = shtns_lmidx(plan%planrc%planrxy,l,m)
      END FUNCTION lm

! g suave y no periódica en [0,1], y su derivada
      REAL(KIND=GP) FUNCTION g(z)
      REAL(KIND=GP), INTENT(IN) :: z
      g = SIN(2.3_GP*z+0.4_GP)+z**2
      END FUNCTION g

      REAL(KIND=GP) FUNCTION dg(z)
      REAL(KIND=GP), INTENT(IN) :: z
      dg = 2.3_GP*COS(2.3_GP*z+0.4_GP)+2*z
      END FUNCTION dg

      END PROGRAM prueba_3d
