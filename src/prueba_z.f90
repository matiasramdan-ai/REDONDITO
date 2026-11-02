!=================================================================
! Prueba de la FFT en z de sphp (fftp1d_*_z), sobre arreglos
! (lm,z): cada tarea tiene sus lm y todos los planos, así que la
! FFT en z no comunica. Corre con cualquier cantidad de tareas.
! (a) z periódica: c_lm cos(2 pi q z/L) + d_lm sin(2 pi p z/L)
!     tiene cuatro coeficientes exactos en kz, y la ida y vuelta
!     devuelve el campo (la directa normaliza).
! (b) FC-Gram en z (C = 25, d = 5, tablas de SPECTER): la derivada
!     de c_lm g(z), con g suave y no periódica en [0,1], contra la
!     exacta en los puntos físicos, para nz = 64, 128 y 256.
!=================================================================

      PROGRAM prueba_z
      USE fprecision
      USE commtypes
      USE mpivars
      USE fcgram
      IMPLICIT NONE

      INTEGER, PARAMETER :: lmax = 7, nlat = 3*(lmax+1)/2
      INTEGER, PARAMETER :: nphi = 2*nlat, nlm = (lmax+1)*(lmax+2)/2
      INTEGER, PARAMETER :: q = 3, p = 5, C = 25, od = 5
      REAL(KIND=GP), PARAMETER :: pi = 3.14159265358979323846_GP
      CHARACTER(len=*), PARAMETER :: tdir = &
               '/home/mramdan/Doctorado/REDONDO/SPECTER/tables'

      TYPE(FCPLAN) :: plan
      COMPLEX(KIND=GP), ALLOCATABLE :: f(:,:),e(:,:)
      COMPLEX(KIND=GP) :: cl(nlm),dl(nlm)
      COMPLEX(KIND=GP), PARAMETER :: im = (0.0_GP,1.0_GP)
      REAL(KIND=GP)    :: err(2),emax(2),dz,zj,kn,eant
      INTEGER :: nz,lm,j,ic

! Coeficientes distintos para cada lm
      DO lm = 1,nlm
         cl(lm) = CMPLX(1.0_GP*lm/nlm,-0.5_GP,KIND=GP)
         dl(lm) = CMPLX(0.3_GP,1.0_GP-1.0_GP*lm/nlm,KIND=GP)
      END DO

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,nlm,nprocs,myrank,lmsta,lmend)

! (a) z periódica
      nz = 16
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL fcgram_create_plan(plan,(/nlat,nphi,nz/),(/0,0,0/), &
                              (/0,0,0/),(/1.0_GP,1.0_GP,1.0_GP/), &
                              '.',FFTW_ESTIMATE,lmax)
      ALLOCATE( f(lmsta:lmend,nz),e(lmsta:lmend,nz) )
      DO j = 1,nz
         f(:,j) = cl(lmsta:lmend)*COS(2*pi*q*(j-1)/nz)+ &
                  dl(lmsta:lmend)*SIN(2*pi*p*(j-1)/nz)
      END DO
      e = 0.0_GP
      e(:,q+1)    = cl(lmsta:lmend)/2
      e(:,nz-q+1) = cl(lmsta:lmend)/2
      e(:,p+1)    = -im*dl(lmsta:lmend)/2
      e(:,nz-p+1) = im*dl(lmsta:lmend)/2
      CALL fftp1d_real_to_complex_z(plan,f,MPI_COMM_WORLD)
      err(1) = MAXVAL(ABS(f-e))
      CALL fftp1d_complex_to_real_z(plan,f,MPI_COMM_WORLD)
      DO j = 1,nz
         e(:,j) = cl(lmsta:lmend)*COS(2*pi*q*(j-1)/nz)+ &
                  dl(lmsta:lmend)*SIN(2*pi*p*(j-1)/nz)
      END DO
      err(2) = MAXVAL(ABS(f-e))
      CALL MPI_REDUCE(err,emax,2,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         PRINT*, 'tareas:', nprocs
         PRINT*, '(a) z periódica, nz =', nz
         PRINT '(3X,A,ES10.2)', 'coeficientes, max|F-F_exacto|   :', emax(1)
         PRINT '(3X,A,ES10.2)', 'ida y vuelta, max|f-f_original|:', emax(2)
      ENDIF
      DEALLOCATE( f,e )
      CALL fcgram_destroy_plan(plan)

! (b) FC-Gram en z: derivada de c_lm g(z)
      IF (myrank.eq.0) THEN
         PRINT*, '(b) FC-Gram en z (C = 25, d = 5): derivada, ', &
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
         ALLOCATE( f(lmsta:lmend,nz),e(lmsta:lmend,nz) )
! Los puntos físicos son los primeros nz-C; los demás los pone la
! continuación
         f = 0.0_GP
         DO j = 1,nz-C
            zj = (j-1)*dz
            f(:,j) = cl(lmsta:lmend)*g(zj)
            e(:,j) = cl(lmsta:lmend)*dg(zj)
         END DO
         CALL fftp1d_real_to_complex_z(plan,f,MPI_COMM_WORLD)
! Derivada: i kz, con kz = 2 pi k/(nz dz); el de Nyquist se anula
         DO j = 1,nz
            kn = j-1
            IF (j-1.gt.nz/2) kn = j-1-nz
            IF (j-1.eq.nz/2) kn = 0.0_GP
            f(:,j) = im*2*pi*kn/(nz*dz)*f(:,j)
         END DO
         CALL fftp1d_complex_to_real_z(plan,f,MPI_COMM_WORLD)
         err(1) = MAXVAL(ABS(f(:,1:nz-C)-e(:,1:nz-C)))/ &
                  MAXVAL(ABS(e(:,1:nz-C)))
         CALL MPI_REDUCE(err,emax,1,GC_REAL,MPI_MAX,0, &
                         MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) THEN
            IF (eant.gt.0) THEN
               PRINT '(3X,A,I4,A,ES10.2,A,F5.2)', 'nz =',nz,':', &
                     emax(1),'   orden',LOG(eant/emax(1))/LOG(2.0_GP)
            ELSE
               PRINT '(3X,A,I4,A,ES10.2)', 'nz =',nz,':',emax(1)
            ENDIF
         ENDIF
         eant = emax(1)
         DEALLOCATE( f,e )
         CALL fcgram_destroy_plan(plan)
      END DO

      CALL MPI_FINALIZE(ierr)

      CONTAINS

! g suave y no periódica en [0,1], y su derivada
      REAL(KIND=GP) FUNCTION g(z)
      REAL(KIND=GP), INTENT(IN) :: z
      g = SIN(2.3_GP*z+0.4_GP)+z**2
      END FUNCTION g

      REAL(KIND=GP) FUNCTION dg(z)
      REAL(KIND=GP), INTENT(IN) :: z
      dg = 2.3_GP*COS(2.3_GP*z+0.4_GP)+2*z
      END FUNCTION dg

      END PROGRAM prueba_z
