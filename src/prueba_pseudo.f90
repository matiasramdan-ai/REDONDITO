!=================================================================
! Prueba de los operadores de pseudospec_hd.f90 en la cáscara
! 0.5 < r < 1.5, con lmax = 7, nz = 128 y FC-Gram en r (C = 25,
! d = 5, tablas de SPECTER). Sobre f = c_lm g(r), con
! g = exp(r) sin(3r) y un c_lm complejo distinto en cada lm, compara
! con las cuentas exactas en los puntos físicos:
!   derivr de orden 1 y 2 en (lm,r), y en (lm,kz) contra (lm,r),
!   laplang, multr, dlr, saxpby_c y copy, y prodvec en el espacio
!   físico.
! De fc_filter mide cuánto cambia un campo suave. Corre con
! cualquier cantidad de tareas y tiene que dar lo mismo.
!=================================================================

      PROGRAM prueba_pseudo
      USE fprecision
      USE commtypes
      USE mpivars
      USE fcgram
      USE fft
      USE grid
      USE boxsize
      USE kes
      USE class_GWorkspace3D
      IMPLICIT NONE

      CHARACTER(len=*), PARAMETER :: tdir = &
               '/home/mramdan/Doctorado/REDONDO/SPECTER/tables'
      TYPE(GWorkspace), TARGET :: ws
      COMPLEX(KIND=GP), ALLOCATABLE :: a(:,:),b(:,:),c(:,:),e(:,:)
      REAL(KIND=GP), ALLOCATABLE :: u(:,:,:,:),v(:,:,:,:),w(:,:,:,:)
      REAL(KIND=GP) :: err(9),emax(9),ct,st
      COMPLEX(KIND=GP) :: cl
      INTEGER :: n,i,j,k,lm

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL iniciar(7,128)
      n = nz-Cz
      ALLOCATE( a(lmsta:lmend,nz),b(lmsta:lmend,nz) )
      ALLOCATE( c(lmsta:lmend,nz),e(lmsta:lmend,nz) )
      DO k = 1,nz
         DO lm = lmsta,lmend
            cl = CMPLX(1.0_GP+0.1_GP*lm,0.5_GP-0.05_GP*lm,KIND=GP)
            a(lm,k) = cl*g(z(k),0)
         END DO
      END DO
      err = 0.0_GP

! derivr de orden 1 y 2 en (lm,r), contra la derivada exacta
      CALL derivr(a,b,1,'r')
      err(1) = relerr(b,dexacta(1))
      CALL derivr(a,e,2,'r')
      err(2) = relerr(e,dexacta(2))
! derivr en (lm,kz), con la ida y la vuelta a mano
      CALL copy(a,c)
      CALL fftp1d_real_to_complex_z(planfc,c,MPI_COMM_WORLD)
      CALL derivr(c,e,1,'k')
      CALL fftp1d_complex_to_real_z(planfc,e,MPI_COMM_WORLD)
      err(3) = relerr(e,b)
! laplang y multr son exactos
      CALL laplang(a,b)
      DO k = 1,nz
         DO lm = lmsta,lmend
            c(lm,k) = -ell(lm)*(ell(lm)+1)*a(lm,k)
         END DO
      END DO
      err(4) = relerr(b,c)
      CALL multr(a,b,-2)
      DO k = 1,nz
         c(:,k) = a(:,k)/z(k)**2
      END DO
      err(5) = relerr(b,c)
! dlr, contra d2g/dr2 - l(l+1) g/r**2
      CALL dlr(a,b)
      c = dexacta(2)
      DO k = 1,nz
         DO lm = lmsta,lmend
            c(lm,k) = c(lm,k)-ell(lm)*(ell(lm)+1)*a(lm,k)/z(k)**2
         END DO
      END DO
      err(6) = relerr(b,c)
! saxpby_c y copy
      CALL saxpby_c(c,a,2.0_GP,b,-3.0_GP)
      err(7) = relerr(c,2*a-3*b)
      CALL copy(b,c)
      err(7) = MAX(err(7),relerr(c,b))
! fc_filter sobre un campo suave
      CALL copy(a,c)
      CALL fc_filter(c,'r')
      err(8) = relerr(c,a)
! prodvec, contra el producto vectorial escrito a mano
      ALLOCATE( u(nlat,nphi,ksta:kend,3),v(nlat,nphi,ksta:kend,3) )
      ALLOCATE( w(nlat,nphi,ksta:kend,3) )
      DO k = ksta,kend
         DO j = 1,nphi
            DO i = 1,nlat
               ct = COS(theta(i))
               st = SIN(theta(i))
               u(i,j,k,:) = (/ st*COS(phi(j)),z(k),ct /)
               v(i,j,k,:) = (/ z(k)**2,SIN(phi(j)),st*ct /)
            END DO
         END DO
      END DO
      CALL prodvec(u(:,:,:,1),u(:,:,:,2),u(:,:,:,3),v(:,:,:,1), &
                   v(:,:,:,2),v(:,:,:,3),w(:,:,:,1),w(:,:,:,2), &
                   w(:,:,:,3))
      DO k = ksta,pkend
         DO j = 1,nphi
            DO i = 1,nlat
               err(9) = MAX(err(9),MAXVAL(ABS(w(i,j,k,:)- &
                        cruz(u(i,j,k,:),v(i,j,k,:)))))
            END DO
         END DO
      END DO

      CALL MPI_REDUCE(err,emax,9,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         PRINT*, 'tareas:', nprocs, '  lmax =', lmax, '  nz =', nz
         PRINT*, 'max|f - f_exacta|/max|f_exacta| en los puntos físicos:'
         PRINT '(3X,A,ES10.2)', 'derivr, orden 1           :', emax(1)
         PRINT '(3X,A,ES10.2)', 'derivr, orden 2           :', emax(2)
         PRINT '(3X,A,ES10.2)', 'derivr en kz contra en r  :', emax(3)
         PRINT '(3X,A,ES10.2)', 'laplang                   :', emax(4)
         PRINT '(3X,A,ES10.2)', 'multr                     :', emax(5)
         PRINT '(3X,A,ES10.2)', 'dlr                       :', emax(6)
         PRINT '(3X,A,ES10.2)', 'saxpby_c y copy           :', emax(7)
         PRINT '(3X,A,ES10.2)', 'fc_filter, cambio         :', emax(8)
         PRINT '(3X,A,ES10.2)', 'prodvec, max|error|       :', emax(9)
      ENDIF
      CALL terminar()
      CALL MPI_FINALIZE(ierr)

      CONTAINS

! g = exp(r) sin(3r) y sus derivadas primera y segunda
      REAL(KIND=GP) FUNCTION g(r,p)
      REAL(KIND=GP), INTENT(IN) :: r
      INTEGER, INTENT(IN) :: p
      SELECT CASE (p)
      CASE (0)
         g = EXP(r)*SIN(3*r)
      CASE (1)
         g = EXP(r)*(SIN(3*r)+3*COS(3*r))
      CASE DEFAULT
         g = EXP(r)*(-8*SIN(3*r)+6*COS(3*r))
      END SELECT
      END FUNCTION g

! c_lm por la derivada de orden p de g
      FUNCTION dexacta(p) RESULT(d)
      INTEGER, INTENT(IN) :: p
      COMPLEX(KIND=GP) :: d(lmsta:lmend,nz)
      INTEGER :: k,lm
      DO k = 1,nz
         DO lm = lmsta,lmend
            d(lm,k) = CMPLX(1.0_GP+0.1_GP*lm,0.5_GP-0.05_GP*lm, &
                            KIND=GP)*g(z(k),p)
         END DO
      END DO
      END FUNCTION dexacta

! Producto vectorial de dos vectores en (r,theta,phi)
      FUNCTION cruz(x,y) RESULT(c)
      REAL(KIND=GP), INTENT(IN) :: x(3),y(3)
      REAL(KIND=GP) :: c(3)
      c = (/ x(2)*y(3)-x(3)*y(2),x(3)*y(1)-x(1)*y(3), &
             x(1)*y(2)-x(2)*y(1) /)
      END FUNCTION cruz

! Error relativo en los puntos físicos, con los máximos de todas
! las tareas
      REAL(KIND=GP) FUNCTION relerr(x,y)
      COMPLEX(KIND=GP), INTENT(IN) :: x(lmsta:lmend,nz),y(lmsta:lmend,nz)
      REAL(KIND=GP) :: loc(2),glo(2)
      loc(1) = MAXVAL(ABS(x(:,1:nz-Cz)-y(:,1:nz-Cz)))
      loc(2) = MAXVAL(ABS(y(:,1:nz-Cz)))
      CALL MPI_ALLREDUCE(loc,glo,2,GC_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      relerr = glo(1)/glo(2)
      END FUNCTION relerr

! La namelist, la grilla, el plan y el pool
      SUBROUTINE iniciar(lm0,nz0)
      INTEGER, INTENT(IN) :: lm0,nz0
      IF (myrank.eq.0) THEN
         OPEN(1,FILE='prueba.inp',ACTION='write')
         WRITE(1,'(A,I0,A,I0,A)') '&grid lmax = ',lm0,', nz = ',nz0, &
                                  ', Cz = 25, oz = 5 /'
         WRITE(1,'(A)') '&boxparams rint = 0.5, rext = 1.5 /'
         CLOSE(1)
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      CALL grid_init('prueba.inp')
      CALL range(1,nlm,nprocs,myrank,lmsta,lmend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL box_init('prueba.inp')
      CALL fcgram_create_plan(planfc,(/nlat,nphi,nz/),(/0,0,Cz/), &
               (/0,0,oz/),(/1.0_GP,1.0_GP,z(2)-z(1)/),tdir, &
               FFTW_ESTIMATE,lmax)
      CALL ws%initialize_pool(1,3)
      END SUBROUTINE iniciar

! Borra la namelist y las grillas que escribe box_init
      SUBROUTINE terminar()
      CHARACTER(len=10), PARAMETER :: arch(4) = (/ 'prueba.inp', &
               'theta.txt ','phi.txt   ','z.txt     ' /)
      INTEGER :: i
      IF (myrank.eq.0) THEN
         DO i = 1,4
            OPEN(1,FILE=TRIM(arch(i)))
            CLOSE(1,STATUS='delete')
         END DO
      ENDIF
      END SUBROUTINE terminar

      END PROGRAM prueba_pseudo
