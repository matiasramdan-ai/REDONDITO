!=================================================================
! Prueba de boundary en la cáscara 0.5 < r < 1.5, con lmax = 7,
! nz = 128 y FC-Gram en r (C = 25, d = 5, tablas de SPECTER).
! setup_bc arma las matrices B_l con no deslizamiento en las dos
! paredes, y v_imposebc se aplica a lo que dejaría una etapa:
!   Z = c_lm (1 + exp(r) sin(3r)), que no se anula en las paredes;
!   phi = D_l W0, con W0 = c_lm (r-r1)**2 (r-r2)**2 exp(r), que
!   cumple las cuatro condiciones, con valores agregados en las
!   cuatro filas de borde y en el modo l = 0.
! Mide las condiciones, que lo que no son paredes no cambie, que
! aplicarla dos veces dé lo mismo, y W contra W0 (error de truncado).
! Corre con cualquier cantidad de tareas y tiene que dar lo mismo.
!=================================================================

      PROGRAM prueba_boundary
      USE fprecision
      USE commtypes
      USE mpivars
      USE fcgram
      USE fft
      USE grid
      USE boxsize
      USE kes
      USE class_GWorkspace3D
      USE boundary
      IMPLICIT NONE

      CHARACTER(len=*), PARAMETER :: tdir = &
               '/home/mramdan/Doctorado/REDONDO/SPECTER/tables'
      TYPE(GWorkspace), TARGET :: ws
      TYPE(BCPLAN) :: bcp
      CHARACTER(len=15) :: bck(2)
      COMPLEX(KIND=GP), ALLOCATABLE :: z0(:,:),w0(:,:),f0(:,:)
      COMPLEX(KIND=GP), ALLOCATABLE :: tor(:,:),fi(:,:),pol(:,:),pol2(:,:)
      REAL(KIND=GP), ALLOCATABLE :: d1(:,:)
      REAL(KIND=GP) :: err(8),emax(8),dmax(2),gmax(2),r1,r2,wm
      COMPLEX(KIND=GP) :: cl
      INTEGER :: n,k,lm

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL iniciar(7,128)
      n = nz-Cz
      r1 = z(1)
      r2 = z(n)
      bck(1) = 'noslip'
      bck(2) = 'noslip'
      CALL setup_bc(bcp,planfc,bck,'v')
      ALLOCATE( z0(lmsta:lmend,nz),w0(lmsta:lmend,nz),f0(lmsta:lmend,nz) )
      ALLOCATE( tor(lmsta:lmend,nz),fi(lmsta:lmend,nz) )
      ALLOCATE( pol(lmsta:lmend,nz),pol2(lmsta:lmend,nz),d1(n,n) )
      DO k = 1,nz
         DO lm = lmsta,lmend
            cl = CMPLX(1.0_GP+0.1_GP*lm,0.5_GP-0.05_GP*lm,KIND=GP)
            z0(lm,k) = cl*(1.0_GP+EXP(z(k))*SIN(3*z(k)))
            w0(lm,k) = cl*(z(k)-r1)**2*(z(k)-r2)**2*EXP(z(k))
         END DO
      END DO
      CALL dlr(w0,f0)
! Lo que dejaría una etapa en las filas de borde y en l = 0
      tor = z0
      fi = f0
      fi(:,1) = fi(:,1)+1.0_GP
      fi(:,2) = fi(:,2)-2.0_GP
      fi(:,n-1) = fi(:,n-1)+3.0_GP
      fi(:,n) = fi(:,n)-4.0_GP
      DO lm = lmsta,lmend
         IF (ell(lm).eq.0) fi(lm,:) = (5.0_GP,1.0_GP)
      END DO
      CALL v_imposebc(tor,fi,pol)

      err = 0.0_GP
      CALL derivr_mat(1,d1)
      err(1) = MAX(MAXVAL(ABS(tor(:,1))),MAXVAL(ABS(tor(:,n))))
      err(2) = MAX(MAXVAL(ABS(tor(:,2:n-1)-z0(:,2:n-1))), &
                   MAXVAL(ABS(tor(:,n+1:nz)-z0(:,n+1:nz))))
      DO lm = lmsta,lmend
         IF (ell(lm).eq.0) THEN
            err(5) = MAX(err(5),MAXVAL(ABS(fi(lm,:))),MAXVAL(ABS(pol(lm,:))))
            CYCLE
         ENDIF
         wm = MAXVAL(ABS(pol(lm,1:n)))
         err(3) = MAX(err(3),MAX(ABS(pol(lm,1)),ABS(pol(lm,n)))/wm)
         err(4) = MAX(err(4),MAX(ABS(SUM(d1(1,:)*pol(lm,1:n))), &
                  ABS(SUM(d1(n,:)*pol(lm,1:n))))*(r2-r1)/wm)
         err(6) = MAX(err(6),MAXVAL(ABS(fi(lm,3:n-2)-f0(lm,3:n-2)))/ &
                  MAXVAL(ABS(f0(lm,3:n-2))))
         err(8) = MAX(err(8),MAXVAL(ABS(pol(lm,1:n)-w0(lm,1:n)))/ &
                  MAXVAL(ABS(w0(lm,1:n))))
      END DO
! Aplicarla de nuevo no cambia nada
      CALL v_imposebc(tor,fi,pol2)
      dmax(1) = MAXVAL(ABS(pol2-pol))
      dmax(2) = MAXVAL(ABS(pol))
      CALL MPI_ALLREDUCE(dmax,gmax,2,GC_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      err(7) = gmax(1)/gmax(2)

      CALL MPI_REDUCE(err,emax,8,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         PRINT*, 'tareas:', nprocs, '  lmax =', lmax, '  nz =', nz
         PRINT*, 'tipos en las paredes:', bcp%bcint, bcp%bcext
         PRINT '(3X,A,ES10.2)', 'Z en las paredes                 :', emax(1)
         PRINT '(3X,A,ES10.2)', 'Z fuera de las paredes, cambio   :', emax(2)
         PRINT '(3X,A,ES10.2)', 'W en las paredes, relativo       :', emax(3)
         PRINT '(3X,A,ES10.2)', 'dW/dr en las paredes, relativo   :', emax(4)
         PRINT '(3X,A,ES10.2)', 'phi y W en l = 0                 :', emax(5)
         PRINT '(3X,A,ES10.2)', 'phi en el interior, cambio       :', emax(6)
         PRINT '(3X,A,ES10.2)', 'aplicarla dos veces, cambio de W :', emax(7)
         PRINT '(3X,A,ES10.2)', 'W contra W0 (truncado)           :', emax(8)
      ENDIF
      CALL terminar()
      CALL MPI_FINALIZE(ierr)

      CONTAINS

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

      END PROGRAM prueba_boundary
