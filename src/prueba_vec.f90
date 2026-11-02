!=================================================================
! Prueba de las transformadas vectoriales de sphp (sphp_vec.F90):
! análisis y síntesis de campos vectoriales conocidos, con sus
! coeficientes Q, S, T exactos en la convención de SHTns,
!    v = Q Y e_r + S r grad(Y) - T r x grad(Y),
!    v_r = Q,  v_th = dS/dth + (1/sin th) dT/dph,
!    v_ph = (1/sin th) dS/dph - dT/dth,
! donde Q, S y T son también los campos escalares que arman sus
! coeficientes. Como en prueba_sphp, el campo va multiplicado por
! k en el plano k y corre con cualquier cantidad de tareas hasta nz.
!=================================================================

      PROGRAM prueba_vec
      USE fprecision
      USE commtypes
      USE mpivars
      USE fcgram
      IMPLICIT NONE

! Grilla angular de la regla de 3/2 y nz planos periódicos
      INTEGER, PARAMETER :: lmax = 15, nz = 6
      INTEGER, PARAMETER :: nlat = 3*(lmax+1)/2, nphi = 2*nlat
      INTEGER, PARAMETER :: nlm = (lmax+1)*(lmax+2)/2
      INTEGER, PARAMETER :: ncasos = 8
      REAL(KIND=GP), PARAMETER :: pi = 3.14159265358979323846_GP

      TYPE(FCPLAN) :: plan
      TYPE(shtns_info), POINTER :: info
      REAL(C_DOUBLE), POINTER   :: ct(:),st(:)
! cth, sth y ph: cos(theta), sin(theta) y phi en la grilla
      REAL(KIND=GP), DIMENSION(nlat,nphi) :: cth,sth,ph,fr,ft,fp
      REAL(KIND=GP), ALLOCATABLE    :: vr(:,:,:),vt(:,:,:),vp(:,:,:)
      COMPLEX(KIND=GP), ALLOCATABLE :: q(:,:),s(:,:),t(:,:)
      COMPLEX(KIND=GP) :: aq(nlm),as(nlm),at(nlm)
      REAL(KIND=GP)    :: xr(nlm),xi(nlm),err(2),emax(2),c
      CHARACTER(len=36):: nombre
      INTEGER :: ic,j,k

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL range(1,nlm,nprocs,myrank,lmsta,lmend)
      ALLOCATE( vr(nlat,nphi,ksta:kend),vt(nlat,nphi,ksta:kend) )
      ALLOCATE( vp(nlat,nphi,ksta:kend) )
      ALLOCATE( q(lmsta:lmend,nz),s(lmsta:lmend,nz),t(lmsta:lmend,nz) )

      CALL fcgram_create_plan(plan,(/nlat,nphi,nz/),(/0,0,0/), &
                              (/0,0,0/),(/1.0_GP,1.0_GP,1.0_GP/), &
                              '.',FFTW_ESTIMATE,lmax)

! La grilla: los nodos de Gauss en theta los da SHTns
      CALL C_F_POINTER(plan%planrc%planrxy,info)
      CALL C_F_POINTER(info%ct,ct,(/nlat/))
      CALL C_F_POINTER(info%st,st,(/nlat/))
      DO j = 1,nphi
         cth(:,j) = ct
         sth(:,j) = st
         ph(:,j)  = 2*pi*(j-1)/nphi
      END DO

      IF (myrank.eq.0) THEN
         PRINT*, 'tareas:', nprocs
         PRINT*, 'caso, análisis max|QST-exactos|, ', &
                 'síntesis max|v-v_exacto|'
      ENDIF
      DO ic = 1,ncasos
         CALL caso(ic)
! Análisis: la directa destruye v, así que se arma cada vez
         DO k = ksta,kend
            vr(:,:,k) = k*fr
            vt(:,:,k) = k*ft
            vp(:,:,k) = k*fp
         END DO
         CALL fftp2d_vec_to_qst(plan,vr,vt,vp,q,s,t,MPI_COMM_WORLD)
         err(1) = 0.0_GP
         DO k = 1,nz
            err(1) = MAX(err(1),MAXVAL(ABS(q(:,k)-k*aq(lmsta:lmend))), &
                     MAXVAL(ABS(s(:,k)-k*as(lmsta:lmend))), &
                     MAXVAL(ABS(t(:,k)-k*at(lmsta:lmend))))
         END DO
! Síntesis desde los coeficientes exactos
         DO k = 1,nz
            q(:,k) = k*aq(lmsta:lmend)
            s(:,k) = k*as(lmsta:lmend)
            t(:,k) = k*at(lmsta:lmend)
         END DO
         CALL fftp2d_qst_to_vec(plan,q,s,t,vr,vt,vp,MPI_COMM_WORLD)
         err(2) = 0.0_GP
         DO k = ksta,kend
            err(2) = MAX(err(2),MAXVAL(ABS(vr(:,:,k)-k*fr)), &
                     MAXVAL(ABS(vt(:,:,k)-k*ft)), &
                     MAXVAL(ABS(vp(:,:,k)-k*fp)))
         END DO
         CALL MPI_REDUCE(err,emax,2,GC_REAL,MPI_MAX,0, &
                         MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) PRINT '(1X,A36,2ES10.2)', nombre,emax
      END DO

! Ida y vuelta de coeficientes al azar, los mismos en todas las
! tareas: los m = 0 reales, y S y T sin l = 0 (no tienen sentido)
      IF (myrank.eq.0) THEN
         CALL RANDOM_NUMBER(xr)
         CALL RANDOM_NUMBER(xi)
      ENDIF
      CALL MPI_BCAST(xr,nlm,GC_REAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_BCAST(xi,nlm,GC_REAL,0,MPI_COMM_WORLD,ierr)
      aq = CMPLX(xr-0.5_GP,xi-0.5_GP,KIND=GP)
      aq(1:lmax+1) = REAL(aq(1:lmax+1),KIND=GP)
      as = 2*aq
      at = -aq
      as(1) = 0.0_GP
      at(1) = 0.0_GP
      DO k = 1,nz
         q(:,k) = k*aq(lmsta:lmend)
         s(:,k) = k*as(lmsta:lmend)
         t(:,k) = k*at(lmsta:lmend)
      END DO
      CALL fftp2d_qst_to_vec(plan,q,s,t,vr,vt,vp,MPI_COMM_WORLD)
      CALL fftp2d_vec_to_qst(plan,vr,vt,vp,q,s,t,MPI_COMM_WORLD)
      err(1) = 0.0_GP
      DO k = 1,nz
         err(1) = MAX(err(1),MAXVAL(ABS(q(:,k)-k*aq(lmsta:lmend))), &
                  MAXVAL(ABS(s(:,k)-k*as(lmsta:lmend))), &
                  MAXVAL(ABS(t(:,k)-k*at(lmsta:lmend))))
      END DO
      CALL MPI_REDUCE(err,emax,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      nombre = 'ida y vuelta al azar'
      IF (myrank.eq.0) PRINT '(1X,A36,ES10.2)', nombre,emax(1)

      CALL fcgram_destroy_plan(plan)
      CALL MPI_FINALIZE(ierr)

      CONTAINS

!-----------------------------------------------------------------
! Arma las componentes fr, ft, fp del caso ic en la grilla y sus
! coeficientes exactos aq, as, at
      SUBROUTINE caso(ic)
      INTEGER, INTENT(IN) :: ic
      REAL(KIND=GP), DIMENSION(nlat,nphi) :: p0,p1,p2
      INTEGER :: k

      aq = 0.0_GP
      as = 0.0_GP
      at = 0.0_GP
      fr = 0.0_GP
      ft = 0.0_GP
      fp = 0.0_GP
      SELECT CASE (ic)
      CASE (1)
! Radial: Q = cos(th)
         nombre = 'radial, Q = cos(th)'
         fr = cth
         aq(lm(1,0)) = SQRT(4*pi/3)
      CASE (2)
! Esferoidal: S = cos(th)
         nombre = 'esferoidal, S = cos(th)'
         ft = -sth
         as(lm(1,0)) = SQRT(4*pi/3)
      CASE (3)
! Toroidal, la rotación rígida: T = cos(th)
         nombre = 'toroidal, T = cos(th)'
         fp = sth
         at(lm(1,0)) = SQRT(4*pi/3)
      CASE (4)
         nombre = 'esferoidal, S = sin(th) cos(ph)'
         ft = cth*COS(ph)
         fp = -SIN(ph)
         as(lm(1,1)) = SQRT(2*pi/3)
      CASE (5)
         nombre = 'toroidal, T = sin(th) cos(ph)'
         ft = -SIN(ph)
         fp = -cth*COS(ph)
         at(lm(1,1)) = SQRT(2*pi/3)
      CASE (6)
! P_lmax(cos(th)), con la recurrencia de Legendre
         nombre = 'radial, Q = P_lmax(cos(th))'
         p0 = 1.0_GP
         p1 = cth
         DO k = 1,lmax-1
            p2 = ((2*k+1)*cth*p1-k*p0)/(k+1)
            p0 = p1
            p1 = p2
         END DO
         fr = p1
         aq(lm(lmax,0)) = SQRT(4*pi/(2*lmax+1))
      CASE (7)
! Y_ll = c sin(th)^l exp(i l ph), con
! c = sqrt((2l+1)/(4 pi)) (2l-1)!!/sqrt((2l)!)
         nombre = 'toroidal, sin(th)^lmax cos(lmax ph)'
         ft = -lmax*sth**(lmax-1)*SIN(lmax*ph)
         fp = -lmax*sth**(lmax-1)*cth*COS(lmax*ph)
         c = SQRT((2*lmax+1)/(4*pi))
         DO k = 1,lmax
            c = c*SQRT((2*k-1.0_GP)/(2*k))
         END DO
         at(lm(lmax,lmax)) = 1/(2*c)
      CASE (8)
         nombre = 'suma de los casos 1, 4 y 5'
         fr = cth
         ft = cth*COS(ph)-SIN(ph)
         fp = -SIN(ph)-cth*COS(ph)
         aq(lm(1,0)) = SQRT(4*pi/3)
         as(lm(1,1)) = SQRT(2*pi/3)
         at(lm(1,1)) = SQRT(2*pi/3)
      END SELECT

      END SUBROUTINE caso

!-----------------------------------------------------------------
! Lugar del modo (l,m) en los arreglos de coeficientes
      INTEGER FUNCTION lm(l,m)
      INTEGER, INTENT(IN) :: l,m
      lm = shtns_lmidx(plan%planrc%planrxy,l,m)
      END FUNCTION lm

      END PROGRAM prueba_vec
