!=================================================================
! Prueba de sphp en REDONDITO, que crece con cada paso: análisis y
! síntesis de campos conocidos, con sus coeficientes exactos en los
! armónicos ortonormales sin la fase de Condon-Shortley. Como el
! campo es real y solo se guardan los m >= 0,
!    f = sum_l a_l0 Y_l0 + sum_(m>0) 2 Re(a_lm Y_lm).
! Las dos direcciones pasan por fftp2d_real_to_complex_xy y
! fftp2d_complex_to_real_xy, con el campo multiplicado por k en el
! plano k. Corre con cualquier cantidad de tareas hasta nz: cada
! una tiene sus planos en el espacio físico y sus lm en el
! espectral, y los errores son el máximo sobre todas.
!=================================================================

      PROGRAM prueba_sphp
      USE fprecision
      USE commtypes
      USE mpivars
      USE fcgram
      IMPLICIT NONE

! Grilla angular de la regla de 3/2 y nz planos periódicos
      INTEGER, PARAMETER :: lmax = 15, nz = 6
      INTEGER, PARAMETER :: nlat = 3*(lmax+1)/2, nphi = 2*nlat
      INTEGER, PARAMETER :: nlm = (lmax+1)*(lmax+2)/2
      INTEGER, PARAMETER :: ncasos = 9
      REAL(KIND=GP), PARAMETER :: pi = 3.14159265358979323846_GP

      TYPE(FCPLAN) :: plan
      TYPE(shtns_info), POINTER :: info
      REAL(C_DOUBLE), POINTER   :: ct(:),st(:)
! cth, sth y ph: cos(theta), sin(theta) y phi en la grilla
      REAL(KIND=GP), DIMENSION(nlat,nphi) :: cth,sth,ph,f
      REAL(KIND=GP), ALLOCATABLE    :: f3(:,:,:)
      COMPLEX(KIND=GP), ALLOCATABLE :: b3(:,:)
      REAL(KIND=GP)    :: xr(nlm),xi(nlm),err(2),emax(2)
      COMPLEX(KIND=GP) :: a(nlm)
      CHARACTER(len=36):: nombre
      INTEGER :: ic,j,k

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,nlat/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL range(1,nlm,nprocs,myrank,lmsta,lmend)
      ALLOCATE( f3(nlat,nphi,ksta:kend),b3(lmsta:lmend,nz) )

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
         PRINT*, 'caso, análisis max|b-a|, síntesis max|g-f|'
      ENDIF
      DO ic = 1,ncasos
         CALL caso(ic)
! Análisis: la directa destruye f3, así que se arma cada vez
         DO k = ksta,kend
            f3(:,:,k) = k*f
         END DO
         CALL fftp2d_real_to_complex_xy(plan,f3,b3,MPI_COMM_WORLD)
         err(1) = 0.0_GP
         DO k = 1,nz
            err(1) = MAX(err(1),MAXVAL(ABS(b3(:,k)-k*a(lmsta:lmend))))
         END DO
! Síntesis desde los coeficientes exactos
         DO k = 1,nz
            b3(:,k) = k*a(lmsta:lmend)
         END DO
         CALL fftp2d_complex_to_real_xy(plan,b3,f3,MPI_COMM_WORLD)
         err(2) = 0.0_GP
         DO k = ksta,kend
            err(2) = MAX(err(2),MAXVAL(ABS(f3(:,:,k)-k*f)))
         END DO
         CALL MPI_REDUCE(err,emax,2,GC_REAL,MPI_MAX,0, &
                         MPI_COMM_WORLD,ierr)
         IF (myrank.eq.0) PRINT '(1X,A36,2ES10.2)', nombre,emax
      END DO

! Ida y vuelta de coeficientes al azar (los m = 0, reales), los
! mismos en todas las tareas
      IF (myrank.eq.0) THEN
         CALL RANDOM_NUMBER(xr)
         CALL RANDOM_NUMBER(xi)
         a = CMPLX(xr-0.5_GP,xi-0.5_GP,KIND=GP)
         a(1:lmax+1) = REAL(a(1:lmax+1),KIND=GP)
      ENDIF
      CALL MPI_BCAST(a,nlm,GC_COMPLEX,0,MPI_COMM_WORLD,ierr)
      DO k = 1,nz
         b3(:,k) = k*a(lmsta:lmend)
      END DO
      CALL fftp2d_complex_to_real_xy(plan,b3,f3,MPI_COMM_WORLD)
      CALL fftp2d_real_to_complex_xy(plan,f3,b3,MPI_COMM_WORLD)
      err(1) = 0.0_GP
      DO k = 1,nz
         err(1) = MAX(err(1),MAXVAL(ABS(b3(:,k)-k*a(lmsta:lmend))))
      END DO
      CALL MPI_REDUCE(err,emax,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      nombre = 'ida y vuelta al azar'
      IF (myrank.eq.0) PRINT '(1X,A36,ES10.2)', nombre,emax(1)

      CALL fcgram_destroy_plan(plan)
      CALL MPI_FINALIZE(ierr)

      CONTAINS

!-----------------------------------------------------------------
! Arma el campo f del caso ic en la grilla y sus coeficientes
! exactos a
      SUBROUTINE caso(ic)
      INTEGER, INTENT(IN) :: ic
      REAL(KIND=GP), DIMENSION(nlat,nphi) :: p0,p1,p2
      REAL(KIND=GP) :: c
      INTEGER :: k

      a = 0.0_GP
      SELECT CASE (ic)
      CASE (1)
         nombre = '1'
         f = 1.0_GP
         a(lm(0,0)) = SQRT(4*pi)
      CASE (2)
         nombre = 'cos(th)'
         f = cth
         a(lm(1,0)) = SQRT(4*pi/3)
      CASE (3)
         nombre = 'sin(th) cos(ph)'
         f = sth*COS(ph)
         a(lm(1,1)) = SQRT(2*pi/3)
      CASE (4)
         nombre = 'sin(th) sin(ph)'
         f = sth*SIN(ph)
         a(lm(1,1)) = CMPLX(0.0_GP,-SQRT(2*pi/3),KIND=GP)
      CASE (5)
         nombre = '3 cos(th)^2 - 1'
         f = 3*cth**2-1
         a(lm(2,0)) = 4*SQRT(pi/5)
      CASE (6)
         nombre = 'sin(th)^2 cos(2 ph)'
         f = sth**2*COS(2*ph)
         a(lm(2,2)) = 2*SQRT(2*pi/15)
      CASE (7)
! P_lmax(cos(th)), con la recurrencia de Legendre
         nombre = 'P_lmax(cos(th))'
         p0 = 1.0_GP
         p1 = cth
         DO k = 1,lmax-1
            p2 = ((2*k+1)*cth*p1-k*p0)/(k+1)
            p0 = p1
            p1 = p2
         END DO
         f = p1
         a(lm(lmax,0)) = SQRT(4*pi/(2*lmax+1))
      CASE (8)
! Y_ll = c sin(th)^l exp(i l ph), con
! c = sqrt((2l+1)/(4 pi)) (2l-1)!!/sqrt((2l)!)
         nombre = 'sin(th)^lmax cos(lmax ph)'
         f = sth**lmax*COS(lmax*ph)
         c = SQRT((2*lmax+1)/(4*pi))
         DO k = 1,lmax
            c = c*SQRT((2*k-1.0_GP)/(2*k))
         END DO
         a(lm(lmax,lmax)) = 1/(2*c)
      CASE (9)
         nombre = 'suma: 2(1) - 3(2) + (4) + (6)/2'
         f = 2-3*cth+sth*SIN(ph)+0.5_GP*sth**2*COS(2*ph)
         a(lm(0,0)) = 2*SQRT(4*pi)
         a(lm(1,0)) = -3*SQRT(4*pi/3)
         a(lm(1,1)) = CMPLX(0.0_GP,-SQRT(2*pi/3),KIND=GP)
         a(lm(2,2)) = SQRT(2*pi/15)
      END SELECT

      END SUBROUTINE caso

!-----------------------------------------------------------------
! Lugar del modo (l,m) en los arreglos de coeficientes
      INTEGER FUNCTION lm(l,m)
      INTEGER, INTENT(IN) :: l,m
      lm = shtns_lmidx(plan%planrc%planrxy,l,m)
      END FUNCTION lm

      END PROGRAM prueba_sphp
