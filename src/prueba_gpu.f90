!=================================================================
! Prueba de sphp-gpu en REDONDITO: compara el camino en la GPU con
! el de la CPU, ya validado en prueba_sphp, prueba_vec, prueba_z y
! prueba_3d, para todas las transformadas: las 2D escalares y
! vectoriales, la de z con FC-Gram (C = 25, d = 5) y las 3D. Los
! coeficientes angulares son al azar (los m = 0, reales; S y T sin
! l = 0), los mismos en todas las tareas, multiplicados por k en el
! plano k; los datos de z son al azar en cada tarea. Una GPU por
! tarea (device_init) y cualquier cantidad de tareas hasta nz; los
! errores son el máximo sobre todas. Con GHOST_GPU_EXCHANGE=staged
! el intercambio pasa por la CPU. Cada arreglo de la GPU tiene su
! propia copia en la GPU (gdev_alloc), como los campos de SPECTER.
!=================================================================

      PROGRAM prueba_gpu
      USE fprecision
      USE commtypes
      USE mpivars
      USE fcgram
      USE gdevice
      IMPLICIT NONE

! Grilla angular de la regla de 3/2 y nz planos con FC-Gram en z
      INTEGER, PARAMETER :: lmax = 15, nz = 64, C = 25, od = 5
      INTEGER, PARAMETER :: nlat = 3*(lmax+1)/2, nphi = 2*nlat
      INTEGER, PARAMETER :: nlm = (lmax+1)*(lmax+2)/2
      CHARACTER(len=*), PARAMETER :: tdir = &
               '/home/mramdan/Doctorado/REDONDO/SPECTER/tables'

      TYPE(FCPLAN) :: plan
! Campos y coeficientes de la CPU (f, b, una columna por componente)
! y de la GPU (vr, vt, vp y q, s, t; los escalares van en vr y q)
      REAL(KIND=GP), ALLOCATABLE    :: f(:,:,:,:)
      COMPLEX(KIND=GP), ALLOCATABLE :: b(:,:,:)
      REAL(KIND=GP), ALLOCATABLE, TARGET    :: vr(:,:,:),vt(:,:,:),vp(:,:,:)
      COMPLEX(KIND=GP), ALLOCATABLE, TARGET :: q(:,:),s(:,:),t(:,:)
      REAL(KIND=GP), ALLOCATABLE    :: zr(:,:),zi(:,:)
      REAL(KIND=GP)     :: xr(nlm),xi(nlm),err(10),emax(10),dz
      COMPLEX(KIND=GP)  :: a(nlm,3)
      CHARACTER(len=36) :: nombre(10)
      INTEGER(C_SIZE_T) :: nbf,nbb
      INTEGER :: ic,k

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
      CALL range(1,nlat/2+1,nprocs,myrank,ista,iend)
      CALL range(1,nz,nprocs,myrank,ksta,kend)
      CALL range(1,nlm,nprocs,myrank,lmsta,lmend)
      CALL device_init(myrank)
      ALLOCATE( f(nlat,nphi,ksta:kend,3),b(lmsta:lmend,nz,3) )
      ALLOCATE( vr(nlat,nphi,ksta:kend),vt(nlat,nphi,ksta:kend) )
      ALLOCATE( vp(nlat,nphi,ksta:kend) )
      ALLOCATE( q(lmsta:lmend,nz),s(lmsta:lmend,nz),t(lmsta:lmend,nz) )
      ALLOCATE( zr(lmsta:lmend,nz),zi(lmsta:lmend,nz) )
      nbf = INT(SIZE(vr),C_SIZE_T)*INT(STORAGE_SIZE(vr)/8,C_SIZE_T)
      nbb = INT(SIZE(q),C_SIZE_T)*INT(STORAGE_SIZE(q)/8,C_SIZE_T)

      dz = 1.0_GP/(nz-C-1)
      CALL fcgram_create_plan(plan,(/nlat,nphi,nz/),(/0,0,C/), &
                              (/0,0,od/),(/1.0_GP,1.0_GP,dz/), &
                              tdir,FFTW_ESTIMATE,lmax)
      CALL gdev_alloc(C_LOC(vr),nbf)
      CALL gdev_alloc(C_LOC(vt),nbf)
      CALL gdev_alloc(C_LOC(vp),nbf)
      CALL gdev_alloc(C_LOC(q),nbb)
      CALL gdev_alloc(C_LOC(s),nbb)
      CALL gdev_alloc(C_LOC(t),nbb)

      IF (myrank.eq.0) THEN
         DO ic = 1,3
            CALL RANDOM_NUMBER(xr)
            CALL RANDOM_NUMBER(xi)
            a(:,ic) = CMPLX(xr-0.5_GP,xi-0.5_GP,KIND=GP)
            a(1:lmax+1,ic) = REAL(a(1:lmax+1,ic),KIND=GP)
         END DO
         a(1,2:3) = 0.0_GP
      ENDIF
      CALL MPI_BCAST(a,3*nlm,GC_COMPLEX,0,MPI_COMM_WORLD,ierr)

! 2D escalar: síntesis en la CPU y en la GPU, desde los mismos
! coeficientes
      nombre(1) = 'escalar, síntesis max|GPU-CPU|'
      CALL coeficientes()
      q = b(:,:,1)
      CALL gdev_update_to(C_LOC(q),nbb)
      CALL fftp2d_complex_to_real_xy(plan,b(:,:,1),f(:,:,:,1),MPI_COMM_WORLD)
      gdev_active = .TRUE.
      CALL fftp2d_complex_to_real_xy(plan,q,vr,MPI_COMM_WORLD)
      gdev_active = .FALSE.
      CALL gdev_update_from(C_LOC(vr),nbf)
      err(1) = MAXVAL(ABS(vr-f(:,:,:,1)))

! 2D escalar: análisis en la CPU y en la GPU, del mismo campo (la
! directa destruye su entrada, así que la de la GPU se carga antes)
      nombre(2) = 'escalar, análisis max|GPU-CPU|'
      vr = f(:,:,:,1)
      CALL gdev_update_to(C_LOC(vr),nbf)
      CALL fftp2d_real_to_complex_xy(plan,f(:,:,:,1),b(:,:,1),MPI_COMM_WORLD)
      gdev_active = .TRUE.
      CALL fftp2d_real_to_complex_xy(plan,vr,q,MPI_COMM_WORLD)
      gdev_active = .FALSE.
      CALL gdev_update_from(C_LOC(q),nbb)
      err(2) = MAXVAL(ABS(q-b(:,:,1)))

! 2D escalar: ida y vuelta en la GPU, contra los coeficientes de
! partida
      nombre(3) = 'escalar, ida y vuelta en la GPU'
      err(3) = 0.0_GP
      DO k = 1,nz
         err(3) = MAX(err(3),MAXVAL(ABS(q(:,k)-k*a(lmsta:lmend,1))))
      END DO

! 2D vectorial: lo mismo, de (Q,S,T) a (v_r,v_theta,v_phi) y al
! revés
      nombre(4) = 'vectorial, síntesis max|GPU-CPU|'
      CALL coeficientes()
      q = b(:,:,1)
      s = b(:,:,2)
      t = b(:,:,3)
      CALL gdev_update_to(C_LOC(q),nbb)
      CALL gdev_update_to(C_LOC(s),nbb)
      CALL gdev_update_to(C_LOC(t),nbb)
      CALL fftp2d_qst_to_vec(plan,b(:,:,1),b(:,:,2),b(:,:,3), &
              f(:,:,:,1),f(:,:,:,2),f(:,:,:,3),MPI_COMM_WORLD)
      gdev_active = .TRUE.
      CALL fftp2d_qst_to_vec(plan,q,s,t,vr,vt,vp,MPI_COMM_WORLD)
      gdev_active = .FALSE.
      CALL gdev_update_from(C_LOC(vr),nbf)
      CALL gdev_update_from(C_LOC(vt),nbf)
      CALL gdev_update_from(C_LOC(vp),nbf)
      err(4) = MAX(MAXVAL(ABS(vr-f(:,:,:,1))),MAXVAL(ABS(vt-f(:,:,:,2))), &
                   MAXVAL(ABS(vp-f(:,:,:,3))))

      nombre(5) = 'vectorial, análisis max|GPU-CPU|'
      vr = f(:,:,:,1)
      vt = f(:,:,:,2)
      vp = f(:,:,:,3)
      CALL gdev_update_to(C_LOC(vr),nbf)
      CALL gdev_update_to(C_LOC(vt),nbf)
      CALL gdev_update_to(C_LOC(vp),nbf)
      CALL fftp2d_vec_to_qst(plan,f(:,:,:,1),f(:,:,:,2),f(:,:,:,3), &
              b(:,:,1),b(:,:,2),b(:,:,3),MPI_COMM_WORLD)
      gdev_active = .TRUE.
      CALL fftp2d_vec_to_qst(plan,vr,vt,vp,q,s,t,MPI_COMM_WORLD)
      gdev_active = .FALSE.
      CALL gdev_update_from(C_LOC(q),nbb)
      CALL gdev_update_from(C_LOC(s),nbb)
      CALL gdev_update_from(C_LOC(t),nbb)
      err(5) = MAX(MAXVAL(ABS(q-b(:,:,1))),MAXVAL(ABS(s-b(:,:,2))), &
                   MAXVAL(ABS(t-b(:,:,3))))

      nombre(6) = 'vectorial, ida y vuelta en la GPU'
      err(6) = 0.0_GP
      DO k = 1,nz
         err(6) = MAX(err(6),MAXVAL(ABS(q(:,k)-k*a(lmsta:lmend,1))), &
                  MAXVAL(ABS(s(:,k)-k*a(lmsta:lmend,2))),            &
                  MAXVAL(ABS(t(:,k)-k*a(lmsta:lmend,3))))
      END DO

! z: directa (continuación de FC-Gram, FFT y división por nz) e
! inversa, de los mismos datos al azar en la CPU y en la GPU
      nombre(7) = 'z, directa max|GPU-CPU|'
      CALL RANDOM_NUMBER(zr)
      CALL RANDOM_NUMBER(zi)
      b(:,:,1) = CMPLX(zr-0.5_GP,zi-0.5_GP,KIND=GP)
      q = b(:,:,1)
      CALL gdev_update_to(C_LOC(q),nbb)
      CALL fftp1d_real_to_complex_z(plan,b(:,:,1),MPI_COMM_WORLD)
      gdev_active = .TRUE.
      CALL fftp1d_real_to_complex_z(plan,q,MPI_COMM_WORLD)
      gdev_active = .FALSE.
      CALL gdev_update_from(C_LOC(q),nbb)
      err(7) = MAXVAL(ABS(q-b(:,:,1)))

      nombre(8) = 'z, inversa max|GPU-CPU|'
      CALL RANDOM_NUMBER(zr)
      CALL RANDOM_NUMBER(zi)
      b(:,:,1) = CMPLX(zr-0.5_GP,zi-0.5_GP,KIND=GP)
      q = b(:,:,1)
      CALL gdev_update_to(C_LOC(q),nbb)
      CALL fftp1d_complex_to_real_z(plan,b(:,:,1),MPI_COMM_WORLD)
      gdev_active = .TRUE.
      CALL fftp1d_complex_to_real_z(plan,q,MPI_COMM_WORLD)
      gdev_active = .FALSE.
      CALL gdev_update_from(C_LOC(q),nbb)
      err(8) = MAXVAL(ABS(q-b(:,:,1)))

! 3D: inversa desde los mismos coeficientes y directa del mismo
! campo (las dos destruyen su entrada)
      nombre(9) = '3D, inversa max|GPU-CPU|'
      CALL coeficientes()
      q = b(:,:,1)
      CALL gdev_update_to(C_LOC(q),nbb)
      CALL fftp3d_complex_to_real(plan,b(:,:,1),f(:,:,:,1),MPI_COMM_WORLD)
      gdev_active = .TRUE.
      CALL fftp3d_complex_to_real(plan,q,vr,MPI_COMM_WORLD)
      gdev_active = .FALSE.
      CALL gdev_update_from(C_LOC(vr),nbf)
      err(9) = MAXVAL(ABS(vr-f(:,:,:,1)))

      nombre(10) = '3D, directa max|GPU-CPU|'
      vr = f(:,:,:,1)
      CALL gdev_update_to(C_LOC(vr),nbf)
      CALL fftp3d_real_to_complex(plan,f(:,:,:,1),b(:,:,1),MPI_COMM_WORLD)
      gdev_active = .TRUE.
      CALL fftp3d_real_to_complex(plan,vr,q,MPI_COMM_WORLD)
      gdev_active = .FALSE.
      CALL gdev_update_from(C_LOC(q),nbb)
      err(10) = MAXVAL(ABS(q-b(:,:,1)))

      CALL MPI_REDUCE(err,emax,10,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         PRINT*, 'tareas:', nprocs
         DO k = 1,10
            PRINT '(1X,A36,ES10.2)', nombre(k),emax(k)
         END DO
      ENDIF

      CALL gdev_free(C_LOC(vr))
      CALL gdev_free(C_LOC(vt))
      CALL gdev_free(C_LOC(vp))
      CALL gdev_free(C_LOC(q))
      CALL gdev_free(C_LOC(s))
      CALL gdev_free(C_LOC(t))
      CALL fcgram_destroy_plan(plan)
      CALL MPI_FINALIZE(ierr)

      CONTAINS

!-----------------------------------------------------------------
! Los coeficientes de la CPU, b(:,k,ic) = k a(:,ic) en mis lm
      SUBROUTINE coeficientes()
      INTEGER :: ic,k
      DO ic = 1,3
         DO k = 1,nz
            b(:,k,ic) = k*a(lmsta:lmend,ic)
         END DO
      END DO
      END SUBROUTINE coeficientes

      END PROGRAM prueba_gpu
