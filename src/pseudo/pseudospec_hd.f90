!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines to compute spatial derivatives and nonlinear
! terms in Navier-Stokes equations in 3D using a 
! pseudospectral-FCGram method. You should use the FFTPLANS
! and MPIVARS modules (see the file 'fftp_mod.f90') in each
! program that calls any of the subroutines in this file.
!
! NOTATION: index 'i' is 'x'
!           index 'j' is 'y'
!           index 'k' is 'z'
!
! 2019 Mauro Fontana & Pablo D. Mininni.
!      Department of Physics,
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mfontana@df.uba.ar
!
! 17 Oct 2026: Enstrophy in energy (kin=0) sums the three
!              components of the curl (J. Cullen & G. Marchetti)
!=================================================================

!******************************************************************
!  TODO: benchmark populating with zeros after nz-Cz and summing up
! to nz to see if it helps vectorization. (Energy and related
! subroutines)
!******************************************************************

!*****************************************************************
      SUBROUTINE derivr(a,b,p,esp)
!-----------------------------------------------------------------
!
! Derivada radial de orden p de la matriz 'a' (z es el radio), en
! el espacio que indica esp: llama a derivr_r en (lm,r) o a
! derivr_k en (lm,kz). Las derivadas angulares salen de SHTns, en
! la síntesis.
!
! Parameters
!     a   : input matrix
!     b   : at the output contains the derivative d^p a/dz^p, en el
!           mismo espacio que a
!     p   : orden de la derivada
!     esp : espacio de a y b, 'r' para (lm,r) o 'k' para (lm,kz)
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(lmsta:lmend,nz) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz) :: b
      INTEGER, INTENT(IN)          :: p
      CHARACTER(len=1), INTENT(IN) :: esp

      SELECT CASE (esp)
      CASE ('r')
         CALL derivr_r(a,b,p)
      CASE ('k')
         CALL derivr_k(a,b,p)
      CASE DEFAULT
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Unknown space ", &
             esp, " in derivr. Aborting..."
         CALL MPI_FINALIZE(ierr)
         STOP
      END SELECT

      RETURN
      END SUBROUTINE derivr

!*****************************************************************
      SUBROUTINE derivr_k(a,b,p)
!-----------------------------------------------------------------
!
! Derivada radial de orden p de la matriz 'a' en (lm,kz):
! multiplica por (i kz)**p (z es el radio).
!
! Parameters
!     a  : input matrix in (lm,kz)
!     b  : at the output contains d^p a/dz^p, in (lm,kz)
!     p  : orden de la derivada
!
      USE kes
      USE var
      USE grid
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(lmsta:lmend,nz) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz) :: b
      INTEGER, INTENT(IN) :: p
      COMPLEX(KIND=GP) :: ip
      INTEGER :: lm,k

! i**p se calcula una vez, fuera del lazo
      ip = im**p
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active)
#else
!$omp parallel do private (lm)
#endif
      DO k = 1,nz
         DO lm = lmsta,lmend
            b(lm,k) = ip*kz(k)**p*a(lm,k)
         END DO
      END DO

      RETURN
      END SUBROUTINE derivr_k

!*****************************************************************
      SUBROUTINE derivr_r(a,b,p)
!-----------------------------------------------------------------
!
! Derivada radial de orden p de la matriz 'a' en (lm,r): ida a
! (lm,kz) con FC-Gram, derivr_k y vuelta. La ida se hace sobre una
! copia, así la entrada no cambia.
!
! Parameters
!     a  : input matrix in (lm,r)
!     b  : at the output contains d^p a/dz^p, in (lm,r)
!     p  : orden de la derivada
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      USE class_GWorkspace3D, ONLY: gws
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(lmsta:lmend,nz) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz) :: b
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:) :: c1
      INTEGER, INTENT(IN) :: p
      LOGICAL :: bret

      CALL gws%get_complex_tmp(c1,bret)
      CALL copy(a,c1)
      CALL fftp1d_real_to_complex_z(planfc,c1,MPI_COMM_WORLD)
      CALL derivr_k(c1,b,p)
      CALL fftp1d_complex_to_real_z(planfc,b,MPI_COMM_WORLD)
      CALL gws%free_complex_tmp(c1)

      RETURN
      END SUBROUTINE derivr_r

!*****************************************************************
      SUBROUTINE laplang(a,b)
!-----------------------------------------------------------------
!
! Laplaciano angular de la matriz 'a': el laplaciano sobre la
! esfera de radio uno, que multiplica cada modo por -l(l+1). Vale
! en (lm,r) y en (lm,kz), porque no depende de r. El laplaciano
! horizontal es este dividido por r**2.
!
! Parameters
!     a: input matrix
!     b: at the output contains the angular Laplacian
!
      USE kes
      USE grid
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(lmsta:lmend,nz) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz) :: b
      INTEGER :: lm,k

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active)
#else
!$omp parallel do private (lm)
#endif
      DO k = 1,nz
         DO lm = lmsta,lmend
            b(lm,k) = -ell(lm)*(ell(lm)+1)*a(lm,k)
         END DO
      END DO

      RETURN
      END SUBROUTINE laplang

!*****************************************************************
      SUBROUTINE multr(a,b,p)
!-----------------------------------------------------------------
!
! Multiplica la matriz 'a' en (lm,r) por r**p, con p entero
! (z es el radio).
!
! Parameters
!     a  : input matrix
!     b  : at the output contains r**p a
!     p  : potencia de r
!
      USE grid
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(lmsta:lmend,nz) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz) :: b
      INTEGER, INTENT(IN) :: p
      INTEGER :: lm,k

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active)
#else
!$omp parallel do private (lm)
#endif
      DO k = 1,nz
         DO lm = lmsta,lmend
            b(lm,k) = z(k)**p*a(lm,k)
         END DO
      END DO

      RETURN
      END SUBROUTINE multr

!*****************************************************************
      SUBROUTINE prodvec(ar,at,ap,br,bt,bp,cr,ct,cp)
!-----------------------------------------------------------------
!
! Producto vectorial C = A x B punto a punto en el espacio físico
! (theta,phi,r), con las componentes en (r,theta,phi), en el orden
! en que salen de fftp2d_qst_to_vec. No hace falta normalizar: las
! transformadas de sphp ya salen normalizadas.
!
! Parameters
!     ar,at,ap : componentes r, theta y phi de A
!     br,bt,bp : componentes r, theta y phi de B
!     cr,ct,cp : componentes r, theta y phi de A x B [output]
!
      USE fprecision
      USE mpivars
      USE grid
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      REAL(KIND=GP), INTENT (IN), DIMENSION(nlat,nphi,ksta:kend) :: ar,at,ap
      REAL(KIND=GP), INTENT (IN), DIMENSION(nlat,nphi,ksta:kend) :: br,bt,bp
      REAL(KIND=GP), INTENT(OUT), DIMENSION(nlat,nphi,ksta:kend) :: cr,ct,cp
      INTEGER :: i,j,k

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,pkend
         DO j = 1,nphi
            DO i = 1,nlat
               cr(i,j,k) = at(i,j,k)*bp(i,j,k)-ap(i,j,k)*bt(i,j,k)
               ct(i,j,k) = ap(i,j,k)*br(i,j,k)-ar(i,j,k)*bp(i,j,k)
               cp(i,j,k) = ar(i,j,k)*bt(i,j,k)-at(i,j,k)*br(i,j,k)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE prodvec

!*****************************************************************
      SUBROUTINE dlr(a,b)
!-----------------------------------------------------------------
!
! Operador D_l = d^2/dr^2 - l(l+1)/r**2 de la matriz 'a' en (lm,r)
! (z es el radio): la parte radial con derivr y la angular con
! laplang y multr, sumadas con saxpby_c. Solo vale en (lm,r),
! porque 1/r**2 no es diagonal en kz. Usa dos temporarios
! complejos del pool (derivr usa uno, antes de pedir estos).
!
! Parameters
!     a  : input matrix in (lm,r)
!     b  : at the output contains D_l a, in (lm,r)
!
      USE fprecision
      USE mpivars
      USE grid
      USE class_GWorkspace3D, ONLY: gws
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(lmsta:lmend,nz) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz) :: b
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:) :: c1,c2
      LOGICAL :: bret

      CALL derivr(a,b,2,'r')
      CALL gws%get_complex_tmp(c1,bret)
      CALL gws%get_complex_tmp(c2,bret)
      CALL laplang(a,c1)
      CALL multr(c1,c2,-2)
      CALL saxpby_c(b,b,1.0_GP,c2,1.0_GP)
      CALL gws%free_complex_tmp(c1)
      CALL gws%free_complex_tmp(c2)

      RETURN
      END SUBROUTINE dlr

!*****************************************************************
      SUBROUTINE curlk(a,b,c,dir)
!-----------------------------------------------------------------
!
! Computes the curl of the vector field A in Fourier
! space. The needed components of the field A are
! given by the matrixes a and b, and the order must
! follow the right hand convention [(x,y), (x,z) or (y,z)].
!
! Parameters
!     a  : input matrix
!     b  : input matrix
!     c  : at the output contains curl(A)_dir
!     dir: =1 computes the x-component
!          =2 computes the y-component
!          =3 computes the z-component
!
      USE fprecision
      USE kes
      USE var
      USE grid
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: c
      INTEGER, INTENT(IN) :: dir
      INTEGER             :: i,j,k

!
! Computes the x-component: c = d_y(b) - d_z(a)
!
      IF (dir.eq.1) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  c(k,j,i) = im*(ky(j)*b(k,j,i)-kz(k)*a(k,j,i))
               END DO
            END DO
         END DO
!
! Computes the y-component: c = d_z(a) - d_x(b)
!
      ELSE IF (dir.eq.2) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  c(k,j,i) = im*(kz(k)*a(k,j,i)-kx(i)*b(k,j,i))
               END DO
            END DO
         END DO
!
! Computes the z-component: c = d_x(b) - d_y(a)
!
      ELSE
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (k)
#endif
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz
                  c(k,j,i) = im*(kx(i)*b(k,j,i)-ky(j)*a(k,j,i))
               END DO
            END DO
         END DO
      ENDIF

      RETURN
      END SUBROUTINE curlk

!*****************************************************************
      SUBROUTINE gradre(a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Three-dimensional inner product A.grad(A) in
! real space. The components of the field A are
! given by the matrixes a, b and c.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: product (A.grad)A_x in Fourier space [output]
!     e: product (A.grad)A_y in Fourier space [output]
!     f: product (A.grad)A_z in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      USE class_GWorkspace3D, ONLY: gws
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: d,e,f
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3,c4
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:)    :: r1,r2,r3,r4
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:)    :: rx,ry,rz
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k,dir
      LOGICAL :: bret

      CALL gws%get_complex_tmp(c1,bret)
      CALL gws%get_complex_tmp(c2,bret)
      CALL gws%get_complex_tmp(c3,bret)
      CALL gws%get_complex_tmp(c4,bret)
      CALL gws%get_real_tmp(r1,bret)
      CALL gws%get_real_tmp(r2,bret)
      CALL gws%get_real_tmp(r3,bret)
      CALL gws%get_real_tmp(r4,bret)
      CALL gws%get_real_tmp(rx,bret)
      CALL gws%get_real_tmp(ry,bret)
      CALL gws%get_real_tmp(rz,bret)
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
!
! Computes (A_x.dx)A_dir, (A_y.dy)A_dir and (A_z.dz)A_dir, and
! accumulates them in rx, ry, rz (only in the physical domain,
! up to pkend)
!
      DO dir = 1,3
         IF (dir.eq.1) THEN
            CALL copy(a,c1)
         ELSE IF (dir.eq.2) THEN
            CALL copy(b,c1)
         ELSE
            CALL copy(c,c1)
         ENDIF
         CALL derivk(a,c2,dir)
         CALL derivk(b,c3,dir)
         CALL derivk(c,c4,dir)
         CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(planfc,c3,r3,MPI_COMM_WORLD)
         CALL fftp3d_complex_to_real(planfc,c4,r4,MPI_COMM_WORLD)
         IF (dir.eq.1) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
            DO k = ksta,pkend
               DO j = 1,ny
                  DO i = 1,nx
                     rx(i,j,k) = r1(i,j,k)*r2(i,j,k)
                     ry(i,j,k) = r1(i,j,k)*r3(i,j,k)
                     rz(i,j,k) = r1(i,j,k)*r4(i,j,k)
                  END DO
               END DO
            END DO
         ELSE IF (dir.eq.2) THEN
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
            DO k = ksta,pkend
               DO j = 1,ny
                  DO i = 1,nx
                     rx(i,j,k) = rx(i,j,k)+r1(i,j,k)*r2(i,j,k)
                     ry(i,j,k) = ry(i,j,k)+r1(i,j,k)*r3(i,j,k)
                     rz(i,j,k) = rz(i,j,k)+r1(i,j,k)*r4(i,j,k)
                  END DO
               END DO
            END DO
         ELSE
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
            DO k = ksta,pkend
               DO j = 1,ny
                  DO i = 1,nx
                     rx(i,j,k) = (rx(i,j,k)+r1(i,j,k)*r2(i,j,k))*tmp
                     ry(i,j,k) = (ry(i,j,k)+r1(i,j,k)*r3(i,j,k))*tmp
                     rz(i,j,k) = (rz(i,j,k)+r1(i,j,k)*r4(i,j,k))*tmp
                  END DO
               END DO
            END DO
         ENDIF
      END DO

      CALL fftp3d_real_to_complex(planfc,rx,d,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,ry,e,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,rz,f,MPI_COMM_WORLD)
      CALL gws%free_complex_tmp(c1)
      CALL gws%free_complex_tmp(c2)
      CALL gws%free_complex_tmp(c3)
      CALL gws%free_complex_tmp(c4)
      CALL gws%free_real_tmp(r1)
      CALL gws%free_real_tmp(r2)
      CALL gws%free_real_tmp(r3)
      CALL gws%free_real_tmp(r4)
      CALL gws%free_real_tmp(rx)
      CALL gws%free_real_tmp(ry)
      CALL gws%free_real_tmp(rz)

      RETURN
      END SUBROUTINE gradre

!*****************************************************************
      SUBROUTINE prodre(a,b,c,d,e,f)
!-----------------------------------------------------------------
!
! Three-dimensional cross product curl(A)xA in
! real space. The components of the field A are
! given by the matrixes a, b and c
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : product [curl(A)xA]_x in Fourier space [output]
!     e  : product [curl(A)xA]_y in Fourier space [output]
!     f  : product [curl(A)xA]_z in Fourier space [output]
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      USE class_GWorkspace3D, ONLY: gws
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(nz,ny,ista:iend) :: d,e,f
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:) :: r1,r2,r3,r4,r5,r6,r7
      REAL(KIND=GP)    :: tmp
      INTEGER :: i,j,k
      LOGICAL :: bret

      CALL gws%get_real_tmp(r1,bret)
      CALL gws%get_real_tmp(r2,bret)
      CALL gws%get_real_tmp(r3,bret)
      CALL gws%get_real_tmp(r4,bret)
      CALL gws%get_real_tmp(r5,bret)
      CALL gws%get_real_tmp(r6,bret)
      CALL gws%get_real_tmp(r7,bret)
!
! Computes curl(A)
!
      CALL curlk(b,c,d,1)
      CALL curlk(a,c,e,2)
      CALL curlk(a,b,f,3)
      CALL fftp3d_complex_to_real(planfc,d,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,e,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,f,r3,MPI_COMM_WORLD)
!
! Computes A
!
      CALL copy(a,d)
      CALL copy(b,e)
      CALL copy(c,f)
      CALL fftp3d_complex_to_real(planfc,d,r4,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,e,r5,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,f,r6,MPI_COMM_WORLD)
!
! Computes curl(A)xA
!
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2
#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,pkend
         DO j = 1,ny
            DO i = 1,nx
               r7(i,j,k) = (r2(i,j,k)*r6(i,j,k)-r5(i,j,k)*r3(i,j,k))*tmp
               r3(i,j,k) = (r3(i,j,k)*r4(i,j,k)-r6(i,j,k)*r1(i,j,k))*tmp
               r1(i,j,k) = (r1(i,j,k)*r5(i,j,k)-r4(i,j,k)*r2(i,j,k))*tmp
            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planfc,r7,d,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,r3,e,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planfc,r1,f,MPI_COMM_WORLD)
      CALL gws%free_real_tmp(r1)
      CALL gws%free_real_tmp(r2)
      CALL gws%free_real_tmp(r3)
      CALL gws%free_real_tmp(r4)
      CALL gws%free_real_tmp(r5)
      CALL gws%free_real_tmp(r6)
      CALL gws%free_real_tmp(r7)

      RETURN
      END SUBROUTINE prodre

!*****************************************************************
      SUBROUTINE energy(a,b,c,d,kin)
!-----------------------------------------------------------------
!
! Computes the mean energy of a vector field utilizing Parseval's
! theorem to compute sum(u^2) in (kx,ky,z) domain.
! The output is only valid in the first node.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : at the output contains the energy
!     kin: =0 computes the energy of the curl of the field
!          =1 computes the energy of the field
!          =2 computes the energy of the curl^2 of the field
!
      USE fprecision
      USE commtypes
      USE grid
      USE fft
      USE mpivars
      USE class_GWorkspace3D, ONLY: gws
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      INTEGER, INTENT(IN)           :: kin
      DOUBLE PRECISION, INTENT(OUT) :: d

      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: C1,C2,C3,C4
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: V1,V2,V3

      DOUBLE PRECISION  :: dloc
      REAL(KIND=GP)     :: tmp
      INTEGER           :: i,j,k
      LOGICAL           :: bret

      CALL gws%get_complex_htmp(C1,bret)
      CALL gws%get_complex_htmp(C2,bret)
      CALL gws%get_complex_htmp(C3,bret)
      CALL gws%get_complex_htmp(C4,bret)
      dloc = 0.0D0
      ! One nx*ny from parsevals, an nx*ny*(nz-Cz) from averaging over
      ! physical domain and a nz**2 due to fft normalization.
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2/ &
            real((nz-Cz),kind=GP)

      ! Computes the energy of field (a,b,c)
      IF (kin .eq. 1) THEN
         CALL copy(a,C1)
         CALL copy(b,C2)
         CALL copy(c,C3)
         V1 => C1; V2 => C2; V3 => C3

      ! Compute the energy of the curl of the field (a,b,c)
      ELSE IF (kin .eq. 0) THEN
         CALL curlk(b,c,C1,1)
         CALL curlk(a,c,C2,2)
         CALL curlk(a,b,C3,3)
         V1 => C1; V2 => C2; V3 => C3

      ! Compute the energy of the curl^2 of the field (a,b,c)
      ELSE IF (kin .eq. 2) THEN
         CALL curlk(b,c,C1,1)
         CALL curlk(a,c,C2,2)
         CALL curlk(a,b,C3,3)

         CALL curlk(C2,C3,C4,1)
         CALL curlk(C1,C3,C3,2)
         CALL curlk(C1,C2,C1,3)
         V1 => C4; V2 => C3; V3 => C1
      ENDIF
      CALL fftp1d_complex_to_real_z(planfc,V1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,V2,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,V3,MPI_COMM_WORLD)

      ! Compute mean value
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               dloc = dloc + (real(V1(k,j,1),kind=GP)**2+aimag(V1(k,j,1))**2+ &
                              real(V2(k,j,1),kind=GP)**2+aimag(V2(k,j,1))**2+ &
                              real(V3(k,j,1),kind=GP)**2+aimag(V3(k,j,1))**2)*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc + 2*(real(V1(k,j,i),kind=GP)**2+aimag(V1(k,j,i))**2+ &
                                   real(V2(k,j,i),kind=GP)**2+aimag(V2(k,j,i))**2+ &
                                   real(V3(k,j,i),kind=GP)**2+aimag(V3(k,j,i))**2)*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc + 2*(real(V1(k,j,i),kind=GP)**2+aimag(V1(k,j,i))**2+ &
                                   real(V2(k,j,i),kind=GP)**2+aimag(V2(k,j,i))**2+ &
                                   real(V3(k,j,i),kind=GP)**2+aimag(V3(k,j,i))**2)*tmp
               END DO
            END DO
         END DO
      ENDIF

! Computes the reduction between nodes
      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      CALL gws%free_complex_htmp(C1)
      CALL gws%free_complex_htmp(C2)
      CALL gws%free_complex_htmp(C3)
      CALL gws%free_complex_htmp(C4)
      RETURN
      END SUBROUTINE energy

!*****************************************************************
      SUBROUTINE helicity(a,b,c,d)
!-----------------------------------------------------------------
!
! Computes the mean helicity of a vector field.
! The output is only valid in the first node.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: at the output contains the helicity
!
      USE fprecision
      USE commtypes
      USE fft
      USE grid
      USE mpivars
      USE class_GWorkspace3D, ONLY: gws
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      DOUBLE PRECISION, INTENT(OUT) :: d

      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c4,c5,c6

      DOUBLE PRECISION :: dloc
      REAL(KIND=GP)    :: tmp
      INTEGER          :: i,j,k
      LOGICAL          :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
      CALL gws%get_complex_htmp(c4,bret)
      CALL gws%get_complex_htmp(c5,bret)
      CALL gws%get_complex_htmp(c6,bret)
      dloc = 0.0D0
      ! One nx*ny from parsevals, an nx*ny*(nz-Cz) from averaging over
      ! physical domain and a nz**2 due to fft normalization.
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2 / &
            real(nz-Cz,kind=GP)

      ! Computes the field and its curl in (z,ky,kx)
      CALL copy(a,c1)
      CALL copy(b,c2)
      CALL copy(c,c3)
      CALL curlk(b,c,c4,1)
      CALL curlk(a,c,c5,2)
      CALL curlk(a,b,c6,3)
      CALL fftp1d_complex_to_real_z(planfc,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c2,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c3,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c4,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c5,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c6,MPI_COMM_WORLD)

      ! Compute average values
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               dloc = dloc + (real(c1(k,j,1)*conjg(c4(k,j,1)))+ &
                              real(c2(k,j,1)*conjg(c5(k,j,1)))+ &
                              real(c3(k,j,1)*conjg(c6(k,j,1))))*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc + 2*(real(c1(k,j,i)*conjg(c4(k,j,i)))+ &
                                   real(c2(k,j,i)*conjg(c5(k,j,i)))+ &
                                   real(c3(k,j,i)*conjg(c6(k,j,i))))*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc + 2*(real(c1(k,j,i)*conjg(c4(k,j,i)))+ &
                                   real(c2(k,j,i)*conjg(c5(k,j,i)))+ &
                                   real(c3(k,j,i)*conjg(c6(k,j,i))))*tmp
               END DO
            END DO
         END DO
      ENDIF

      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

      CALL gws%free_complex_htmp(c1)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c3)
      CALL gws%free_complex_htmp(c4)
      CALL gws%free_complex_htmp(c5)
      CALL gws%free_complex_htmp(c6)
      RETURN
      END SUBROUTINE helicity

!*****************************************************************
      SUBROUTINE cross(a,b,c,d,e,f,g,kin)
!-----------------------------------------------------------------
!
! Computes the cross helicity or averaged inner
! product of two vector fields. The output is
! only valid in the first node.
!
! Parameters
!     a  : first field x-component
!     b  : first field y-component
!     c  : first field z-component
!     d  : second field x-component
!     e  : second field y-component
!     f  : second field z-component
!     g  : at the output contains the inner product
!     kin: =0 computes the inner product of the curls
!          =1 computes the inner product of the fields
!
      USE fprecision
      USE commtypes
      USE kes
      USE fft
      USE grid
      USE mpivars
      USE class_GWorkspace3D, ONLY: gws
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      DOUBLE PRECISION, INTENT(OUT) :: g

      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c4,c5,c6

      DOUBLE PRECISION              :: gloc
      REAL(KIND=GP)                 :: tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      LOGICAL             :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
      CALL gws%get_complex_htmp(c4,bret)
      CALL gws%get_complex_htmp(c5,bret)
      CALL gws%get_complex_htmp(c6,bret)
      gloc = 0.0D0
      ! One nx*ny from parsevals, an nx*ny*(nz-Cz) from averaging over
      ! physical domain and a nz**2 due to fft normalization.
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2/ &
            real(nz-Cz, kind=GP)

      ! Computes the fields (or their curls) in (z,ky,kx)
      IF (kin.eq.1) THEN
         CALL copy(a,c1)
         CALL copy(b,c2)
         CALL copy(c,c3)
         CALL copy(d,c4)
         CALL copy(e,c5)
         CALL copy(f,c6)
      ELSE
         CALL curlk(b,c,c1,1)
         CALL curlk(a,c,c2,2)
         CALL curlk(a,b,c3,3)
         CALL curlk(e,f,c4,1)
         CALL curlk(d,f,c5,2)
         CALL curlk(d,e,c6,3)
      ENDIF
      CALL fftp1d_complex_to_real_z(planfc,c1,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c2,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c3,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c4,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c5,MPI_COMM_WORLD)
      CALL fftp1d_complex_to_real_z(planfc,c6,MPI_COMM_WORLD)

      ! Compute average value
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:gloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               gloc = gloc + (DBLE(real(c1(k,j,1)*conjg(c4(k,j,1))))+ &
                              real(c2(k,j,1)*conjg(c5(k,j,1)))+ &
                              real(c3(k,j,1)*conjg(c6(k,j,1))))*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:gloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:gloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  gloc = gloc + 2*(DBLE(real(c1(k,j,i)*conjg(c4(k,j,i))))+ &
                                   real(c2(k,j,i)*conjg(c5(k,j,i)))+ &
                                   real(c3(k,j,i)*conjg(c6(k,j,i))))*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:gloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:gloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  gloc = gloc + 2*(DBLE(real(c1(k,j,i)*conjg(c4(k,j,i))))+ &
                                   real(c2(k,j,i)*conjg(c5(k,j,i)))+ &
                                   real(c3(k,j,i)*conjg(c6(k,j,i))))*tmp
               END DO
            END DO
         END DO
      ENDIF

      ! Computes the reduction between nodes
      CALL MPI_REDUCE(gloc,g,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                          MPI_COMM_WORLD,ierr)

      CALL gws%free_complex_htmp(c1)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c3)
      CALL gws%free_complex_htmp(c4)
      CALL gws%free_complex_htmp(c5)
      CALL gws%free_complex_htmp(c6)
      RETURN
      END SUBROUTINE cross

!*****************************************************************
      SUBROUTINE hdcheck(a,b,c,d,e,f,t,dt,hel)
!-----------------------------------------------------------------
!
!  Computes the kinetic energy and helicity of the fluid.
!
! Output files contain:
! 'balance.txt':  time, <v^2>, <omega^2>, mechanic injection rate
! 'helicity.txt': time, kinetic helicity
!
! Parameters
!     a  : velocity field in the x-direction
!     b  : velocity field in the y-direction
!     c  : velocity field in the z-direction
!     d  : force in the x-direction
!     e  : force in the y-direction
!     f  : force in the z-direction
!     t  : number of time steps made
!     dt : time step
!     hel: =0 skips kinetic helicity computation
!          =1 computes the kinetic helicity
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: d,e,f
      DOUBLE PRECISION    :: eng,ens,pot,khe
      REAL(KIND=GP)       :: dt
      INTEGER, INTENT(IN) :: hel
      INTEGER, INTENT(IN) :: t
      INTEGER             :: i,j,k

! Computes the mean energy, enstrophy, and kinetic helicity
      CALL energy(a,b,c,eng,1)
      CALL energy(a,b,c,ens,0)
      IF (hel.eq.1) THEN
         CALL helicity(a,b,c,khe)
      ENDIF

! Computes the energy injection rate
      CALL cross(a,b,c,d,e,f,pot,1)

! Creates external files to store the results
      IF (myrank.eq.0) THEN
         OPEN(1,file='balance.txt',position='append')
   10       FORMAT( 1P E13.6,2E23.16,E24.16 )
            WRITE(1,10) (t-1)*dt,eng,ens,pot
         CLOSE(1)
         IF (hel.eq.1) THEN
            OPEN(1,file='helicity.txt',position='append')
   11          FORMAT( 1P E13.6,E24.16 )
               WRITE(1,11) (t-1)*dt,khe
            CLOSE(1)
         ENDIF
      ENDIF

      RETURN
      END SUBROUTINE hdcheck

!*****************************************************************
      SUBROUTINE maxabs(a,b,c,d,kin)
!-----------------------------------------------------------------
!
! Computes the maximum absolute value of the field 
! vorticity, current density, or of the original field. 
! The output is only valid in the first node.
!
! Parameters
!     a  : field x-component
!     b  : field y-component
!     c  : field z-component
!     d  : at the output contains the maximum value
!     kin: =0 computes the maximum of vorticity
!          =1 computes the maximum of current density
!          =2 computes the maximum of the field
!
      USE fprecision
      USE commtypes
      USE fft
      USE grid
      USE mpivars
      USE class_GWorkspace3D, ONLY: gws
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: c1,c2,c3
      REAL(KIND=GP), POINTER, DIMENSION(:,:,:)    :: r1,r2,r3
      REAL(KIND=GP), INTENT(OUT)   :: d
      REAL(KIND=GP)                :: dloc
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: i,j,k
      LOGICAL             :: bret

      CALL gws%get_complex_htmp(c1,bret)
      CALL gws%get_complex_htmp(c2,bret)
      CALL gws%get_complex_htmp(c3,bret)
      CALL gws%get_real_htmp(r1,bret)
      CALL gws%get_real_htmp(r2,bret)
      CALL gws%get_real_htmp(r3,bret)
      IF (kin.eq.0) THEN
         CALL curlk(b,c,c1,1)
         CALL curlk(a,c,c2,2)
         CALL curlk(a,b,c3,3)
      ELSE IF (kin.eq.1) THEN
         CALL laplak(a,c1)
         CALL laplak(b,c2)
         CALL laplak(c,c3)
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  c1(k,j,i) = a(k,j,i)
                  c2(k,j,i) = b(k,j,i)
                  c3(k,j,i) = c(k,j,i)
               END DO
            END DO
         END DO
      ENDIF
      CALL fftp3d_complex_to_real(planfc,c1,r1,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c2,r2,MPI_COMM_WORLD)
      CALL fftp3d_complex_to_real(planfc,c3,r3,MPI_COMM_WORLD)
      dloc = 0.0_GP
!$omp parallel do if (pkend-ksta.ge.nth) private (j,i) reduction(max:dloc)
      DO k = ksta,pkend
!$omp parallel do if (pkend-ksta.lt.nth) private (i) reduction(max:dloc)
         DO j = 1,ny
            DO i = 1,nx
               dloc = max(dloc, r1(i,j,k)**2+r2(i,j,k)**2+r3(i,j,k)**2)
            END DO
         END DO
      END DO
      dloc = sqrt(dloc)/(real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))
      CALL MPI_REDUCE(dloc,d,1,GC_REAL,MPI_MAX,0,MPI_COMM_WORLD,ierr)

      CALL gws%free_complex_htmp(c1)
      CALL gws%free_complex_htmp(c2)
      CALL gws%free_complex_htmp(c3)
      CALL gws%free_real_htmp(r1)
      CALL gws%free_real_htmp(r2)
      CALL gws%free_real_htmp(r3)
      RETURN
      END SUBROUTINE maxabs

!*****************************************************************
      SUBROUTINE saxpby_c(z,x,a,y,b)
!-----------------------------------------------------------------
!
! Computes z = ax + by for complex x, y. Vale en (lm,r) y en
! (lm,kz). z puede ser el mismo arreglo que x o que y, porque la
! cuenta es punto a punto.
!
! Parameters
!     z     : output
!     x,y   : input arrays
!     a, b  : constants
!
      USE fprecision
      USE grid, ONLY: nz
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,nz) :: z
      COMPLEX(KIND=GP), INTENT   (IN), DIMENSION(lmsta:lmend,nz) :: x,y
      REAL   (KIND=GP), INTENT   (IN)                            :: a,b
      INTEGER                                                    :: lm,k

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active)
#else
!$omp parallel do private (lm)
#endif
      DO k = 1,nz
         DO lm = lmsta,lmend
            z(lm,k) = a*x(lm,k) + b*y(lm,k)
         END DO
      END DO

      RETURN
      END SUBROUTINE saxpby_c

!*****************************************************************
      SUBROUTINE saxpby_r(z,x,a,y,b)
!-----------------------------------------------------------------
!
! Computes z = ax + by for real x, y, en el espacio físico
! (theta,phi,r). z puede ser el mismo arreglo que x o que y,
! porque la cuenta es punto a punto.
!
! Parameters
!     z     : output
!     x,y   : input arrays
!     a, b  : constants
!
      USE fprecision
      USE grid, ONLY: nlat,nphi
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      REAL(KIND=GP), INTENT(INOUT), DIMENSION(nlat,nphi,ksta:kend) :: z
      REAL(KIND=GP), INTENT   (IN), DIMENSION(nlat,nphi,ksta:kend) :: x,y
      REAL(KIND=GP), INTENT   (IN)                                 :: a,b
      INTEGER                                                      :: i,j,k

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(3) if(target: gdev_active)
#else
!$omp parallel do collapse(2) private (i)
#endif
      DO k = ksta,kend
         DO j = 1,nphi
            DO i = 1,nlat
               z(i,j,k) = a*x(i,j,k) + b*y(i,j,k)
            END DO
         END DO
      END DO

      RETURN
      END SUBROUTINE saxpby_r

!*****************************************************************
      SUBROUTINE fc_filter(a,esp)
!-----------------------------------------------------------------
!
! Filters a field along z (the radius) using exponential tails,
! with the formula a*exp(-alpha(2 kz/ N)**2p), en el espacio que
! indica esp: llama a fc_filter_r en (lm,r) o a fc_filter_k en
! (lm,kz). En los ángulos no hace falta: el no lineal se trunca en
! lmax, con la regla de 3/2 de la grilla.
!
! Parameters
!     a    : input/output array
!     esp  : espacio de a, 'r' para (lm,r) o 'k' para (lm,kz)
!
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,nz) :: a
      CHARACTER(len=1), INTENT(IN) :: esp

      SELECT CASE (esp)
      CASE ('r')
         CALL fc_filter_r(a)
      CASE ('k')
         CALL fc_filter_k(a)
      CASE DEFAULT
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Unknown space ", &
             esp, " in fc_filter. Aborting..."
         CALL MPI_FINALIZE(ierr)
         STOP
      END SELECT

      RETURN
      END SUBROUTINE fc_filter

!*****************************************************************
      SUBROUTINE fc_filter_k(a)
!-----------------------------------------------------------------
!
! Filters a field in (lm,kz) along z (the radius) using exponential
! tails, with the formula a*exp(-alpha(2 kz/ N)**2p).
!
! Parameters
!     a    : input/output array

      USE fprecision
      USE kes
      USE mpivars
      USE grid
      USE boxsize
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,nz) :: a
      DOUBLE PRECISION, PARAMETER  :: alpha=16*log(10d0), p=50d0
      INTEGER :: lm,k

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active)
#else
!$omp parallel do private (lm)
#endif
      DO k = 1,nz
         DO lm = lmsta,lmend
            a(lm,k) = a(lm,k)*exp(-alpha*(2*kz(k)/nz/Dkz)**(2*p))
         END DO
      END DO

      RETURN
      END SUBROUTINE fc_filter_k

!*****************************************************************
      SUBROUTINE fc_filter_r(a)
!-----------------------------------------------------------------
!
! Filters a field in (lm,r) along z (the radius): ida a (lm,kz)
! con FC-Gram, fc_filter_k y vuelta, en el mismo arreglo (el
! filtro se aplica en el lugar, así que no hace falta copia).
!
! Parameters
!     a    : input/output array in (lm,r)

      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,nz) :: a

      CALL fftp1d_real_to_complex_z(planfc,a,MPI_COMM_WORLD)
      CALL fc_filter_k(a)
      CALL fftp1d_complex_to_real_z(planfc,a,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE fc_filter_r

!*****************************************************************
      SUBROUTINE copy(a,b)
!-----------------------------------------------------------------
!
! Copies the complex matrix 'a' into 'b', b = a. Vale en (lm,r) y
! en (lm,kz).
!
      USE fprecision
      USE grid
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT (IN), DIMENSION(lmsta:lmend,nz) :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz) :: b
      INTEGER :: lm,k

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do collapse(2) if(target: gdev_active)
#else
!$omp parallel do private (lm)
#endif
      DO k = 1,nz
         DO lm = lmsta,lmend
            b(lm,k) = a(lm,k)
         END DO
      END DO
      RETURN
      END SUBROUTINE copy

!*****************************************************************
      SUBROUTINE setmode(a,k,j,i,val)
!-----------------------------------------------------------------
!
! Sets one element of the matrix 'a', a(k,j,i) = val, in the
! copy of the matrix being worked on (the device copy while
! gdev_active is set)
!
      USE fprecision
      USE grid
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a
      COMPLEX(KIND=GP), INTENT(IN) :: val
      INTEGER, INTENT(IN) :: k,j,i

#if defined(GHOST_GPU)
!$omp target if(target: gdev_active)
      a(k,j,i) = val
!$omp end target
#else
      a(k,j,i) = val
#endif
      RETURN
      END SUBROUTINE setmode

!*****************************************************************
      SUBROUTINE divergence(a,b,c,d)
!-----------------------------------------------------------------
!
! Computes the mean square divergence of a vector field.
! The output is only valid in the first node.
!
! Parameters
!     a: input matrix in the x-direction
!     b: input matrix in the y-direction
!     c: input matrix in the z-direction
!     d: at the output contains the divergence
!
      USE var
      USE grid
      USE fft
      USE mpivars
      USE commtypes
      USE class_GWorkspace3D, ONLY: gws
!$    USE threads
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(nz,ny,ista:iend) :: a,b,c
      DOUBLE PRECISION, INTENT(OUT) :: d

      COMPLEX(KIND=GP), POINTER, DIMENSION(:,:,:) :: C1,C2

      DOUBLE PRECISION  :: dloc
      REAL(KIND=GP)     :: tmp

      INTEGER :: i,j,k
      LOGICAL :: bret

      CALL gws%get_complex_htmp(C1,bret)
      CALL gws%get_complex_htmp(C2,bret)
      ! One nx*ny from parsevals, an nx*ny*(nz-Cz) from averaging over
      ! physical domain and a nz**2 due to fft normalization.
      tmp = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2/ &
            real(nz-Cz,kind=GP)
      d = 0d0
      dloc = 0d0

      ! Get pointwise divergence
      CALL derivk(a,C1,1)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C2(k,j,i) = C1(k,j,i)
            ENDDO
         ENDDO
      ENDDO

      CALL derivk(b,C1,2)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C2(k,j,i) = C2(k,j,i) + C1(k,j,i)
            ENDDO
         ENDDO
      ENDDO

      CALL derivk(c,C1,3)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               C2(k,j,i) = C2(k,j,i) + C1(k,j,i)
            ENDDO
         ENDDO
      ENDDO

      CALL fftp1d_complex_to_real_z(planfc,C2,MPI_COMM_WORLD)

      ! Compute mean value
      IF (ista.eq.1) THEN
!$omp parallel do private (k) reduction(+:dloc)
         DO j = 1,ny
            DO k = 1,nz-Cz
               dloc = dloc + (real(C2(k,j,1),kind=GP)**2 + aimag(C2(k,j,1))**2)*tmp
            END DO
         END DO
!$omp parallel do if (iend-2.ge.nth) private (j,k) reduction(+:dloc)
         DO i = 2,iend
!$omp parallel do if (iend-2.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc + 2*(real(C2(k,j,i),kind=GP)**2 + aimag(C2(k,j,i))**2)*tmp
               END DO
            END DO
         END DO
      ELSE
!$omp parallel do if (iend-ista.ge.nth) private (j,k) reduction(+:dloc)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k) reduction(+:dloc)
            DO j = 1,ny
               DO k = 1,nz-Cz
                  dloc = dloc + 2*(real(C2(k,j,i),kind=GP)**2 + aimag(C2(k,j,i))**2)*tmp
               END DO
            END DO
         END DO
      ENDIF

      CALL MPI_REDUCE(dloc,d,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,&
          MPI_COMM_WORLD,ierr)

      CALL gws%free_complex_htmp(C1)
      CALL gws%free_complex_htmp(C2)
      RETURN
      END SUBROUTINE divergence

!*****************************************************************
      SUBROUTINE normvec(a,b,c,d,kin)
!-----------------------------------------------------------------
!
! Normalizes a vector field.
!
! Parameters
!     a  : input matrix in the x-direction
!     b  : input matrix in the y-direction
!     c  : input matrix in the z-direction
!     d  : input specifying the desired value
!     kin: =0 normalize to specified magnetic energy
!          =1 normalize to specified kinetic energy
!          =2 normalize to specified magnetic enstrophy
!

      USE fprecision
      USE grid
      USE commtypes
      USE mpivars
!$    USE threads

      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(nz,ny,ista:iend) :: a,b,c
      REAL(KIND=GP), INTENT(IN) :: d
      INTEGER, INTENT(IN)       :: kin

      DOUBLE PRECISION  :: tmp
      REAL(KIND=GP)     :: rmp
      INTEGER :: i,j,k

      CALL energy(a,b,c,tmp,kin)
      CALL MPI_BCAST(tmp,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      rmp = sqrt(d/tmp)

!$omp parallel do if (iend-ista.ge.nth) private (j,k)
      DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
            DO k = 1,nz
               a(k,j,i) = a(k,j,i)*rmp
               b(k,j,i) = b(k,j,i)*rmp
               c(k,j,i) = c(k,j,i)*rmp
            ENDDO
        ENDDO
      ENDDO

      RETURN
      END SUBROUTINE normvec
