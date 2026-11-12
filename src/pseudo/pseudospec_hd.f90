!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Rutinas pseudoespectrales de la cáscara: derivadas radiales con
! FC-Gram, el laplaciano angular, potencias de r, el operador D_l,
! el producto vectorial en el espacio físico, combinaciones
! lineales, el filtro radial y copias. Los campos están en (lm,r) o
! en (lm,kz), con forma (lmsta:lmend,nz), o en el espacio físico,
! con forma (nlat,nphi,ksta:kend). Hay que usar los módulos de
! sphp_mod.F90 y pseudospec_mod.f90 en cada programa que las llame.
!
! Adaptado de pseudospec_hd de SPECTER, 2019 Mauro Fontana & Pablo
! D. Mininni. Department of Physics, Facultad de Ciencias Exactas y
! Naturales, Universidad de Buenos Aires.
!=================================================================

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
