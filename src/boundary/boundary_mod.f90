!======================================================================
! Condiciones de borde de la cáscara, en las paredes r = rint y
! r = rext.
!
! Adaptado del boundary de SPECTER, 2020 Mauro Fontana. DF-UBA.
!======================================================================

MODULE boundary
      USE fprecision

      IMPLICIT NONE

      ! Tipo de condición en cada pared, rint y rext (ver v_parsebc)
      TYPE BCPLAN
         INTEGER  :: bcint,bcext
      ENDTYPE BCPLAN

      ! Factorizaciones LU de las matrices B_l de la recuperación de W
      ! (invDl), una por cada l de 1 a lmax, con sus pivotes y la
      ! escala de cada fila
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:) :: blu
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:)   :: besc
      INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: bpiv

      CONTAINS

      INCLUDE 'vboundary.f90'

!***************************************************************************
      SUBROUTINE setup_bc(this,planfc,bckind,field)
!---------------------------------------------------------------------------
!  Add field and boundary conditions to the plan
!  ARGUMENTS:
!    this   : 'this' instance of a BCPLAN. [INOUT]
!    planfc : FCPLAN that will be used for the boundary conditions. [INOUT]
!    bckind : tipo de condición en cada pared, rint y rext. [IN]
!    field  : Field whose boundary conditions must be set up. Only
!             'v' (the velocity) is available.
!---------------------------------------------------------------------------
      USE fcgram

      IMPLICIT NONE

      TYPE(BCPLAN), INTENT(INOUT)          :: this
      TYPE(FCPLAN), INTENT(INOUT)          :: planfc
      CHARACTER(len=15), INTENT(IN)        :: bckind(2)
      CHARACTER(len=*), INTENT(IN)         :: field

      IF ( field .eq. 'v' ) THEN
         CALL v_setup(this,planfc,bckind)
      ENDIF

      RETURN
      ENDSUBROUTINE setup_bc

!**********************************************************************
      SUBROUTINE dirichlet(a)
!----------------------------------------------------------------------
!  Anula el campo 'a' en (lm,r) en las dos paredes, r = rint (z(1)) y
! r = rext (z(nz-Cz)), en todos los lm de la tarea. Es la condición
! de Z, el no deslizamiento de la parte toroidal, y sirve para
! cualquier campo que valga cero en las paredes. Los puntos de la
! continuación no se tocan: la FFT en z los vuelve a calcular a
! partir de los puntos físicos.
!  ARGUMENTS:
!    a : campo en (lm,r). [INOUT]
!----------------------------------------------------------------------
      USE fprecision
      USE grid
      USE mpivars
      USE gdevice, ONLY: gdev_active
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,nz) :: a
      INTEGER :: lm

#if defined(GHOST_GPU)
!$omp target teams distribute parallel do if(target: gdev_active)
#else
!$omp parallel do
#endif
      DO lm = lmsta,lmend
         a(lm,1)     = 0.0_GP
         a(lm,nz-Cz) = 0.0_GP
      END DO

      RETURN
      END SUBROUTINE dirichlet

!**********************************************************************
      SUBROUTINE derivr_mat(p,d)
!----------------------------------------------------------------------
!  Matriz de la derivada radial de orden p de FC-Gram en los N = nz-Cz
! puntos físicos: d(i,j) es la derivada en z(i) del vector canónico
! e_j, así que d aplicada a los valores de un campo en los puntos
! físicos da lo mismo que derivr. Las columnas se calculan con
! derivr, de a tantas como lm tiene la tarea (una por fila de un
! arreglo en (lm,r)), y cada tarea arma la matriz entera, sin
! comunicar. Se llama una vez, al armar las matrices de borde, con
! gdev_active falso.
!  ARGUMENTS:
!    p : orden de la derivada. [IN]
!    d : la matriz, de N x N. [OUT]
!----------------------------------------------------------------------
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: p
      REAL(KIND=GP), INTENT(OUT), DIMENSION(nz-Cz,nz-Cz) :: d
      COMPLEX(KIND=GP), ALLOCATABLE, DIMENSION(:,:) :: c1,c2
      INTEGER :: n,j0,j,lm

      n = nz-Cz
      ALLOCATE( c1(lmsta:lmend,nz),c2(lmsta:lmend,nz) )
      DO j0 = 1,n,lmend-lmsta+1
         c1 = 0.0_GP
         DO lm = lmsta,lmend
            j = j0+lm-lmsta
            IF ( j .le. n ) c1(lm,j) = 1.0_GP
         END DO
         CALL derivr(c1,c2,p,'r')
         DO lm = lmsta,lmend
            j = j0+lm-lmsta
            IF ( j .le. n ) d(:,j) = REAL(c2(lm,1:n),KIND=GP)
         END DO
      END DO
      DEALLOCATE( c1,c2 )

      RETURN
      END SUBROUTINE derivr_mat

!**********************************************************************
      SUBROUTINE invDl_setup()
!----------------------------------------------------------------------
!  Arma y factoriza, para cada l de 1 a lmax, la matriz B_l de la
! recuperación de W a partir de phi = D_l W, en los N = nz-Cz puntos
! físicos (cuenta 4, ec. 14):
!    fila 1        : W = 0 en rint
!    fila 2        : dW/dr = 0 en rint (fila 1 de la matriz de d/dr)
!    filas 3 a N-2 : D_l W = phi (filas de la matriz de d2/dr2, menos
!                    l(l+1)/r**2 en la diagonal)
!    fila N-1      : dW/dr = 0 en rext (fila N de la matriz de d/dr)
!    fila N        : W = 0 en rext
! Cada fila se divide por su máximo valor absoluto (besc) y la matriz
! se factoriza con LU (blu, bpiv). Cada tarea guarda todos los l,
! porque sus lm tienen casi todos. Se llama una vez, con gdev_active
! falso.
!----------------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE grid
      USE mpivars
      IMPLICIT NONE

      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:) :: d1,d2
      INTEGER :: n,l,i,info

      n = nz-Cz
      ALLOCATE( d1(n,n),d2(n,n) )
      ALLOCATE( blu(n,n,lmax),besc(n,lmax),bpiv(n,lmax) )
      CALL derivr_mat(1,d1)
      CALL derivr_mat(2,d2)

      DO l = 1,lmax
         blu(:,:,l) = 0.0_GP
         blu(1,1,l) = 1.0_GP
         blu(2,:,l) = d1(1,:)
         DO i = 3,n-2
            blu(i,:,l) = d2(i,:)
            blu(i,i,l) = blu(i,i,l) - l*(l+1)/z(i)**2
         END DO
         blu(n-1,:,l) = d1(n,:)
         blu(n,n,l) = 1.0_GP
         DO i = 1,n
            besc(i,l) = MAXVAL(ABS(blu(i,:,l)))
            blu(i,:,l) = blu(i,:,l)/besc(i,l)
         END DO
#if defined(GDOUBLE_PRECISION)
         CALL DGETRF(n,n,blu(:,:,l),n,bpiv(:,l),info)
#else
         CALL SGETRF(n,n,blu(:,:,l),n,bpiv(:,l),info)
#endif
         IF ( info .ne. 0 ) THEN
            IF ( myrank .eq. 0 ) PRINT*, "[ERROR] LU of B_l failed ", &
                "for l = ", l, " in invDl_setup. Aborting..."
            CALL MPI_FINALIZE(ierr)
            STOP
         ENDIF
      END DO
      DEALLOCATE( d1,d2 )

      RETURN
      END SUBROUTINE invDl_setup

!**********************************************************************
      SUBROUTINE invDl(a,b)
!----------------------------------------------------------------------
!  Recupera W a partir de phi = D_l W en (lm,r): para cada lm de la
! tarea resuelve B_l W = P phi con la LU de invDl_setup, donde P phi
! es phi con ceros en las filas 1, 2, N-1 y N (las de las condiciones
! W = dW/dr = 0 en las paredes) y cada fila dividida por besc. Las
! partes real e imaginaria van como dos columnas del lado derecho,
! porque B_l es real. El modo l = 0 queda en cero (libertad de gauge),
! y los puntos de la continuación también: la FFT en z los vuelve a
! calcular. Corre en la CPU.
!  ARGUMENTS:
!    a : phi en (lm,r). [IN]
!    b : W en (lm,r). [OUT]
!----------------------------------------------------------------------
      USE fprecision
      USE grid
      USE kes
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(lmsta:lmend,nz)  :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz) :: b
      REAL(KIND=GP), DIMENSION(nz-Cz,2) :: rhs
      INTEGER :: n,lm,l,info

      n = nz-Cz
      b = 0.0_GP
      DO lm = lmsta,lmend
         l = ell(lm)
         IF ( l .eq. 0 ) CYCLE
         rhs(:,1) = REAL(a(lm,1:n),KIND=GP)/besc(:,l)
         rhs(:,2) = AIMAG(a(lm,1:n))/besc(:,l)
         rhs(1:2,:) = 0.0_GP
         rhs(n-1:n,:) = 0.0_GP
#if defined(GDOUBLE_PRECISION)
         CALL DGETRS('N',n,2,blu(:,:,l),n,bpiv(:,l),rhs,n,info)
#else
         CALL SGETRS('N',n,2,blu(:,:,l),n,bpiv(:,l),rhs,n,info)
#endif
         b(lm,1:n) = CMPLX(rhs(:,1),rhs(:,2),KIND=GP)
      END DO

      RETURN
      END SUBROUTINE invDl

!***************************************************************************
      FUNCTION preprocess(string) RESULT(output)
!---------------------------------------------------------------------------
!  Returns a string preprocessed and ready for comparison with another
! alphabetic lower case string.
!---------------------------------------------------------------------------
      IMPLICIT NONE
         CHARACTER(len=*), INTENT(in) :: string
         CHARACTER(len=len(string))   :: output

         output = lowcase(trim(adjustl(string)))

      END FUNCTION preprocess

!***************************************************************************
      FUNCTION lowcase(string) RESULT(lower)
!---------------------------------------------------------------------------
!  Returns a lower case version of a string. Assumes all characters are 
!  ASCII encoded. Based on: https://www.star.le.ac.uk/%7ecgp/fortran.html
!---------------------------------------------------------------------------
         IMPLICIT NONE
         CHARACTER(len=*), INTENT(in) :: string
         CHARACTER(len=len(string))   :: lower
         INTEGER :: j

         DO j = 1,len(string)
         IF(string(j:j) >= "A" .AND. string(j:j) <= "Z") then
            lower(j:j) = achar(iachar(string(j:j)) + 32)
         ELSE
            lower(j:j) = string(j:j)
         END IF
         END DO
      END FUNCTION lowcase

ENDMODULE boundary
