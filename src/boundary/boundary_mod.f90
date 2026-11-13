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

      ! Para cada l de 1 a lmax, M_l: las columnas 3 a N-2 de B_l^-1
      ! (binv), con B_l la matriz de la recuperación de W, así que
      ! W = M_l phi(3:N-2) (invDl)
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:,:) :: binv
      ! Los lm de la tarea ordenados por l (lmord), con los de cada l
      ! entre lmini(l) y lmini(l+1)-1, y los arreglos de trabajo de invDl
      INTEGER, ALLOCATABLE, DIMENSION(:)         :: lmini,lmord
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:) :: bphi,bw

      ! Las filas 1, 2, N-1 y N de la matriz de d2/dr2 (v_imposebc_fast)
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:) :: d2bor

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
!  Arma, para cada l de 1 a lmax, la matriz B_l de la recuperación de W
! a partir de phi = D_l W, en los N = nz-Cz puntos físicos (cuenta 4,
! ec. 14):
!    fila 1        : W = 0 en rint
!    fila 2        : dW/dr = 0 en rint (fila 1 de la matriz de d/dr)
!    filas 3 a N-2 : D_l W = phi (filas de la matriz de d2/dr2, menos
!                    l(l+1)/r**2 en la diagonal)
!    fila N-1      : dW/dr = 0 en rext (fila N de la matriz de d/dr)
!    fila N        : W = 0 en rext
! Divide cada fila por su máximo valor absoluto (la diagonal S),
! factoriza S^-1 B_l con LU y con ella calcula M_l, las columnas 3 a
! N-2 de B_l^-1 (binv), las únicas que usa P phi: así invDl es un
! producto de matrices. El escalado solo mejora el redondeo de la LU,
! porque (S^-1 B_l)^-1 S^-1 = B_l^-1. También arma la tabla de los lm
! de la tarea por l y guarda las filas de borde de d2/dr2. Cada tarea
! guarda todos los l, porque sus lm tienen casi todos. Se llama una
! vez, con gdev_active falso.
!----------------------------------------------------------------------
      USE fprecision
      USE commtypes
      USE grid
      USE kes
      USE mpivars
      IMPLICIT NONE

      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:,:) :: d1,d2,b
      REAL(KIND=GP), ALLOCATABLE, DIMENSION(:)   :: esc
      INTEGER, ALLOCATABLE, DIMENSION(:)         :: piv
      INTEGER :: n,l,i,info,j,lm

      n = nz-Cz
      ALLOCATE( d1(n,n),d2(n,n),b(n,n),esc(n),piv(n) )
      ALLOCATE( binv(n,n-4,lmax) )
      CALL derivr_mat(1,d1)
      CALL derivr_mat(2,d2)

      DO l = 1,lmax
         b = 0.0_GP
         b(1,1) = 1.0_GP
         b(2,:) = d1(1,:)
         DO i = 3,n-2
            b(i,:) = d2(i,:)
            b(i,i) = b(i,i) - l*(l+1)/z(i)**2
         END DO
         b(n-1,:) = d1(n,:)
         b(n,n) = 1.0_GP
         DO i = 1,n
            esc(i) = MAXVAL(ABS(b(i,:)))
            b(i,:) = b(i,:)/esc(i)
         END DO
#if defined(GDOUBLE_PRECISION)
         CALL DGETRF(n,n,b,n,piv,info)
#else
         CALL SGETRF(n,n,b,n,piv,info)
#endif
         IF ( info .ne. 0 ) THEN
            IF ( myrank .eq. 0 ) PRINT*, "[ERROR] LU of B_l failed ", &
                "for l = ", l, " in invDl_setup. Aborting..."
            CALL MPI_FINALIZE(ierr)
            STOP
         ENDIF
! Las columnas 3 a N-2 de B_l^-1: S^-1 B_l M_l = esas columnas de S^-1
         binv(:,:,l) = 0.0_GP
         DO j = 1,n-4
            binv(j+2,j,l) = 1.0_GP/esc(j+2)
         END DO
#if defined(GDOUBLE_PRECISION)
         CALL DGETRS('N',n,n-4,b,n,piv,binv(:,:,l),n,info)
#else
         CALL SGETRS('N',n,n-4,b,n,piv,binv(:,:,l),n,info)
#endif
      END DO

! Los lm de la tarea ordenados por l, para resolver juntos los que
! comparten B_l
      ALLOCATE( lmini(0:lmax+1),lmord(lmend-lmsta+1) )
      ALLOCATE( bphi(n-4,2*(lmax+1)),bw(n,2*(lmax+1)) )
      j = 0
      DO l = 0,lmax
         lmini(l) = j+1
         DO lm = lmsta,lmend
            IF ( ell(lm) .eq. l ) THEN
               j = j+1
               lmord(j) = lm
            ENDIF
         END DO
      END DO
      lmini(lmax+1) = j+1

! Las filas de d2/dr2 de los cuatro puntos de borde, para recalcular phi
! solo ahí (v_imposebc_fast)
      ALLOCATE( d2bor(4,n) )
      d2bor(1,:) = d2(1,:)
      d2bor(2,:) = d2(2,:)
      d2bor(3,:) = d2(n-1,:)
      d2bor(4,:) = d2(n,:)
      DEALLOCATE( d1,d2,b,esc,piv )

      RETURN
      END SUBROUTINE invDl_setup

!**********************************************************************
      SUBROUTINE invDl(a,b)
!----------------------------------------------------------------------
!  Recupera W a partir de phi = D_l W en (lm,r): W = M_l phi(3:N-2),
! con M_l las columnas 3 a N-2 de B_l^-1 (invDl_setup), cumple
! D_l W = phi en los puntos 3 a N-2 y W = dW/dr = 0 en las paredes.
! Los lm de la tarea con el mismo l comparten M_l y van juntos en un
! solo DGEMM, con dos columnas por modo: las partes real e imaginaria,
! porque M_l es real.
! El modo l = 0 queda en cero (libertad de gauge), y los puntos de la
! continuación también: la FFT en z los vuelve a calcular. Corre en la
! CPU.
!  ARGUMENTS:
!    a : phi en (lm,r). [IN]
!    b : W en (lm,r). [OUT]
!----------------------------------------------------------------------
      USE fprecision
      USE grid
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(lmsta:lmend,nz)  :: a
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz) :: b
      INTEGER :: n,l,i,j,nm

      n = nz-Cz
      b = 0.0_GP
      DO l = 1,lmax
         nm = lmini(l+1)-lmini(l)
         IF ( nm .eq. 0 ) CYCLE
         DO j = 1,nm
            i = lmord(lmini(l)+j-1)
            bphi(:,2*j-1) = REAL(a(i,3:n-2),KIND=GP)
            bphi(:,2*j)   = AIMAG(a(i,3:n-2))
         END DO
#if defined(GDOUBLE_PRECISION)
         CALL DGEMM('N','N',n,2*nm,n-4,1.0_GP,binv(:,:,l),n,bphi,n-4, &
                    0.0_GP,bw,n)
#else
         CALL SGEMM('N','N',n,2*nm,n-4,1.0_GP,binv(:,:,l),n,bphi,n-4, &
                    0.0_GP,bw,n)
#endif
         DO j = 1,nm
            i = lmord(lmini(l)+j-1)
            b(i,1:n) = CMPLX(bw(:,2*j-1),bw(:,2*j),KIND=GP)
         END DO
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
