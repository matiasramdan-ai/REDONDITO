!======================================================================
! Condiciones de borde de la velocidad en las paredes de la cáscara.
! Tipos:
!   0 : no deslizamiento, W = dW/dr = Z = 0 (cuenta 4).
!
! Adaptado del boundary de SPECTER, 2020 Mauro Fontana. DF-UBA.
!======================================================================

!***********************************************************************
      FUNCTION v_parsebc(str) RESULT(bc)
!-----------------------------------------------------------------------
!  Function to map a given string value denoting a BC type into an 
!  integer for internal module usage.
!-----------------------------------------------------------------------
      USE mpivars
      USE commtypes
      IMPLICIT NONE

      CHARACTER(len=*), INTENT(IN) :: str
      INTEGER                      :: bc

      IF ( str .eq. 'noslip' ) THEN
         bc = 0
      ELSE
         IF ( myrank .eq. 0 ) PRINT*, "[ERROR] Unknown boundary ",& 
             "condition type", str, " Aborting..."
         CALL MPI_FINALIZE(ierr)
         STOP
      ENDIF
      RETURN
      ENDFUNCTION v_parsebc

!***************************************************************************
      SUBROUTINE v_setup(this,planfc,bckind)
!---------------------------------------------------------------------------
!  Guarda en el plan el tipo de condición de cada pared y arma lo que
! necesita v_imposebc. Por ahora el único tipo es el no deslizamiento,
! que en las dos paredes usa las matrices B_l de invDl (invDl_setup).
!  ARGUMENTS:
!    this   : 'this' instance of a BCPLAN. [INOUT]
!    planfc : FCPLAN that will be used for the boundary conditions. [INOUT]
!    bckind : tipo de condición en cada pared, rint y rext. [IN]
!---------------------------------------------------------------------------
      USE fcgram
 
      IMPLICIT NONE
 
      TYPE(BCPLAN), INTENT(INOUT)   :: this
      TYPE(FCPLAN), INTENT(INOUT)   :: planfc
      CHARACTER(len=15), INTENT(IN) :: bckind(2)

      this%bcint = v_parsebc(preprocess(bckind(1)))
      this%bcext = v_parsebc(preprocess(bckind(2)))
      CALL invDl_setup()

      END SUBROUTINE v_setup

!***********************************************************************
      SUBROUTINE v_imposebc(tor,phi,pol)
!-----------------------------------------------------------------------
!  Impone el no deslizamiento en las dos paredes al final de cada
! etapa del Runge-Kutta (cuenta 4, sección 6), sobre los campos en
! (lm,r):
!   1. Z = 0 en las paredes (dirichlet).
!   2. W = B_l^-1 P phi, con W = dW/dr = 0 en las paredes (invDl).
!   3. phi = D_l W (dlr), que cambia phi solo en los puntos 1, 2, N-1
!      y N, y deja en cero el modo l = 0.
! Corre en la CPU, porque invDl corre en la CPU.
!  ARGUMENTS:
!    tor : el potencial toroidal Z, en (lm,r). [INOUT]
!    phi : phi = D_l W, en (lm,r). [INOUT]
!    pol : el potencial poloidal W, en (lm,r). [OUT]
!-----------------------------------------------------------------------
      USE fprecision
      USE grid, ONLY: nz
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,nz) :: tor,phi
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz)   :: pol

      CALL dirichlet(tor)
      CALL invDl(phi,pol)
      CALL dlr(pol,phi)

      RETURN
      END SUBROUTINE v_imposebc

!***********************************************************************
      SUBROUTINE v_imposebc_fast(tor,phi,pol)
!-----------------------------------------------------------------------
!  Lo mismo que v_imposebc, pero en el paso 3 recalcula phi = D_l W
! solo en los puntos 1, 2, N-1 y N, con las filas de la matriz de
! d2/dr2 que guarda invDl_setup (d2bor), y pone en cero el modo l = 0.
! En los demás puntos deja el phi que entró, que ya es D_l W salvo el
! redondeo, porque esas filas de B_l son las de D_l. W sale completo,
! igual que en v_imposebc. Se ahorra las dos FFT en z de dlr.
!  ARGUMENTS:
!    tor : el potencial toroidal Z, en (lm,r). [INOUT]
!    phi : phi = D_l W, en (lm,r). [INOUT]
!    pol : el potencial poloidal W, en (lm,r). [OUT]
!-----------------------------------------------------------------------
      USE fprecision
      USE grid, ONLY: nz,Cz,z
      USE kes, ONLY: ell
      USE mpivars
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,nz) :: tor,phi
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,nz)   :: pol
      INTEGER :: n,lm,k,i,fila(4)

      n = nz-Cz
      fila = (/ 1,2,n-1,n /)
      CALL dirichlet(tor)
      CALL invDl(phi,pol)
      DO lm = lmsta,lmend
         IF ( ell(lm) .eq. 0 ) THEN
            phi(lm,:) = 0.0_GP
            CYCLE
         ENDIF
         DO k = 1,4
            i = fila(k)
            phi(lm,i) = SUM(d2bor(k,:)*pol(lm,1:n)) - &
                        ell(lm)*(ell(lm)+1)*pol(lm,i)/z(i)**2
         END DO
      END DO

      RETURN
      END SUBROUTINE v_imposebc_fast
