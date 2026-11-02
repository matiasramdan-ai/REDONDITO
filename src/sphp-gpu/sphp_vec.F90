!=================================================================
! FFTP3D v3
! Parallel Fast Fourier Transform in 3D
!
! Performs parallel real-to-complex and complex-to-real FFTs 
! using MPI and the FFTW 3.x library in each node. You should use 
! the FFTPLANS and MPIVARS modules (see the file 'sphp_mod.F90') 
! in each program that calls any of the subroutines in this 
! file. Also, you must create plans for the parallel FFT using 
! the 'fftp3d_create_plan' subroutine, which creates in turn 
! derived data types for message passing using the subroutine 
! 'fftp3d_create_block'.
!
! When dealing with non-periodic data, an FC-Gram scheme
! is utilized to eliminate Gibbs phenomenon.
!
! Versión sphp, transformadas vectoriales: de un vector
! (v_r,v_theta,v_phi)(theta,phi,z) a sus coeficientes (Q,S,T)(lm,z)
! de SHTns, v = Q Y e_r + S r grad(Y) - T r x grad(Y), plano por
! plano, y al revés. Usan el plan, la FFT en z y los cambios de
! reparto de sphp_sca.F90; el reparto es el mismo: por planos
! (ksta:kend) en el espacio físico y por lm (lmsta:lmend) en el
! espectral.
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!
! 16 Feb 2004: Performs complex FFTs in place.
!  8 Jul 2004: itype pointers only used to store datatypes of 
!              blocks in the row and column each processor is.
!  9 Jul 2004: Transposition uses data cache blocking.
! 13 Feb 2007: Transposition uses strip mining (rreddy@psc.edu)
! 25 Aug 2009: Hybrid MPI/OpenMP support (D. Rosenberg & P. Mininni)
! 30 Aug 2009: SINGLE/DOUBLE precision (D. Rosenberg & P. Mininni)
!  3 Jan 2017: Anisotropic boxes (P. Mininni)
! 29 Aug 2019: FC-Gram (M. Fontana)
!    Dec 2020: Updated calls to MPI3 and MPI4 practices. (P. Mininni)
!
! References:
! Lyon, M.; SIAM J. Sci. Comp., 33(6) (2011). DOI:10.1137/11082436X
! Mininni PD, Gomez DO, Mahajan SM; Astrophys. J. 619, 1019 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Phys. Scripta T116, 123 (2005)
! Gomez DO, Mininni PD, Dmitruk P; Adv. Sp. Res. 35, 899 (2005)
!=================================================================
#include "fftw_wrappers.h"

!=================================================================
! MODULE fftp3d_vec_gpu
!
! Device path helpers of the vector transforms: execution of the
! vector transforms of SHTns on the device. The exchanges between
! tasks are the ones of the scalar transforms (module fftp3d_gpu in
! sphp_sca.F90), once per component.
!=================================================================
      MODULE fftp3d_vec_gpu
      IMPLICIT NONE
      CONTAINS

!*****************************************************************
      SUBROUTINE fftp3d_exec_vec_to_qst_dev(plan,vr,vt,vp,q,s,t)
!-----------------------------------------------------------------
!
! 2D vector spherical harmonic transform of the (theta,phi) planes
! on the device, plane by plane with SHTns: from (v_r,v_theta,v_phi)
! to (Q,S,T). SHTns does the FFT in phi in place, so the input data
! is destroyed
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN)   , INTENT(IN)            :: plan
      REAL(KIND=GP)   , INTENT(INOUT), TARGET :: vr(plan%nlat,plan%nphi,ksta:kend)
      REAL(KIND=GP)   , INTENT(INOUT), TARGET :: vt(plan%nlat,plan%nphi,ksta:kend)
      REAL(KIND=GP)   , INTENT(INOUT), TARGET :: vp(plan%nlat,plan%nphi,ksta:kend)
      COMPLEX(KIND=GP), INTENT(OUT)  , TARGET :: q(plan%nlm,ksta:kend)
      COMPLEX(KIND=GP), INTENT(OUT)  , TARGET :: s(plan%nlm,ksta:kend)
      COMPLEX(KIND=GP), INTENT(OUT)  , TARGET :: t(plan%nlm,ksta:kend)
      INTEGER :: k
!$omp target data use_device_addr(vr,vt,vp,q,s,t)
      DO k = ksta,kend
         CALL cu_spat_to_SHqst(plan%planrxy,vr(1,1,k),vt(1,1,k),vp(1,1,k), &
                               q(1,k),s(1,k),t(1,k),plan%lmax)
      END DO
!$omp end target data
      CALL device_sync()
      END SUBROUTINE fftp3d_exec_vec_to_qst_dev

!*****************************************************************
      SUBROUTINE fftp3d_exec_qst_to_vec_dev(plan,q,s,t,vr,vt,vp)
!-----------------------------------------------------------------
!
! Inverse 2D vector spherical harmonic transform of the (theta,phi)
! planes on the device, plane by plane with SHTns: from (Q,S,T) to
! (v_r,v_theta,v_phi). The input is INOUT because the interface of
! SHTns declares the coefficients so
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN)   , INTENT(IN)            :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET :: q(plan%nlm,ksta:kend)
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET :: s(plan%nlm,ksta:kend)
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET :: t(plan%nlm,ksta:kend)
      REAL(KIND=GP)   , INTENT(OUT)  , TARGET :: vr(plan%nlat,plan%nphi,ksta:kend)
      REAL(KIND=GP)   , INTENT(OUT)  , TARGET :: vt(plan%nlat,plan%nphi,ksta:kend)
      REAL(KIND=GP)   , INTENT(OUT)  , TARGET :: vp(plan%nlat,plan%nphi,ksta:kend)
      INTEGER :: k
!$omp target data use_device_addr(q,s,t,vr,vt,vp)
      DO k = ksta,kend
         CALL cu_SHqst_to_spat(plan%planrxy,q(1,k),s(1,k),t(1,k),       &
                               vr(1,1,k),vt(1,1,k),vp(1,1,k),plan%lmax)
      END DO
!$omp end target data
      CALL device_sync()
      END SUBROUTINE fftp3d_exec_qst_to_vec_dev

      END MODULE fftp3d_vec_gpu


!*****************************************************************
      SUBROUTINE fftp2d_vec_to_qst(plan,vr,vt,vp,q,s,t,comm)
!-----------------------------------------------------------------
! Computes the 2D vector spherical harmonic transform with SHTns in
! the first two indices, plane by plane: from the components of a
! vector, (v_r,v_theta,v_phi), to its coefficients (Q,S,T), with
! v = Q Y e_r + S r grad(Y) - T r x grad(Y). If the input arrays
! have shape (nlat,nphi,nz) the output arrays have shape (nlm,nz),
! con los lm en el orden de SHTns, after the transposition between
! tasks. The input data is destroyed during the computation.
! On the device path (gdev_active) the transform and the exchange
! run on the device.
!
! Parameters
!     plan     : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     vr,vt,vp : real input arrays v_r, v_theta, v_phi [INOUT]
!     q,s,t    : complex output arrays Q, S, T [OUT]
!     comm     : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fftplans
      USE fcgram
      USE gtimer
      USE gdevice
      USE fftp3d_gpu
      USE fftp3d_vec_gpu
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,plan%z%n)    :: q,s,t
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: vr,vt,vp

      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: k

      CALL GTStart(htot)

      IF (gdev_active) THEN
! gcarr1, gcarr2, gcarr3 and the exchange tables have the layout of
! the first plan created (see fftp3d_create_plan)
      IF (SIZE(gcarr1,1).ne.plan%planrc%nlm) &
         STOP 'fftp3d: the device path supports only plans with the layout of the first plan created'
!
! Device path: 2D vector transform of the (theta,phi) planes on the
! device
!
         CALL GTStart(hfft)
         CALL fftp3d_exec_vec_to_qst_dev(plan%planrc,vr,vt,vp,gcarr1,gcarr2,gcarr3)
         CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
!
! Exchange between tasks, with the blocks packed on the device, once
! per component
!
         CALL GTStart(hcom)
         CALL fftp3d_exchange_fwd(plan%planrc,gcarr1,q,comm)
         CALL fftp3d_exchange_fwd(plan%planrc,gcarr2,s,comm)
         CALL fftp3d_exchange_fwd(plan%planrc,gcarr3,t,comm)
         CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)

      ELSE
!
! Host path (sphp): 2D vector transform in each plane using the
! SHTns library. SHTns does the FFT in phi in place, so the input
! data is destroyed
!
      CALL GTStart(hfft)
      DO k = ksta,kend
         CALL spat_to_SHqst(plan%planrc%planrxy,vr(:,:,k),vt(:,:,k), &
                 vp(:,:,k),plan%planrc%carr1(:,k),plan%planrc%carr2(:,k), &
                 plan%planrc%carr3(:,k))
      END DO
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft); 
   
!
! Cambio de reparto entre las tareas, de mis planos a mis lm, una
! vez por componente, directamente en q, s y t
!
      CALL fftp3d_redist_k_to_lm(plan,plan%planrc%carr1,q,comm)
      CALL fftp3d_redist_k_to_lm(plan,plan%planrc%carr2,s,comm)
      CALL fftp3d_redist_k_to_lm(plan,plan%planrc%carr3,t,comm)

      ENDIF

      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_vec_to_qst

!*****************************************************************
      SUBROUTINE fftp2d_qst_to_vec(plan,q,s,t,vr,vt,vp,comm)
!-----------------------------------------------------------------
! Computes the inverse 2D vector spherical harmonic transform with
! SHTns, plane by plane, after the transposition between tasks:
! from the coefficients (Q,S,T) to the components of the vector,
! (v_r,v_theta,v_phi), with v = Q Y e_r + S r grad(Y) - T r x grad(Y).
! If the input arrays have shape (nlm,nz), con los lm en el orden de
! SHTns, the output arrays have shape (nlat,nphi,nz). The input data
! is not destroyed: los bloques se mandan directamente desde q, s y
! t, y SHqst_to_spat no cambia los coeficientes.
! On the device path (gdev_active) the exchange and the transform
! run on the device.
!
! Parameters
!     plan     : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     q,s,t    : complex input arrays Q, S, T [IN]
!     vr,vt,vp : real output arrays v_r, v_theta, v_phi [OUT]
!     comm     : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE fcgram
      USE gtimer
      USE gdevice
      USE fftp3d_gpu
      USE fftp3d_vec_gpu
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(lmsta:lmend,plan%z%n) :: q,s,t 
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: vr,vt,vp

      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: k

      CALL GTStart(htot)

      IF (gdev_active) THEN
! gcarr1, gcarr2, gcarr3 and the exchange tables have the layout of
! the first plan created (see fftp3d_create_plan)
      IF (SIZE(gcarr1,1).ne.plan%plancr%nlm) &
         STOP 'fftp3d: the device path supports only plans with the layout of the first plan created'
!
! Device path: exchange between tasks, with the blocks unpacked on
! the device, once per component
!
         CALL GTStart(hcom)
         CALL fftp3d_exchange_bwd(plan%plancr,q,gcarr1,comm)
         CALL fftp3d_exchange_bwd(plan%plancr,s,gcarr2,comm)
         CALL fftp3d_exchange_bwd(plan%plancr,t,gcarr3,comm)
         CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)
!
! 2D vector transform of the (theta,phi) planes on the device
!
         CALL GTStart(hfft)
         CALL fftp3d_exec_qst_to_vec_dev(plan%plancr,gcarr1,gcarr2,gcarr3,vr,vt,vp)
         CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      ELSE
!
! Host path (sphp). Cambio de reparto entre las tareas, de mis lm a
! mis planos, una vez por componente, directamente desde q, s y t
!
      CALL fftp3d_redist_lm_to_k(plan,q,plan%plancr%carr1,comm)
      CALL fftp3d_redist_lm_to_k(plan,s,plan%plancr%carr2,comm)
      CALL fftp3d_redist_lm_to_k(plan,t,plan%plancr%carr3,comm)
!
! 2D vector transform in each plane using the SHTns library, from
! all the lm of the plane, recibidos en el cambio de reparto
!
      CALL GTStart(hfft)
      DO k = ksta,kend
         CALL SHqst_to_spat(plan%plancr%planrxy,plan%plancr%carr1(:,k), &
                 plan%plancr%carr2(:,k), plan%plancr%carr3(:,k), &
                 vr(:,:,k),vt(:,:,k),vp(:,:,k))
      END DO
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      ENDIF

      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_qst_to_vec
