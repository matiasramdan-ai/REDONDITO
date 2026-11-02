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
! Versión sphp: la FFT 2D de FFTW en (x,y) está reemplazada por la
! transformada en armónicos esféricos de SHTns en (theta,phi),
! plano por plano; la FFT en z sigue siendo la de FFTW, con
! FC-Gram. El espacio físico es f(theta,phi,z), repartido entre
! las tareas por planos (ksta:kend), y el espectral f(lm,kz),
! repartido por lm (lmsta:lmend). Todo sale normalizado: la ida y
! vuelta es la identidad. Los nombres de las rutinas y de los
! módulos son todavía los de SPECTER.
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
#include "gfft_defs.h"

!=================================================================
! MODULE fftp3d_gpu
!
! Device path helpers of FFTP-GPU: execution of the device
! transforms (SHTns and hipFFT), packing of the exchange blocks and
! the exchanges between tasks. The array dummies are explicit-shape
! with TARGET: inside a "target data use_device_addr" region they
! refer to the device copies of the arrays, so C_LOC, or an element
! passed to an assumed-size dummy, yields device addresses for the
! transform libraries and for MPI.
!=================================================================
      MODULE fftp3d_gpu
      IMPLICIT NONE
      CONTAINS

!*****************************************************************
      SUBROUTINE fftp3d_exec_r2c_dev(plan,in,out)
!-----------------------------------------------------------------
!
! 2D spherical harmonic transform of the (theta,phi) planes on the
! device, plane by plane with SHTns. SHTns does the FFT in phi in
! place, so the input data is destroyed
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN)   , INTENT(IN)            :: plan
      REAL(KIND=GP)   , INTENT(INOUT), TARGET :: in (plan%nlat,plan%nphi,ksta:kend)
      COMPLEX(KIND=GP), INTENT(OUT)  , TARGET :: out(plan%nlm,ksta:kend)
      INTEGER :: k
!$omp target data use_device_addr(in,out)
      DO k = ksta,kend
         CALL cu_spat_to_SH(plan%planrxy,in(1,1,k),out(1,k),plan%lmax)
      END DO
!$omp end target data
      CALL device_sync()
      END SUBROUTINE fftp3d_exec_r2c_dev

!*****************************************************************
      SUBROUTINE fftp3d_exec_c2r_dev(plan,in,out)
!-----------------------------------------------------------------
!
! Inverse 2D spherical harmonic transform of the (theta,phi) planes
! on the device, plane by plane with SHTns. The input is INOUT
! because the interface of SHTns declares the coefficients so
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN)   , INTENT(IN)            :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET :: in (plan%nlm,ksta:kend)
      REAL(KIND=GP)   , INTENT(OUT)  , TARGET :: out(plan%nlat,plan%nphi,ksta:kend)
      INTEGER :: k
!$omp target data use_device_addr(in,out)
      DO k = ksta,kend
         CALL cu_SH_to_spat(plan%planrxy,in(1,k),out(1,1,k),plan%lmax)
      END DO
!$omp end target data
      CALL device_sync()
      END SUBROUTINE fftp3d_exec_c2r_dev

!*****************************************************************
      SUBROUTINE fftp3d_exec_c2c_dev(plan,a,dir)
!-----------------------------------------------------------------
!
! 1D complex transform along z, in place, on the device
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN)   , INTENT(IN)            :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), TARGET :: a(lmsta:lmend,plan%nz)
      INTEGER, INTENT(IN) :: dir
      INTEGER(C_INT) :: ier
!$omp target data use_device_addr(a)
      ier = gfft_exec_c2c(plan%dplancz,C_LOC(a),C_LOC(a),INT(dir,C_INT))
!$omp end target data
      IF (ier.ne.0) PRINT *,'fftp3d_exec_c2c_dev: error ',ier
      CALL device_sync()
      END SUBROUTINE fftp3d_exec_c2c_dev

!*****************************************************************
      SUBROUTINE fftp3d_pack_dev(carr,sbuf,nlm,ntot,l0,lr,o)
!-----------------------------------------------------------------
!
! Packs the block of carr with lm in [l0,l0+lr-1] contiguously at
! sbuf(o+1:), on the device
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      IMPLICIT NONE
      INTEGER         , INTENT(IN)    :: nlm,ntot,l0,lr,o
      COMPLEX(KIND=GP), INTENT(IN)    :: carr(nlm,ksta:kend)
      COMPLEX(KIND=GP), INTENT(INOUT) :: sbuf(ntot)
      INTEGER :: l,k
!$omp target teams distribute parallel do collapse(2)
      DO k = ksta,kend
         DO l = 1,lr
            sbuf(o+l+lr*(k-ksta)) = carr(l0+l-1,k)
         END DO
      END DO
      END SUBROUTINE fftp3d_pack_dev

!*****************************************************************
      SUBROUTINE fftp3d_unpack_dev(sbuf,carr,nlm,ntot,l0,lr,o)
!-----------------------------------------------------------------
!
! Unpacks the block at sbuf(o+1:) into carr, lm in [l0,l0+lr-1]
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      IMPLICIT NONE
      INTEGER         , INTENT(IN)    :: nlm,ntot,l0,lr,o
      COMPLEX(KIND=GP), INTENT(IN)    :: sbuf(ntot)
      COMPLEX(KIND=GP), INTENT(INOUT) :: carr(nlm,ksta:kend)
      INTEGER :: l,k
!$omp target teams distribute parallel do collapse(2)
      DO k = ksta,kend
         DO l = 1,lr
            carr(l0+l-1,k) = sbuf(o+l+lr*(k-ksta))
         END DO
      END DO
      END SUBROUTINE fftp3d_unpack_dev

!*****************************************************************
      SUBROUTINE fftp3d_copy_dev(dst,src,n)
!-----------------------------------------------------------------
!
! Contiguous copy of n complex numbers on the device (the block a
! task exchanges with itself never goes through MPI)
!-----------------------------------------------------------------
      USE fprecision
      IMPLICIT NONE
      INTEGER         , INTENT(IN)  :: n
      COMPLEX(KIND=GP), INTENT(OUT) :: dst(n)
      COMPLEX(KIND=GP), INTENT(IN)  :: src(n)
      INTEGER :: i
!$omp target teams distribute parallel do
      DO i = 1,n
         dst(i) = src(i)
      END DO
      END SUBROUTINE fftp3d_copy_dev

!*****************************************************************
      SUBROUTINE fftp3d_exch_contig(sendbuf,nsend,recvbuf,nrecv,so,sc,ro,rc,comm)
!-----------------------------------------------------------------
!
! Exchange of contiguous blocks between tasks with device buffers:
! the block for task r starts at sendbuf(so(r)+1) with sc(r)
! elements, the block from task r is received at recvbuf(ro(r)+1)
! with rc(r) elements. The diagonal block is copied on the device.
! The buffers are entered in a use_device_addr region and handed to
! MPI, which must be GPU-aware.
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      INTEGER         , INTENT(IN)          :: nsend,nrecv
      COMPLEX(KIND=GP), INTENT(IN) , TARGET :: sendbuf(nsend)
      COMPLEX(KIND=GP), INTENT(OUT), TARGET :: recvbuf(nrecv)
      INTEGER, INTENT(IN) :: so(0:nprocs-1),sc(0:nprocs-1)
      INTEGER, INTENT(IN) :: ro(0:nprocs-1),rc(0:nprocs-1)
      TYPE(MPI_Comm)  , INTENT(IN)          :: comm

      CALL fftp3d_copy_dev(recvbuf(ro(myrank)+1),sendbuf(so(myrank)+1),sc(myrank))
!$omp target data use_device_addr(sendbuf,recvbuf)
      CALL fftp3d_exch_contig_do(sendbuf,recvbuf,so,sc,ro,rc,comm)
!$omp end target data
      END SUBROUTINE fftp3d_exch_contig

!*****************************************************************
      SUBROUTINE fftp3d_exch_contig_do(sendbuf,recvbuf,so,sc,ro,rc,comm)
!-----------------------------------------------------------------
!
! The MPI calls of fftp3d_exch_contig, receiving device addresses
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      COMPLEX(KIND=GP), INTENT(IN)  :: sendbuf(*)
      COMPLEX(KIND=GP), INTENT(OUT) :: recvbuf(*)
      INTEGER, INTENT(IN) :: so(0:nprocs-1),sc(0:nprocs-1)
      INTEGER, INTENT(IN) :: ro(0:nprocs-1),rc(0:nprocs-1)
      TYPE(MPI_Comm)  , INTENT(IN)  :: comm
      TYPE(MPI_Request) :: ireq(2*nprocs)
      INTEGER :: irank,isendTo,igetFrom

! All receives and sends are posted before waiting, so that the
! transfers to and from the different tasks proceed at the same time
      DO irank = 1,nprocs-1
         igetFrom = myrank - irank
         IF (igetFrom.lt.0) igetFrom = igetFrom + nprocs
         CALL fftp3d_irecv(recvbuf(ro(igetFrom)+1),rc(igetFrom),igetFrom,comm,ireq(irank))
      END DO
      DO irank = 1,nprocs-1
         isendTo = myrank + irank
         IF (isendTo.ge.nprocs) isendTo = isendTo - nprocs
         CALL fftp3d_isend(sendbuf(so(isendTo)+1),sc(isendTo),isendTo,comm,ireq(nprocs-1+irank))
      END DO
      IF (nprocs.gt.1) CALL MPI_WAITALL(2*(nprocs-1),ireq(1:2*(nprocs-1)),MPI_STATUSES_IGNORE,ierr)
      END SUBROUTINE fftp3d_exch_contig_do

!*****************************************************************
      SUBROUTINE fftp3d_irecv(buf,n,src,comm,req)
!-----------------------------------------------------------------
!
! Receives n complex numbers starting at buf (sequence association)
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n,src
      COMPLEX(KIND=GP) , INTENT(OUT) :: buf(n)
      TYPE(MPI_Comm)   , INTENT(IN)  :: comm
      TYPE(MPI_Request), INTENT(OUT) :: req
      CALL MPI_IRECV(buf,n,GC_COMPLEX,src,1,comm,req,ierr)
      END SUBROUTINE fftp3d_irecv

!*****************************************************************
      SUBROUTINE fftp3d_isend(buf,n,dst,comm,req)
!-----------------------------------------------------------------
!
! Sends n complex numbers starting at buf
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n,dst
      COMPLEX(KIND=GP) , INTENT(IN)  :: buf(n)
      TYPE(MPI_Comm)   , INTENT(IN)  :: comm
      TYPE(MPI_Request), INTENT(OUT) :: req
      CALL MPI_ISEND(buf,n,GC_COMPLEX,dst,1,comm,req,ierr)
      END SUBROUTINE fftp3d_isend

!*****************************************************************
      SUBROUTINE fftp3d_exch_types_host(sendbuf,nsend,recvbuf,nrecv,stype,rtype,comm)
!-----------------------------------------------------------------
!
! Exchange with derived data types on the host copies (staged mode)
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nsend,nrecv
      COMPLEX(KIND=GP)  , INTENT(IN)  :: sendbuf(nsend)
      COMPLEX(KIND=GP)  , INTENT(OUT) :: recvbuf(nrecv)
      TYPE(MPI_Datatype), INTENT(IN)  :: stype(0:nprocs-1),rtype(0:nprocs-1)
      TYPE(MPI_Comm)    , INTENT(IN)  :: comm
      TYPE(MPI_Request) :: ireq1,ireq2
      TYPE(MPI_Status)  :: istatus
      INTEGER :: irank,isendTo,igetFrom

      DO irank = 0,nprocs-1
         isendTo = myrank + irank
         IF (isendTo.ge.nprocs) isendTo = isendTo - nprocs
         igetFrom = myrank - irank
         IF (igetFrom.lt.0) igetFrom = igetFrom + nprocs
         CALL MPI_IRECV(recvbuf,1,rtype(igetFrom),igetFrom,1,comm,ireq2,ierr)
         CALL MPI_ISEND(sendbuf,1,stype(isendTo) ,isendTo ,1,comm,ireq1,ierr)
         CALL MPI_WAIT(ireq1,istatus,ierr)
         CALL MPI_WAIT(ireq2,istatus,ierr)
      END DO
      END SUBROUTINE fftp3d_exch_types_host

!*****************************************************************
      SUBROUTINE fftp3d_exchange_fwd(plan,carr,out,comm)
!-----------------------------------------------------------------
!
! Exchange carr -> out of the device path, from the planes of each
! task to the lm of each task. The blocks received are contiguous
! in out(lmsta:lmend,nz), so only the blocks sent are packed
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN)   , INTENT(IN)          :: plan
      COMPLEX(KIND=GP), INTENT(IN) , TARGET :: carr(plan%nlm,ksta:kend)
      COMPLEX(KIND=GP), INTENT(OUT), TARGET :: out(lmsta:lmend,plan%nz)
      TYPE(MPI_Comm)  , INTENT(IN)          :: comm
      INTEGER :: r,ntot

      ntot = SIZE(gsbuf)
      IF (gexchange.eq.1) THEN
         DO r = 0,nprocs-1
            CALL fftp3d_pack_dev(carr,gsbuf,plan%nlm,ntot,glmoff(r),glmlen(r),gsoff(r))
         END DO
         CALL fftp3d_exch_contig(gsbuf,ntot,out,SIZE(out),gsoff,gscnt,groff,grcnt,comm)
      ELSE
         CALL gdev_update_from(C_LOC(carr),INT(SIZE(carr),C_SIZE_T)*INT(STORAGE_SIZE(carr)/8,C_SIZE_T))
         CALL fftp3d_exch_types_host(carr,SIZE(carr),out,SIZE(out),plan%itype1,plan%itype2,comm)
         CALL gdev_update_to(C_LOC(out),INT(SIZE(out),C_SIZE_T)*INT(STORAGE_SIZE(out)/8,C_SIZE_T))
      ENDIF
      END SUBROUTINE fftp3d_exchange_fwd

!*****************************************************************
      SUBROUTINE fftp3d_exchange_bwd(plan,in,carr,comm)
!-----------------------------------------------------------------
!
! Exchange in -> carr of the device path, from the lm of each task
! to the planes of each task. The blocks sent are contiguous in
! in(lmsta:lmend,nz), so only the blocks received are unpacked
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE gdevice
      IMPLICIT NONE
      TYPE(FFTPLAN)   , INTENT(IN)          :: plan
      COMPLEX(KIND=GP), INTENT(IN) , TARGET :: in(lmsta:lmend,plan%nz)
      COMPLEX(KIND=GP), INTENT(OUT), TARGET :: carr(plan%nlm,ksta:kend)
      TYPE(MPI_Comm)  , INTENT(IN)          :: comm
      INTEGER :: r,ntot

      ntot = SIZE(gsbuf)
      IF (gexchange.eq.1) THEN
         CALL fftp3d_exch_contig(in,SIZE(in),gsbuf,ntot,groff,grcnt,gsoff,gscnt,comm)
         DO r = 0,nprocs-1
            CALL fftp3d_unpack_dev(gsbuf,carr,plan%nlm,ntot,glmoff(r),glmlen(r),gsoff(r))
         END DO
      ELSE
         CALL gdev_update_from(C_LOC(in),INT(SIZE(in),C_SIZE_T)*INT(STORAGE_SIZE(in)/8,C_SIZE_T))
         CALL fftp3d_exch_types_host(in,SIZE(in),carr,SIZE(carr),plan%itype2,plan%itype1,comm)
         CALL gdev_update_to(C_LOC(carr),INT(SIZE(carr),C_SIZE_T)*INT(STORAGE_SIZE(carr)/8,C_SIZE_T))
      ENDIF
      END SUBROUTINE fftp3d_exchange_bwd

!*****************************************************************
      SUBROUTINE fftp3d_cont_z_dev(nz,C,d,dir,a)
!-----------------------------------------------------------------
!
! FC-Gram continuation along z on the device, with the same
! operations as the host loop in fftp1d_real_to_complex_z
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      IMPLICIT NONE
      INTEGER         , INTENT(IN)    :: nz,C,d
      DOUBLE PRECISION, INTENT(IN)    :: dir(C,d)
      COMPLEX(KIND=GP), INTENT(INOUT) :: a(lmsta:lmend,nz)
      INTEGER :: l,ii,jj
!$omp target teams distribute parallel do collapse(2)
      DO ii = 1,C
         DO l = lmsta,lmend
            a(l,nz-C+ii) = dir(ii,1)*a(l,nz-C-d+1) +              &
                           dir(C-ii+1,1)*a(l,d)
            DO jj = 2,d
               a(l,nz-C+ii) = a(l,nz-C+ii) +                      &
                              dir(ii,jj)*a(l,nz-C-d+jj) +         &
                              dir(C-ii+1,jj)*a(l,d-jj+1)
            END DO
         END DO
      END DO
      END SUBROUTINE fftp3d_cont_z_dev

!*****************************************************************
      SUBROUTINE fftp3d_norm_z_dev(nz,a)
!-----------------------------------------------------------------
!
! Normalization of the forward transform along z on the device, as
! in the host path. It multiplies by tmp = 1/nz: the division of a
! complex number has no device implementation in amdflang (the
! link fails with an undefined __divdc3)
!-----------------------------------------------------------------
      USE fprecision
      USE mpivars
      IMPLICIT NONE
      INTEGER         , INTENT(IN)    :: nz
      COMPLEX(KIND=GP), INTENT(INOUT) :: a(lmsta:lmend,nz)
      REAL(KIND=GP)   :: tmp
      INTEGER :: l,k
      tmp = 1.0_GP/REAL(nz,KIND=GP)
!$omp target teams distribute parallel do collapse(2)
      DO k = 1,nz
         DO l = lmsta,lmend
            a(l,k) = a(l,k)*tmp
         END DO
      END DO
      END SUBROUTINE fftp3d_norm_z_dev

      END MODULE fftp3d_gpu

!=================================================================
! General FFTP-GPU routines
!
! Initialization, plan creator, destructor, block creator,
! complex-to-real and real-to-complex parallel transforms for host
! and device paths.
!=================================================================

!*****************************************************************
      SUBROUTINE fftp3d_init_threads(err)
!-----------------------------------------------------------------
! Initializes FFTW threads.
!
! Parameters
!     err : if zero, the initialization failed
!-----------------------------------------------------------------

!$    USE threads
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: err

!$    CALL GFFTW_INIT_THREADS(err)
      IF (err.eq.0) PRINT *,'FFTP threads initialization failed!'

      RETURN
      END SUBROUTINE fftp3d_init_threads

!*****************************************************************
      SUBROUTINE fftp3d_create_plan(plan,n,fftdir,flags,lmax)
!-----------------------------------------------------------------
! Creates the plans in each node: the SHTns plan for the 2D
! transform in (theta,phi), also enabled for the device, the host
! (FFTW) and device (hipFFT) plans for the FFT in z, and the
! shared device-resident buffers of the exchange. The device of
! the task must be selected before (device_init).
!
! Parameters
!     plan   : contains the parallel 3D plan [OUT]
!     n      : the size of the dimensions of the input array,
!              (nlat,nphi,nz) [IN]
!     fftdir : the direction of the Fourier transform [IN]
!              FFTW_FORWARD or FFTW_REAL_TO_COMPLEX (-1)
!              FFTW_BACKWARD or FFTW_COMPLEX_TO_REAL (+1)
!     flags  : flags for the FFTW [IN], kept for compatibility with
!              sphp and ignored: the host plan is always created
!              with FFTW_ESTIMATE, as in fftp-gpu. In offload builds
!              it is used only outside the time step, and measured
!              plans can change from run to run, giving outputs that
!              differ at roundoff level between identical runs. The
!              device plans take no flags.
!     lmax   : grado máximo de los armónicos esféricos [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE mpivars
      USE fftplans
      USE gfft
      USE gdevice
!$    USE threads
      USE gtimer
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3)
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      INTEGER, INTENT(IN) :: lmax
      TYPE(FFTPLAN), INTENT(OUT) :: plan

      INTEGER(C_INT)   :: n1(1),nk,nb,ier
      INTEGER          :: r,ls,le,ks,ke,nlm
      CHARACTER(len=16):: senv

! Host plans, as in sphp but always with FFTW_ESTIMATE (see above)
      ALLOCATE ( plan%ccarr(lmsta:lmend,n(3))    )
      ALLOCATE ( plan%carr1((lmax+1)*(lmax+2)/2,ksta:kend) )
      ALLOCATE ( plan%carr2((lmax+1)*(lmax+2)/2,ksta:kend) )
      ALLOCATE ( plan%carr3((lmax+1)*(lmax+2)/2,ksta:kend) )
!$    CALL GFFTW_PLAN_WITH_NTHREADS(nth)

! Plan de SHTns para la transformada 2D, el mismo en las dos
! direcciones: armónicos ortonormales sin la fase de
! Condon-Shortley, como en MagIC, en una grilla de Gauss con theta
! contiguo, de modo que cada plano es f(theta,phi). La
! optimización polar es 1e-14 y no la 1e-10 de MagIC: con 1e-10 la
! síntesis de modos de m alto pierde hasta 1e-11, y con 1e-14 queda
! en el redondeo y cuesta lo mismo
      plan%planrxy = shtns_create(lmax,lmax,1, &
                                  SHT_ORTHONORMAL+SHT_NO_CS_PHASE)

      CALL shtns_set_grid(plan%planrxy, &
                          SHT_GAUSS+SHT_THETA_CONTIGUOUS,1.0D-14, &
                          n(1),n(2))
      
      CALL GFFTW_PLAN_MANY_DFT(plan%plancz,1,n(3),lmend-lmsta+1,  &
                         plan%ccarr,n(3),lmend-lmsta+1,1,   &
                         plan%ccarr,n(3),lmend-lmsta+1,1,   &
                         fftdir,FFTW_ESTIMATE)
      plan%nlat = n(1)
      plan%nphi = n(2)
      plan%nz   = n(3)
      plan%nlm  = (lmax+1)*(lmax+2)/2
      plan%lmax = lmax
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_block(n,(lmax+1)*(lmax+2)/2,nprocs,myrank, &
              plan%itype1,plan%itype2)

! Device plans. The 2D transform is done by the SHTns plan itself,
! enabled for the device of the task (cushtns_init_gpu), which puts
! its tables on the device for the cu_* transforms. hipFFT takes the
! dimensions in C order: the 1D plan transforms along z, the second
! dimension of ccarr(lmsta:lmend,nz), in place, with stride nb
! between elements and 1 between transforms, as the FFTW plan.
      nlm   = (lmax+1)*(lmax+2)/2
      nk    = kend-ksta+1
      nb    = lmend-lmsta+1
      n1    = [INT(n(3),C_INT)]
      ier = cushtns_init_gpu(plan%planrxy)
      IF (ier.lt.0) THEN
         PRINT *,'fftp3d_create_plan: 2D device plan failed, error ',ier
         STOP
      ENDIF
      ier = gfft_plan_many(plan%dplancz,1_C_INT,n1,n1,nb,1_C_INT,       &
            n1,nb,1_C_INT,INT(GFFT_TYPE_C2C,C_INT),nb)
      IF (ier.ne.0) THEN
         PRINT *,'fftp3d_create_plan: 1D device plan failed, error ',ier
         STOP
      ENDIF

! Shared exchange buffers and block tables (once)
      IF (.not.gbuffers_ready) THEN
         ALLOCATE( gcarr1(nlm,ksta:kend) )
         ALLOCATE( gcarr2(nlm,ksta:kend) )
         ALLOCATE( gcarr3(nlm,ksta:kend) )
         ALLOCATE( gsbuf(nlm*nk)         )
         gcarr1 = 0.0_GP; gcarr2 = 0.0_GP; gcarr3 = 0.0_GP; gsbuf = 0.0_GP
         CALL gdev_alloc(C_LOC(gcarr1),INT(SIZE(gcarr1),C_SIZE_T)*INT(STORAGE_SIZE(gcarr1)/8,C_SIZE_T))
         CALL gdev_alloc(C_LOC(gcarr2),INT(SIZE(gcarr2),C_SIZE_T)*INT(STORAGE_SIZE(gcarr2)/8,C_SIZE_T))
         CALL gdev_alloc(C_LOC(gcarr3),INT(SIZE(gcarr3),C_SIZE_T)*INT(STORAGE_SIZE(gcarr3)/8,C_SIZE_T))
         CALL gdev_alloc(C_LOC(gsbuf) ,INT(SIZE(gsbuf) ,C_SIZE_T)*INT(STORAGE_SIZE(gsbuf)/8 ,C_SIZE_T))
         ALLOCATE( glmoff(0:nprocs-1), glmlen(0:nprocs-1) )
         ALLOCATE( gkoff(0:nprocs-1) , gklen(0:nprocs-1)  )
         ALLOCATE( gsoff(0:nprocs-1) , gscnt(0:nprocs-1)  )
         ALLOCATE( groff(0:nprocs-1) , grcnt(0:nprocs-1)  )
         DO r = 0,nprocs-1
            CALL range(1,nlm ,nprocs,r,ls,le)
            CALL range(1,n(3),nprocs,r,ks,ke)
            glmoff(r) = ls; glmlen(r) = le-ls+1
            gkoff(r)  = ks; gklen(r)  = ke-ks+1
            gscnt(r)  = glmlen(r)*nk
            grcnt(r)  = (lmend-lmsta+1)*gklen(r)
            groff(r)  = (lmend-lmsta+1)*(ks-1)
            IF (r.eq.0) THEN
               gsoff(r) = 0
            ELSE
               gsoff(r) = gsoff(r-1)+gscnt(r-1)
            ENDIF
         END DO
         CALL GET_ENVIRONMENT_VARIABLE('GHOST_GPU_EXCHANGE',senv)
         IF (senv(1:6).eq.'staged') THEN
            gexchange = 2
            IF (myrank.eq.0) PRINT *,'FFTP-GPU: exchange staged through the host'
         ELSE
            gexchange = 1
         ENDIF
         gbuffers_ready = .TRUE.
      ENDIF

      CALL GTStart(hcom,GT_WTIME)
      CALL GTStart(hfft,GT_WTIME)
      CALL GTStart(htra,GT_WTIME)
      CALL GTStart(hcont,GT_WTIME)
      CALL GTStart(htot,GT_WTIME)

      RETURN
      END SUBROUTINE fftp3d_create_plan

!*****************************************************************
      SUBROUTINE fftp3d_destroy_plan(plan)
!-----------------------------------------------------------------
! Destroys the SHTns, FFTW and hipFFT plans in each node, and the
! shared device-resident buffers of the exchange.
!
! Parameters
!     plan : the parallel 3D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE gfft
      USE gdevice
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan
      INTEGER(C_INT) :: ier

! shtns_destroy libera también las tablas de SHTns en la GPU
      CALL shtns_destroy(plan%planrxy)
      CALL GFFTW_DESTROY_PLAN(plan%plancz)
      ier = gfft_destroy(plan%dplancz)

      DEALLOCATE( plan%ccarr  )
      DEALLOCATE( plan%carr1  )
      DEALLOCATE( plan%carr2  )
      DEALLOCATE( plan%carr3  )
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )

      IF (gbuffers_ready) THEN
         CALL gdev_free(C_LOC(gcarr1))
         CALL gdev_free(C_LOC(gcarr2))
         CALL gdev_free(C_LOC(gcarr3))
         CALL gdev_free(C_LOC(gsbuf))
         DEALLOCATE( gcarr1, gcarr2, gcarr3, gsbuf )
         DEALLOCATE( glmoff, glmlen, gkoff, gklen )
         DEALLOCATE( gsoff, gscnt, groff, grcnt )
         gbuffers_ready = .FALSE.
      ENDIF

      CALL GTFree(hcom)
      CALL GTFree(hfft)
      CALL GTFree(htra)
      CALL GTFree(hcont)
      CALL GTFree(htot)

      RETURN
      END SUBROUTINE fftp3d_destroy_plan

!*****************************************************************
      SUBROUTINE fftp3d_create_block(n,nlm,nprocs,myrank,itype1,itype2)
!-----------------------------------------------------------------
! Defines derived data types for sending and receiving 
! blocks of the 3D matrix between processors. The data 
! types are used to transpose the matrix during the FFT.
!
! Versión sphp: la matriz es (lm,z). Antes de la transposición,
! cada tarea tiene los nlm coeficientes de sus planos,
! carr1(nlm,ksta:kend); después, todos los planos de sus lm,
! out(lmsta:lmend,nz). Los bloques se arman con block3d, con la
! dimensión del medio de largo 1.
!
! Parameters
!     n      : the size of the dimensions of the input array [IN]
!     nlm    : cantidad de coeficientes (l,m) por plano [IN]
!     nprocs : the number of processors [IN]
!     myrank : the rank of the processor [IN]
!     itype1 : contains a derived data type for sending [OUT]
!     itype2 : contains a derived data type for receiving [OUT]
!-----------------------------------------------------------------

      USE commtypes
      IMPLICIT NONE

      TYPE(MPI_Datatype), INTENT(OUT), DIMENSION(0:nprocs-1) :: itype1,itype2
      INTEGER, INTENT(IN) :: n(3),nlm,nprocs
      INTEGER, INTENT(IN) :: myrank

      INTEGER :: lmsta,lmend
      INTEGER :: ksta,kend
      INTEGER :: irank,krank
      TYPE(MPI_Datatype) :: itemp1,itemp2

! itype1(irank): en carr1(nlm,ksta:kend), los lm de la tarea irank
! en todos mis planos
      CALL range(1,n(3),nprocs,myrank,ksta,kend)
      DO irank = 0,nprocs-1
         CALL range(1,nlm,nprocs,irank,lmsta,lmend)
         CALL block3d(1,nlm,1,1,ksta,lmsta,lmend,1,1, &
                     ksta,kend,GC_COMPLEX,itemp1)
         itype1(irank) = itemp1
      END DO
! itype2(krank): en out(lmsta:lmend,nz), mis lm en los planos de la
! tarea krank
      CALL range(1,nlm,nprocs,myrank,lmsta,lmend)
      DO krank = 0,nprocs-1
         CALL range(1,n(3),nprocs,krank,ksta,kend)
         CALL block3d(lmsta,lmend,1,1,1,lmsta,lmend,1,1,     &
                     ksta,kend,GC_COMPLEX,itemp2)
         itype2(krank) = itemp2
      END DO

      RETURN
      END SUBROUTINE fftp3d_create_block

!*****************************************************************
      SUBROUTINE fftp3d_real_to_complex(plan,in,out,comm)
!-----------------------------------------------------------------
! Computes the 3D transform in parallel: the 2D spherical harmonic
! transform with SHTns, the transposition between tasks and the
! FFT in z. If the input is non-periodic in z FC-Gram periodic
! extensions are appropriately computed. If the input is
! f(theta,phi,z) the output is f(lm,kz), normalized. The input
! data is destroyed during the computation.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     in   : real input array f(theta,phi,z) [INOUT]
!     out  : complex output array f(lm,kz) [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fcgram
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,plan%z%n) :: out
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: in

      TYPE(MPI_Comm)   , INTENT(IN)            :: comm

      CALL fftp2d_real_to_complex_xy(plan,in,out,comm)
      CALL fftp1d_real_to_complex_z(plan,out,comm)

      RETURN
      END SUBROUTINE fftp3d_real_to_complex

!*****************************************************************
      SUBROUTINE fftp2d_real_to_complex_xy(plan,in,out,comm)
!-----------------------------------------------------------------
! Computes the 2D spherical harmonic transform with SHTns in the
! first two indices, plane by plane. If the input has shape
! (nlat,nphi,nz) the output's shape is (nlm,nz), con los lm en el
! orden de SHTns, after the transposition between tasks. Con los
! armónicos ortonormales no hay que normalizar. The input data is
! destroyed during the computation.
! On the device path (gdev_active) the transform and the exchange
! run on the device.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     in   : real input array f(theta,phi,z) [INOUT]
!     out  : complex output array f(lm,z) [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fftplans
      USE fcgram
      USE gtimer
      USE gdevice
      USE fftp3d_gpu
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,plan%z%n)    :: out
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: IN

      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: k

      CALL GTStart(htot)
      IF (gdev_active) THEN
! gcarr1 and the exchange tables have the layout of the first plan
! created (see fftp3d_create_plan)
      IF (SIZE(gcarr1,1).ne.plan%planrc%nlm) &
         STOP 'fftp3d: the device path supports only plans with the layout of the first plan created'
!
! Device path: 2D transform of the (theta,phi) planes on the device
!
         CALL GTStart(hfft)
         CALL fftp3d_exec_r2c_dev(plan%planrc,in,gcarr1)
         CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
!
! Exchange between tasks, with the blocks packed on the device
!
         CALL GTStart(hcom)
         CALL fftp3d_exchange_fwd(plan%planrc,gcarr1,out,comm)
         CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)

      ELSE
!
! Host path (sphp): 2D transform in each plane using the SHTns
! library. SHTns does the FFT in phi in place, so the input data is
! destroyed
!
      CALL GTStart(hfft)
      DO k = ksta,kend
         CALL spat_to_SH(plan%planrc%planrxy,in(:,:,k),plan%planrc%carr1(:,k))
      END DO
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft); 
   
!
! Transposes the result between nodes using strip mining when
! nstrip>1 (rreddy@psc.edu). Cada tarea manda a las otras los lm
! que les tocan de sus planos y recibe sus lm de los planos de las
! otras, directamente en out
!
      CALL fftp3d_redist_k_to_lm(plan,plan%planrc%carr1,out,comm)
!
! La transposición local de SPECTER no hace falta: los bloques
! llegan en el orden (lm,z) de out
      ENDIF
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_real_to_complex_xy

!*****************************************************************
      SUBROUTINE fftp1d_real_to_complex_z(plan,inout,comm)
!-----------------------------------------------------------------
! Computes the 1D forward FFT along the second index of the input
! array, (lm,z). In the case the input array is non-periodic, a
! suitable FC-Gram periodic extension is computed. The operation is
! performed inplace. El resultado se divide por nz, como el
! análisis de SHTns: los coeficientes son los del campo y la ida y
! vuelta es la identidad, sin normalizar después. On the device
! path (gdev_active) the continuation, the FFT and the
! normalization run on the device.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     in   : contains f(lm,z) as input and returns f(lm,kz) [INOUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE gtimer
      USE fcgram
      USE gdevice
      USE fftp3d_gpu
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,plan%z%n) :: inout
      TYPE(MPI_Comm)   , INTENT(IN)            :: comm

      INTEGER :: n,C,d
      INTEGER :: ii,jj

      n = plan%z%n
      C = plan%z%C
      d = plan%z%d

      CALL GTStart(htot)

      IF (gdev_active) THEN
!
! Device path: FC-Gram continuation, 1D FFT along z and
! normalization on the device
!
         CALL GTStart(hcont)
         IF (C .gt. 0) CALL fftp3d_cont_z_dev(n,C,d,plan%z%dir,inout)
         CALL GTStop(hcont); conttime = conttime + GTGetTime(hcont)
         CALL GTStart(hfft)
         CALL fftp3d_exec_c2c_dev(plan%planrc,inout,GFFT_FORWARD)
         CALL fftp3d_norm_z_dev(n,inout)
         CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      ELSE
!
! Host path (sphp): FC-Gram continuation
!
      CALL GTStart(hcont)

      IF (plan%z%C .ge. 0) THEN
      DO ii=1,C
         inout(:,n-C+ii) = plan%z%dir(ii,1)*inout(:,n-C-d+1) +&
                           plan%z%dir(C-ii+1,1)*inout(:,d)
      ENDDO
      DO jj=2,d
         DO ii=1,C
            inout(:,n-C+ii) = inout(:,n-C+ii) +&
                              plan%z%dir(ii,jj)*inout(:,n-C-d+jj) +&
                              plan%z%dir(C-ii+1,jj)*inout(:,d-jj+1)
         ENDDO
      ENDDO
      ENDIF

      CALL GTStop(hcont); conttime = conttime + GTGetTime(hcont)

!
! 1D FFT in each node using the FFTW library
!
      CALL GTStart(hfft)
      CALL GFFTW_EXECUTE_DFT(plan%planrc%plancz,inout,inout)
! Normalization, so that the coefficients are those of the field
      inout = inout/REAL(n,KIND=GP)
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      ENDIF

      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp1d_real_to_complex_z

!*****************************************************************
      SUBROUTINE fftp3d_complex_to_real(plan,in,out,comm)
!-----------------------------------------------------------------
! Computes the 3D inverse transform in parallel: the inverse FFT in
! z, the transposition between tasks and the inverse 2D spherical
! harmonic transform with SHTns. If the input is f(lm,kz) the
! output is f(theta,phi,z). The input data is destroyed during the
! computation (the FFT in z is done in place).
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     in   : complex input array f(lm,kz) [INOUT]
!     out  : real output array f(theta,phi,z) [OUT]
!     comm : the MPI communicator (handle) [IN]
!
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,plan%z%n) :: in 
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: out

      TYPE(MPI_Comm)   , INTENT(IN)            :: comm

      CALL fftp1d_complex_to_real_z(plan,in,comm)
      CALL fftp2d_complex_to_real_xy(plan,in,out,comm)
     
      RETURN
      END SUBROUTINE fftp3d_complex_to_real

!*****************************************************************
      SUBROUTINE fftp2d_complex_to_real_xy(plan,in,out,comm)
!-----------------------------------------------------------------
! Computes the inverse 2D spherical harmonic transform with SHTns,
! plane by plane, after the transposition between tasks. If the
! input has shape (nlm,nz), con los lm en el orden de SHTns, the
! output's shape is (nlat,nphi,nz). The input data is not
! destroyed: los bloques se mandan directamente desde in, y
! SH_to_spat no cambia los coeficientes.
! On the device path (gdev_active) the exchange and the transform
! run on the device.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     in   : complex input array f(lm,z) [IN]
!     out  : real output array f(theta,phi,z) [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fftplans
      USE fcgram
      USE gtimer
      USE gdevice
      USE fftp3d_gpu
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(lmsta:lmend,plan%z%n) :: in 
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: out

      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: k

      CALL GTStart(htot)
      IF (gdev_active) THEN
! gcarr1 and the exchange tables have the layout of the first plan
! created (see fftp3d_create_plan)
      IF (SIZE(gcarr1,1).ne.plan%plancr%nlm) &
         STOP 'fftp3d: the device path supports only plans with the layout of the first plan created'
!
! Device path: exchange between tasks, with the blocks unpacked on
! the device
!
         CALL GTStart(hcom)
         CALL fftp3d_exchange_bwd(plan%plancr,in,gcarr1,comm)
         CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)
!
! 2D transform of the (theta,phi) planes on the device
!
         CALL GTStart(hfft)
         CALL fftp3d_exec_c2r_dev(plan%plancr,gcarr1,out)
         CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)

      ELSE
!
! Host path (sphp). La transposición local de SPECTER no hace
! falta: in ya está en el orden (lm,z) en que se mandan los bloques
!
! Transposes the result between nodes using strip mining when
! nstrip>1 (rreddy@psc.edu). Cada tarea manda a las otras, desde in,
! sus lm de los planos que les tocan, y recibe en carr1 los lm de las
! otras para sus planos
!
      CALL fftp3d_redist_lm_to_k(plan,in,plan%plancr%carr1,comm)
!
! 2D transform in each plane using the SHTns library, from all the
! lm of the plane, received in carr1
!
      CALL GTStart(hfft)
      DO k = ksta,kend
         CALL SH_to_spat(plan%plancr%planrxy,plan%plancr%carr1(:,k),out(:,:,k))
      END DO
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      ENDIF
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_complex_to_real_xy

!*****************************************************************
      SUBROUTINE fftp1d_complex_to_real_z(plan,inout,comm)
!-----------------------------------------------------------------
! Computes the 1D backward FFT along the second index of the input
! array, (lm,kz). The operation is performed inplace. No normaliza:
! la normalización la hace la directa. On the device path
! (gdev_active) the FFT runs on the device.
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     in   : f(lm,kz) at input, f(lm,z) at output [INOUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE fprecision
      USE mpivars
      USE commtypes
      USE fcgram
      USE gtimer
      USE gdevice
      USE fftp3d_gpu
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,plan%z%n) :: inout 
      TYPE(MPI_Comm)   , INTENT(IN)            :: comm

      INTEGER :: i,j,k

      CALL GTStart(htot)
      CALL GTStart(hfft)
      IF (gdev_active) THEN
!
! Device path: 1D FFT along z on the device
!
         CALL fftp3d_exec_c2c_dev(plan%plancr,inout,GFFT_INVERSE)
      ELSE
!
! Host path (sphp): 1D FFT in each node using the FFTW library
!
         CALL GFFTW_EXECUTE_DFT(plan%plancr%plancz,inout,inout)
      ENDIF
      CALL GTStop(hfft); ffttime = ffttime + GTGetTime(hfft)
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp1d_complex_to_real_z

!*****************************************************************
      SUBROUTINE block3d(imin,imax,jmin,jmax,kmin,ista,iend, &
                        jsta,jend,ksta,kend,ioldtype,inewtype)
!-----------------------------------------------------------------
!
! Soubroutine for defining derived data types in 3D.
!
! Parameters
!     imin : the minimum value in the first dimension [IN]
!     imax : the maximum value in the first dimension [IN]
!     jmin : the minimum value in the second dimension [IN]
!     jmax : the maximum value in the second dimension [IN]
!     kmin : the minimum value in the third dimension [IN]
!     ista : start value of the block in the first dimension [IN]
!     iend : end value of the block in the first dimension [IN]
!     jsta : start value of the block in the second dimension [IN]
!     jend : end value of the block in the second dimension [IN]
!     ksta : start value of the block in the third dimension [IN]
!     kend : end value of the block in the third dimension [IN]
!     ioldtype: data type of the elements in the block [IN]
!     inewtype: the derived data type for the block [OUT]
!-----------------------------------------------------------------

      USE commtypes
      USE fftplans
      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: ista,iend
      INTEGER, INTENT(IN)  :: jsta,jend
      INTEGER, INTENT(IN)  :: ksta,kend
      INTEGER, INTENT(IN)  :: imin,imax
      INTEGER, INTENT(IN)  :: jmin,jmax,kmin
      TYPE(MPI_Datatype), INTENT(IN)  :: ioldtype
      TYPE(MPI_Datatype), INTENT(OUT) :: inewtype

      INTEGER(MPI_ADDRESS_KIND) :: ilb,isize,idist
      TYPE(MPI_Datatype)        :: itemp,itemp2
      INTEGER :: ilen,jlen,klen
      INTEGER :: ierr
      INTEGER :: blocklengths(1)
      INTEGER(MPI_ADDRESS_KIND) :: displs(1)
      TYPE(MPI_Datatype) :: types(1)

      CALL MPI_TYPE_GET_EXTENT(ioldtype,ilb,isize,ierr)
      ilen = iend-ista+1
      jlen = jend-jsta+1
      klen = kend-ksta+1
      CALL MPI_TYPE_VECTOR(jlen,ilen,imax-imin+1,ioldtype,itemp,ierr)
      idist = (imax-imin+1)*(jmax-jmin+1)*isize
      CALL MPI_TYPE_CREATE_HVECTOR(klen,1,idist,itemp,itemp2,ierr)
      CALL MPI_TYPE_FREE(itemp,ierr)
      idist = ((imax-imin+1)*(jmax-jmin+1)*(ksta-kmin) &
              +(imax-imin+1)*(jsta-jmin)+(ista-imin))*isize
      blocklengths(1) = 1
      displs(1) = idist
      types(1) = itemp2
      CALL MPI_TYPE_CREATE_STRUCT(1,blocklengths,displs,types,inewtype,ierr)
      CALL MPI_TYPE_FREE(itemp2,ierr)
      CALL MPI_TYPE_COMMIT(inewtype,ierr)

      RETURN
      END SUBROUTINE block3d

!*****************************************************************
      SUBROUTINE range(n1,n2,nprocs,irank,ista,iend)
!-----------------------------------------------------------------
! Soubroutine for computing the local coordinate range 
! when splitting the original array into the nodes
!
! Parameters
!     n1     : the minimum value in the splitted dimension [IN]
!     n2     : the maximum value in the splitted dimension [IN]
!     nprocs : the number of processors [IN]
!     irank  : the rank of the processor [IN]
!     ista   : start value for the local coordinate [OUT]
!     iend   : end value for the local coordinate [OUT]
!-----------------------------------------------------------------

      IMPLICIT NONE

      INTEGER, INTENT(IN)  :: n1,n2
      INTEGER, INTENT(IN)  :: nprocs,irank
      INTEGER, INTENT(OUT) :: ista,iend

      INTEGER :: myrank
      INTEGER :: iwork1,iwork2

      iwork1 = (n2-n1+1)/nprocs
      iwork2 = MOD(n2-n1+1,nprocs)
      ista = irank*iwork1+n1+MIN(irank,iwork2)
      iend = ista+iwork1-1
      IF (iwork2.gt.irank) iend = iend+1

      RETURN
      END SUBROUTINE range

!*****************************************************************
      SUBROUTINE fftp3d_redist_k_to_lm(plan,in,out,comm)
!-----------------------------------------------------------------
! Transposes the coefficients between nodes, from the planes of
! each node to the lm of each node: if the input has shape
! (nlm,ksta:kend) the output's shape is (lmsta:lmend,nz). Es el
! cambio de reparto de la transformada directa; con una sola tarea
! es una copia. Se manda con itype1 y se recibe con itype2 (see
! fftp3d_create_block).
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     in   : complex input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fcgram
      USE gtimer
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(plan%planrc%nlm,ksta:kend) :: in
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,plan%z%n) :: out

      TYPE(MPI_Request), DIMENSION(0:nprocs-1) :: ireq1,ireq2
      TYPE(MPI_Status)                         :: istatus
      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc

!
! Transposes the data between nodes using strip mining when
! nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs
            CALL MPI_IRECV(out,1,plan%planrc%itype2(igetFrom),igetFrom, &
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(in,1,plan%planrc%itype1(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)

      RETURN
      END SUBROUTINE fftp3d_redist_k_to_lm

!*****************************************************************
      SUBROUTINE fftp3d_redist_lm_to_k(plan,in,out,comm)
!-----------------------------------------------------------------
! Transposes the coefficients between nodes, from the lm of each
! node to the planes of each node: if the input has shape
! (lmsta:lmend,nz) the output's shape is (nlm,ksta:kend). Es el
! cambio de reparto de la transformada inversa, el camino de
! fftp3d_redist_k_to_lm al revés: se manda con itype2 y se recibe con
! itype1 (see fftp3d_create_block).
!
! Parameters
!     plan : the FCPLAN plan (see fcgram_mod.f90) [IN]
!     in   : complex input array [IN]
!     out  : complex output array [OUT]
!     comm : the MPI communicator (handle) [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE fprecision
      USE mpivars
      USE fcgram
      USE gtimer
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(lmsta:lmend,plan%z%n) :: in
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(plan%plancr%nlm,ksta:kend) :: out

      TYPE(MPI_Request), DIMENSION(0:nprocs-1) :: ireq1,ireq2
      TYPE(MPI_Status)                         :: istatus
      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: irank
      INTEGER :: isendTo,igetFrom
      INTEGER :: istrip,iproc

!
! Transposes the data between nodes using strip mining when
! nstrip>1 (rreddy@psc.edu)
!
      CALL GTStart(hcom)
      do iproc = 0, nprocs-1, nstrip
         do istrip=0, nstrip-1
            irank = iproc + istrip

            isendTo = myrank + irank
            if ( isendTo .ge. nprocs ) isendTo = isendTo - nprocs

            igetFrom = myrank - irank
            if ( igetFrom .lt. 0 ) igetFrom = igetFrom + nprocs
            CALL MPI_IRECV(out,1,plan%plancr%itype1(igetFrom),igetFrom, &
                          1,comm,ireq2(irank),ierr)
            CALL MPI_ISEND(in,1,plan%plancr%itype2(isendTo),isendTo, &
                          1,comm,ireq1(irank),ierr)
         enddo

         do istrip=0, nstrip-1
            irank = iproc + istrip
            CALL MPI_WAIT(ireq1(irank),istatus,ierr)
            CALL MPI_WAIT(ireq2(irank),istatus,ierr)
         enddo
      enddo
      CALL GTStop(hcom); comtime = comtime + GTGetTime(hcom)

      RETURN
      END SUBROUTINE fftp3d_redist_lm_to_k
