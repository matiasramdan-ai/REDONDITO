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
! transform in (theta,phi) and the FFTW plan for the FFT in z.
!
! Parameters
!     plan   : contains the parallel 3D plan [OUT]
!     n      : the size of the dimensions of the input array,
!              (nlat,nphi,nz) [IN]
!     fftdir : the direction of the Fourier transform [IN]
!              FFTW_FORWARD or FFTW_REAL_TO_COMPLEX (-1)
!              FFTW_BACKWARD or FFTW_COMPLEX_TO_REAL (+1)
!     flags  : flags for the FFTW [IN]
!              FFTW_ESTIMATE (sub-optimal but faster)
!              FFTW_MEASURE (optimal but slower to create plans)
!              FFTW_PATIENT AND FFTW_EXHAUSTIVE are also available
!              for extra performance, but may take a long time to
!              create plans (specially when using OpenMP)
!     lmax   : grado máximo de los armónicos esféricos [IN]
!-----------------------------------------------------------------

      USE commtypes
      USE mpivars
      USE fftplans
!$    USE threads
      USE gtimer
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: n(3)
      INTEGER, INTENT(IN) :: fftdir
      INTEGER, INTENT(IN) :: flags
      INTEGER, INTENT(IN) :: lmax
      TYPE(FFTPLAN), INTENT(OUT) :: plan

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
                         fftdir,flags)
      plan%nlat = n(1)
      plan%nphi = n(2)
      plan%nz   = n(3)
      plan%nlm  = (lmax+1)*(lmax+2)/2
      ALLOCATE( plan%itype1(0:nprocs-1) )
      ALLOCATE( plan%itype2(0:nprocs-1) )
      CALL fftp3d_create_block(n,(lmax+1)*(lmax+2)/2,nprocs,myrank, &
              plan%itype1,plan%itype2)

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
! Destroys the SHTns and FFTW plans in each node.
!
! Parameters
!     plan : the parallel 3D plan [INOUT]
!-----------------------------------------------------------------

      USE fftplans
      USE gtimer
      IMPLICIT NONE

      TYPE(FFTPLAN), INTENT(INOUT) :: plan

      CALL shtns_destroy(plan%planrxy)
      CALL GFFTW_DESTROY_PLAN(plan%plancz)

      DEALLOCATE( plan%ccarr  )
      DEALLOCATE( plan%carr1  )
      DEALLOCATE( plan%carr2  )
      DEALLOCATE( plan%carr3  )
      DEALLOCATE( plan%itype1 )
      DEALLOCATE( plan%itype2 )

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
      USE fcgram
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(lmsta:lmend,plan%z%n)    :: out
      REAL(KIND=GP), INTENT(INOUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: IN

      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: k

      CALL GTStart(htot)
!
! 2D transform in each plane using the SHTns library. SHTns does
! the FFT in phi in place, so the input data is destroyed
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
! vuelta es la identidad, sin normalizar después.
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
      USE fcgram
      USE gtimer
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(lmsta:lmend,plan%z%n) :: in 
      REAL(KIND=GP), INTENT(OUT), DIMENSION(plan%x%n,plan%y%n,ksta:kend) :: out

      TYPE(MPI_Comm)   , INTENT(IN)            :: comm
      INTEGER :: k

      CALL GTStart(htot)
!
! La transposición local de SPECTER no hace falta: in ya está en el
! orden (lm,z) en que se mandan los bloques
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
      
      CALL GTStop(htot); tottime = tottime + GTGetTime(htot)

      RETURN
      END SUBROUTINE fftp2d_complex_to_real_xy

!*****************************************************************
      SUBROUTINE fftp1d_complex_to_real_z(plan,inout,comm)
!-----------------------------------------------------------------
! Computes the 1D backward FFT along the second index of the input
! array, (lm,kz). The operation is performed inplace. No normaliza:
! la normalización la hace la directa.
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
!$    USE threads
      IMPLICIT NONE

      TYPE(FCPLAN), INTENT(IN) :: plan
      COMPLEX(KIND=GP), INTENT(INOUT), DIMENSION(lmsta:lmend,plan%z%n) :: inout 
      TYPE(MPI_Comm)   , INTENT(IN)            :: comm

      INTEGER :: i,j,k

!
! 1D FFT in each node using the FFTW library
      CALL GTStart(htot)
      CALL GTStart(hfft)
      CALL GFFTW_EXECUTE_DFT(plan%plancr%plancz,inout,inout)
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
