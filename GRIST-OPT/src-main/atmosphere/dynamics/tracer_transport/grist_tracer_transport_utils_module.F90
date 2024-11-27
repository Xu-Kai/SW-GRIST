
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Version 1.0
! Description: Some utilities for tracer transport
! Revision history:
!                  1. fixer should not be used currently
!----------------------------------------------------------------------------

  module grist_tracer_transport_utils_module

    use grist_domain_types, only: global_domain
    use grist_data_types,   only: scalar_2d_field, scalar_3d_field
    use grist_constants,    only: i4, r8, zero, one, half
    use grist_nml_module,   only: nlev, ntracer, write_verbose, physpkg
    use grist_hpe_constants,only: deta_full, deta_face
    use grist_mpi
! data
    use grist_tracer_transport_vars_module, only: scalar_tracer_mxrt_min

    implicit none

    save

    PUBLIC :: tracer_transport_check_mxrt, &
              tracer_transport_fixer_mxrt, & ! low-efficient, needs optimization
              tracer_transport_qneg3_mxrt

  CONTAINS

   subroutine tracer_transport_check_mxrt(mesh, scalar_tracer_mxrt_at_pc_full_level_n,called_env)
! io
     use omp_lib
     type(global_domain),   intent(in) :: mesh
     type(scalar_3d_field), intent(in) :: scalar_tracer_mxrt_at_pc_full_level_n
     character(len=*),      intent(in) :: called_env
! local
     integer(i4)     :: itracer, iv, ilev
     integer(i4)     :: nsum, nsum_global, ierr
     integer(i4)     :: ii

#if defined(CHK_TRACER_TRANS_CHECK_MXRT) || defined(CHK_ALL)
#include "data_check.inc"
      integer(i4)   :: ref_nsum
#endif     

    nsum = 0

call t_startf("tracer_transport_check_mxrt_target")
!$omp target map(tofrom:nsum)
!$omp parallel  private(ii,iv,ilev,itracer) 
!$omp  do  reduction(+:nsum)
        do ii = 1, mesh%nv_compute*nlev*ntracer,1
         iv=ceiling(ii/real(nlev*ntracer,r8))
         ilev=ceiling((ii-(iv-1)*nlev*ntracer)/real(ntracer,r8))
         itracer=ii-(iv-1)*nlev*ntracer-(ilev-1)*ntracer
 !        do iv = 1, mesh%nv_compute
 !           do itracer = 1, ntracer
 !              do ilev = 1, nlev
                  if(scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev,iv).lt.scalar_tracer_mxrt_min(itracer))then
                     nsum = nsum + 1
                  end if
 !              end do
 !           end do
 !        end do
      end do
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("tracer_transport_check_mxrt_target")
#if defined(CHK_TRACER_TRANS_CHECK_MXRT) || defined(CHK_ALL)
      ref_nsum = 0
call  t_startf("tracer_transport_check_mxrt")
!$omp parallel  private(ii,iv,ilev,itracer) 
!$omp  do  reduction(+:ref_nsum)
     do ii = 1, mesh%nv_compute*nlev*ntracer,1
        iv=ceiling(ii/real(nlev*ntracer,r8))
        ilev=ceiling((ii-(iv-1)*nlev*ntracer)/real(ntracer,r8))
        itracer=ii-(iv-1)*nlev*ntracer-(ilev-1)*ntracer
!        do iv = 1, mesh%nv_compute
!           do itracer = 1, ntracer
!              do ilev = 1, nlev
                 if(scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev,iv).lt.scalar_tracer_mxrt_min(itracer))then
                    ref_nsum = ref_nsum + 1
                 end if
!              end do
!           end do
!        end do
     end do
!$omp end do nowait
!$omp end parallel
call t_stopf("tracer_transport_check_mxrt")
     call data_check(nsum, ref_nsum)
#endif

#ifndef SEQ_GRIST
        call reduce(nsum, nsum_global, 'sum')
#else
        nsum_global = nsum
#endif
        if(mpi_rank().eq.0.and.write_verbose) print*,"nsum_global=",nsum_global,"in violation of min",trim(called_env)

      return
    end subroutine tracer_transport_check_mxrt

!----------------------------------------------------------------------------------
! This filling procedure is similar to that used by CAM-FV; Although horizontal
! transport is monotonic and refuses undershooting to machine round off, certain
! small negative values may still be generated in some ill conditions
! (diabatic/hvsplit/boundary)
!----------------------------------------------------------------------------------

    subroutine tracer_transport_fixer_mxrt(ncell,scalar_tracer_mxrt_at_pc_full_level_n,&
                                                 scalar_delhp_at_pc_full_level_n)
! io
     integer(i4),           intent(in)    :: ncell
     type(scalar_3d_field), intent(inout) :: scalar_tracer_mxrt_at_pc_full_level_n
     type(scalar_2d_field), intent(in)    :: scalar_delhp_at_pc_full_level_n
! local
     integer(i4)     :: itracer, iv, ilev
     integer(i4)     :: nsum, nsum_global, ierr
     real(r8)        :: qkold, qfill_up, qfill_lo, qfill

        do iv = 1, ncell
           do itracer = 1, ntracer
              do ilev = 1, nlev
                 if(scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev,iv).lt.scalar_tracer_mxrt_min(itracer))then
             ! as long as tracer_mxrt violates minimum constraint, it will be reset to min value,
             ! regardless whether we borrow from upper or lower levels; when not
             ! borrowing, this will violate mass conservation, but such computational
             ! sources are quite small.
                    qkold = scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev,iv)
                    scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev,iv) = scalar_tracer_mxrt_min(itracer)
             ! borrow from upper or lower levels
                    qfill_up = zero
                    qfill_lo = zero
                    if(ilev.ne.1)    qfill_up = (scalar_tracer_mxrt_min(itracer)-qkold)*&
                                scalar_delhp_at_pc_full_level_n%f(ilev,iv)/scalar_delhp_at_pc_full_level_n%f(ilev-1,iv)
                    if(ilev.ne.nlev) qfill_lo = (scalar_tracer_mxrt_min(itracer)-qkold)*&
                                scalar_delhp_at_pc_full_level_n%f(ilev,iv)/scalar_delhp_at_pc_full_level_n%f(ilev+1,iv)
             ! check how to borrow:
             ! 1) from upper level:
#ifndef NEWFIXER
                    if(ilev.ne.1.and.(scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev-1,iv)-qfill_up).gt.scalar_tracer_mxrt_min(itracer))then
#else
                    if(ilev.ne.1)then
                       if((scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev-1,iv)-qfill_up).gt.scalar_tracer_mxrt_min(itracer))then

#endif
                       scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev-1,iv) = scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev-1,iv)-qfill_up
                       cycle
#ifdef NEWFIXER
                       end if
#endif
             ! 2) from lower level:
#ifndef NEWFIXER
                    else if(ilev.ne.nlev.and.(scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev+1,iv)-qfill_lo).gt.scalar_tracer_mxrt_min(itracer))then
#else
                    else if(ilev.ne.nlev)then
                       if((scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev+1,iv)-qfill_lo).gt.scalar_tracer_mxrt_min(itracer))then
#endif
                       scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev+1,iv) = scalar_tracer_mxrt_at_pc_full_level_n%f(itracer,ilev+1,iv)-qfill_lo
                       cycle
#ifdef NEWFIXER
                       endif
#endif
             ! 3) not borrow:
                    else
                       cycle
                    end if
                 end if
              end do
           end do
        end do
       
      return
    end subroutine tracer_transport_fixer_mxrt
!
! this code is adapted from CAM's qneg3 code
!
 subroutine tracer_transport_qneg3_mxrt(ncell,nlev,ntracer,qqq,called_env)
! io
   use ieee_arithmetic
   integer(i4), intent(in)    :: ncell       ! number of atmospheric columns
   integer(i4), intent(in)    :: nlev        ! number of vertical levels in column
   integer(i4), intent(in)    :: ntracer     ! number of vertical levels in column
   real(r8),    intent(inout) :: qqq(:,:,:)  ! moisture/tracer field
   character*(*), intent(in)  :: called_env  ! name of calling routine
! local
   integer indx(nlev,ncell) ! array of indices of points < qmin
   integer nval(nlev)       ! number of points < qmin for 1 level
   integer nvals            ! number of values found < qmin
   integer nn
   integer iwtmp
   integer iv,ii,k          ! longitude, level indices
   integer m                ! constituent index
   integer iw,kw            ! i,k indices of worst violator
   logical found            ! true => at least 1 minimum violator found
   real(r8) worst           ! biggest violator
   logical amipc
   real(r8) tracer_lower_bound(ntracer), tracer_upper_bound(ntracer), infty
   ! infty = real(z'7ff0000000000000',r8)
!#define CHK_TRACER_QNEG3
#if defined(CHK_TRACER_QNEG3) || defined(CHK_ALL)
#include "data_check.inc"
   real(r8), allocatable :: ref_qqq(:,:,:) 
   allocate(ref_qqq, source=qqq)
#endif
   amipc = trim(physpkg).eq.'AMIPC_PHYSICS'
   infty = ieee_value(infty, ieee_positive_inf)
   do m = 1, ntracer
      tracer_lower_bound(m) = scalar_tracer_mxrt_min(m)
      tracer_upper_bound(m) = infty
   end do
   if (amipc) then
      tracer_lower_bound(4:5) = 1.e-12_r8
      tracer_upper_bound(4:5) = 1.e10_r8
   end if
   !$omp target parallel do private(iv, k, m)
   do iv = 1, ncell
      do k = 1,nlev
         do m = 1, ntracer
            qqq(m,k,iv) = max(tracer_lower_bound(m), qqq(m,k,iv))
            qqq(m,k,iv) = min(tracer_upper_bound(m), qqq(m,k,iv))
         end do
      end do
   end do
   !$omp end target parallel do
#if defined(CHK_TRACER_QNEG3) || defined(CHK_ALL)
   do m = 1, ntracer
      nvals = 0
      found = .false.
      worst = 1.e35_r8
      iw    = -1
!
! Test all field values for being less than minimum value. Set q = qmin
! for all such points. Trace offenders and identify worst one.
!
      if(trim(physpkg).eq.'AMIPC_PHYSICS'  .and. (m == 4 .or. m == 5))then
          ref_qqq(m,:,:) = max(1.e-12_r8,ref_qqq(m,:,:))
          ref_qqq(m,:,:) = min(1.e10_r8, ref_qqq(m,:,:))
      else 
      do k = 1,nlev
         nval(k) = 0
         nn      = 0
         do iv = 1,ncell
            if(ref_qqq(m,k,iv) < scalar_tracer_mxrt_min(m))then
               nn = nn + 1
               indx(k,nn) = iv
            end if
         end do
         nval(k) = nn ! how many for one level
      end do

      do k = 1,nlev
         if (nval(k) > 0) then
            found = .true.
            nvals = nvals + nval(k) ! total count
            iwtmp = -1
! count
            do ii = 1, nval(k)
               iv = indx(k,ii)
               if (ref_qqq(m,k,iv) < worst) then
                  worst = ref_qqq(m,k,iv)
                  iwtmp = ii
               end if
            end do

            if (iwtmp /= -1 ) kw = k
            if (iwtmp /= -1 ) iw = indx(k,iwtmp)
! modify in this way is good
!!!            do ii = 1, nval(k)
!!!               iv = indx(k,ii)
!!!               ref_qqq(m,k,iv) = scalar_tracer_mxrt_min(m)
!!!            end do
            where(ref_qqq(m,:,:) < scalar_tracer_mxrt_min(m)) ref_qqq(m,:,:) = scalar_tracer_mxrt_min(m)

         end if
      end do
  !    if (found .and. abs(worst)>1.e-12_r8) write(6,9000) called_env, m,nvals,iwtmp,scalar_tracer_mxrt_min(m),worst,iw,kw
      end if
   end do
   call data_check(qqq, ref_qqq)
#endif
   return
9000 format(' QNEG3 from ',a,':m=',i3,' ncell=',i5, &
            ' Min. mixing ratio violated at ',i4,' points.  Reset to ', &
            1p,e8.1,' Worst =',e8.1,' at i,k=',i4,i3)

  end subroutine tracer_transport_qneg3_mxrt

  end module grist_tracer_transport_utils_module
