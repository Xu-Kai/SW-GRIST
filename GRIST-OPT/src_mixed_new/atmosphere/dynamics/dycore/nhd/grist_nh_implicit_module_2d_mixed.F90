
!----------------------------------------------------------------------------
! Created on 2018
! Author: Yi Zhang
! Version 1.0
! Description: This module is used for implicitly solving the W equation. 
!              It deals with vertical column data structure
!              There is evidence that the acoustic solver is prone to
!              suffer from an intel optimzation bug, reasons unknown
! Revision history:
!----------------------------------------------------------------------------

  module grist_nh_implicit_module_2d_mixed

   use grist_constants,                only: i4, r8, gravity, cp, cv, rearth, fillvalue, zero, one, r4 => ns
   use grist_nml_module,               only: imbeta, imkesi
   use grist_tridiagonal_solver_mixed, only: solve_tridiagonal
   use grist_domain_types,             only: global_domain

   implicit none

   private
   public   :: grist_nh_solve_www_equation

   CONTAINS

   subroutine grist_nh_solve_www_equation(mesh, nlev, nlevp, &
                                          scalar_pressure_full_level_n,        &
                                          scalar_delhp_full_level_n,           &
                                          scalar_potential_temp_full_level_n,  &
                                          scalar_potential_temp_full_level_np1,&
                                          scalar_delhp_full_level_np1,         &
                                          scalar_phi_face_level_n,             &
                                          scalar_www_face_level_n,             &
                                          scalar_www_face_level_rk,            &
                                          scalar_et_www_face_level_rk,         &
                                          scalar_et_phi_face_level_rk,         &
                                          scalar_delhp_face_level_np1,         &
                                          scalar_delp_face_level_rk,           &
                                          scalar_delhp_face_level_rk,          &
                                          scalar_www_face_level,               &
                                          scalar_mif_face_level,               &
                                          scalar_mif_full_level,               &
                                          dtime,                               &
                                          www_top, www_bot)
      use omp_lib
! io
      type(global_domain),   intent(in)    :: mesh
      integer(i4)        ,   intent(in)    :: nlev, nlevp
      real(r8), allocatable, intent(in)    :: scalar_pressure_full_level_n(:,:)
      real(r8), allocatable, intent(in)    :: scalar_delhp_full_level_n(:,:)
      real(r8), allocatable, intent(in)    :: scalar_potential_temp_full_level_n(:,:)
      real(r8), allocatable, intent(in)    :: scalar_potential_temp_full_level_np1(:,:)
      real(r8), allocatable, intent(in)    :: scalar_delhp_full_level_np1(:,:)
      real(r8), allocatable, intent(in)    :: scalar_phi_face_level_n(:,:)
      real(r8), allocatable, intent(in)    :: scalar_www_face_level_n(:,:)
      real(r4), allocatable, intent(in)    :: scalar_www_face_level_rk(:,:)
      real(r4), allocatable, intent(in)    :: scalar_et_www_face_level_rk(:,:)
      real(r4), allocatable, intent(in)    :: scalar_et_phi_face_level_rk(:,:)
      real(r8), allocatable, intent(in)    :: scalar_delhp_face_level_np1(:,:)
      real(r8), allocatable, intent(in)    :: scalar_delp_face_level_rk(:,:)
      real(r8), allocatable, intent(in)    :: scalar_delhp_face_level_rk(:,:)
      real(r8), allocatable, intent(inout) :: scalar_www_face_level(:,:)
      real(r8), allocatable, intent(in)    :: scalar_mif_face_level(:,:)
      real(r8), allocatable, intent(in)    :: scalar_mif_full_level(:,:)
      real(r8)             , intent(in)    :: dtime
      real(r8)             , intent(in)    :: www_top(mesh%nv_full), www_bot(mesh%nv_full)
! local
      real(r8)                :: scalar_aaa_full_level(nlev ,mesh%nv_full)
      real(r8)                :: scalar_bbb_face_level(nlevp,mesh%nv_full)
      real(r8)                :: scalar_rrr_face_level(nlevp,mesh%nv_full)
      real(r8)                :: scalar_aaa_face_level(nlev-1)
      real(r8)                :: scalar_bbb_fcut_level(nlev-1)
      real(r8)                :: scalar_ccc_face_level(nlev-1)
      real(r8)                :: scalar_rrr_fcut_level(nlev-1)
      integer(i4)             :: ilev, nsize,iv, ii
      real(r8)                :: gama


       !allocate(scalar_aaa_full_level(nlev  ,mesh%nv_full))
       !allocate(scalar_bbb_face_level(nlev+1,mesh%nv_full))
       !allocate(scalar_rrr_face_level(nlev+1,mesh%nv_full))
call t_startf("nh_solve_www_equation_copy")
!$omp target
!$omp parallel workshare
       scalar_aaa_full_level = zero
       scalar_bbb_face_level = zero
       scalar_rrr_face_level = zero
!$omp end parallel workshare 
!$omp end target
       call t_stopf("nh_solve_www_equation_copy")

       !call grist_nh_compute_aaa(mesh,nlev,&
       !                          scalar_pressure_full_level_n ,&
       !                          scalar_delhp_full_level_n    ,&
       !                          scalar_delhp_full_level_np1  ,&
       !                          scalar_phi_face_level_n      ,&
       !                          scalar_mif_full_level        ,&
       !                          scalar_aaa_full_level, dtime) ! out

! inline opt here
       gama       = cp/cv
call t_startf("nh_www_equation_target_1")
!$omp target
!$omp parallel  private(ii,iv,ilev)
!$omp do
      do ii = 1, mesh%nv_halo(1)*nlev,1

         iv=ceiling(ii/real(nlev,r8))
         ilev=ii-(iv-1)*nlev
!      do iv      = 1, mesh%nv_halo(1)
!      do ilev    = 1, nlev

#if !defined NEWLGE && !defined LGE1
         scalar_aaa_full_level(ilev,iv) = scalar_mif_full_level(ilev,iv)*gama*imbeta*imkesi*((gravity*dtime)**2)*&
                                         (scalar_pressure_full_level_n(ilev,iv))/(scalar_phi_face_level_n(ilev+1,iv)-scalar_phi_face_level_n(ilev,iv))
#endif
! #ifdef LGE1
!          scalar_aaa_full_level(ilev,iv) = scalar_mif_full_level(ilev,iv)*gama*imbeta*imkesi*((gravity*dtime)**2)*&
!                                        (scalar_pressure_full_level_n(ilev,iv))/(scalar_phi_face_level_n(ilev+1,iv)-scalar_phi_face_level_n(ilev,iv))*&
!                                         scalar_delhp_full_level_n(ilev,iv)/scalar_delhp_full_level_np1(ilev,iv)
! #endif
!      end do
!      end do
      end do
!$omp end do nowait
!$omp end parallel 
!$omp end target 
      call t_stopf("nh_www_equation_target_1")
! inline opt here

       !call grist_nh_compute_bbb(mesh,nlev,&
       !                          scalar_delhp_face_level_np1 , &
       !                          scalar_aaa_full_level       , &
       !                          scalar_bbb_face_level)        ! out

! inline opt here
call t_startf("nh_www_equation_target_2")
!$omp target
!$omp parallel  private(ii,iv,ilev)
!$omp do
      do ii = 1, mesh%nv_halo(1)*(nlev-1),1
         iv=ceiling(ii/real((nlev-1),r8))
         ilev=1+ii-(iv-1)*(nlev-1)
!       do iv   = 1, mesh%nv_halo(1)
!       do ilev = 2, nlev
          scalar_bbb_face_level(ilev,iv) = scalar_delhp_face_level_np1(ilev,iv)-scalar_aaa_full_level(ilev,iv)-scalar_aaa_full_level(ilev-1,iv)
!       end do
!       end do
      end do
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("nh_www_equation_target_2")

call t_startf("nh_solve_www_equation_copy_2")
!$omp target
!$omp parallel workshare
      scalar_bbb_face_level(1,:)      = fillvalue
      scalar_bbb_face_level(nlev+1,:) = fillvalue
!$omp end parallel workshare
!$omp end target
      call t_stopf("nh_solve_www_equation_copy_2")
! inline opt here

       call grist_nh_compute_rrr(mesh,&
                                 scalar_www_face_level_rk,            &
                                 scalar_et_www_face_level_rk,         &
                                 scalar_et_phi_face_level_rk,         &
                                 scalar_www_face_level_n,             &
                                 scalar_phi_face_level_n,             &
                                 scalar_delhp_face_level_np1,         &
                                 scalar_pressure_full_level_n,        &
                                 scalar_delhp_full_level_n,           &
                                 scalar_potential_temp_full_level_n,  &
                                 scalar_potential_temp_full_level_np1,&
                                 scalar_delhp_full_level_np1,         &
                                 scalar_delp_face_level_rk,           &
                                 scalar_delhp_face_level_rk,          &
                                 scalar_mif_face_level,               &
                                 scalar_rrr_face_level,               & ! out
                                 dtime, mesh%nv_full, nlev, nlevp)

!       call grist_nh_compute_www(mesh,nlev,&
!                                 scalar_aaa_full_level, &   ! in
!                                 scalar_bbb_face_level, &   ! in
!                                 scalar_rrr_face_level, &   ! in
!                                 www_top, www_bot  , &   ! in
!                                 scalar_www_face_level)     ! out

! inline opt

       nsize = nlev-1
call t_startf("nh_www_equation_target_3")
!$omp target
!$omp parallel  private(iv,scalar_aaa_face_level,scalar_ccc_face_level,scalar_bbb_fcut_level,scalar_rrr_fcut_level) 
!$omp do
       do iv = 1, mesh%nv_halo(1)

          scalar_www_face_level(1,iv)      = www_top(iv)
          scalar_www_face_level(nlev+1,iv) = www_bot(iv)

          scalar_aaa_face_level(1)         = zero;  scalar_aaa_face_level(2:nsize)   = scalar_aaa_full_level(2:nsize,iv)
          scalar_ccc_face_level(nsize)     = zero;  scalar_ccc_face_level(1:nsize-1) = scalar_aaa_full_level(2:nsize,iv)
          scalar_bbb_fcut_level(1:nsize)   = scalar_bbb_face_level(2:nsize+1,iv)
          scalar_rrr_fcut_level(1:nsize)   = scalar_rrr_face_level(2:nsize+1,iv)

          call solve_tridiagonal(scalar_aaa_face_level, &
                                 scalar_bbb_fcut_level, &
                                 scalar_ccc_face_level, &
                                 scalar_rrr_fcut_level, &
                                 scalar_www_face_level(2:nsize+1,iv), nsize)
       end do
!$omp end do nowait
!$omp end parallel
!$omp end target 
call t_stopf("nh_www_equation_target_3")
! inline opt

       !deallocate(scalar_aaa_full_level)
       !deallocate(scalar_bbb_face_level)
       !deallocate(scalar_rrr_face_level)

       return
   end subroutine grist_nh_solve_www_equation

!================================================
!         PRIVATE with 1d data structure
!================================================

   subroutine grist_nh_compute_aaa(mesh, nlev, &
                                   scalar_pressure_full_level_n,&
                                   scalar_delhp_full_level_n,   &
                                   scalar_delhp_full_level_np1, &
                                   scalar_phi_face_level_n,     &
                                   scalar_mif_full_level,       &
                                   scalar_aaa_full_level, dtime)
! io
     use omp_lib
     type(global_domain),   intent(in)    :: mesh
     integer(i4),           intent(in)    :: nlev
     real(r8), allocatable, intent(in)    :: scalar_pressure_full_level_n(:,:)  ! n means last direct step
     real(r8), allocatable, intent(in)    :: scalar_delhp_full_level_n(:,:)
     real(r8), allocatable, intent(in)    :: scalar_delhp_full_level_np1(:,:)
     real(r8), allocatable, intent(in)    :: scalar_phi_face_level_n(:,:)
     real(r8), allocatable, intent(in)    :: scalar_mif_full_level(:,:)
     real(r8), allocatable, intent(inout) :: scalar_aaa_full_level(:,:)
     real(r8),              intent(in)    :: dtime
! local
     real(r8)       :: gama
     integer(i4)    :: iv, ilev
     integer(i4)    :: ii

      gama       = cp/cv

!$omp parallel  private(ii,iv,ilev)  
!$omp do schedule(static,100) 
     do ii = 1, mesh%nv_halo(1)*nlev,1
        iv=ceiling(ii/real(nlev,r8))
        ilev=ii-(iv-1)*nlev
!      do iv      = 1, mesh%nv_halo(1)
!      do ilev    = 1, nlev
#if !defined NEWLGE && !defined LGE1
         scalar_aaa_full_level(ilev,iv) = scalar_mif_full_level(ilev,iv)*gama*imbeta*imkesi*((gravity*dtime)**2)*&
                                       (scalar_pressure_full_level_n(ilev,iv))/(scalar_phi_face_level_n(ilev+1,iv)-scalar_phi_face_level_n(ilev,iv))
#endif
#ifdef LGE1
         scalar_aaa_full_level(ilev,iv) = scalar_mif_full_level(ilev,iv)*gama*imbeta*imkesi*((gravity*dtime)**2)*&
                                       (scalar_pressure_full_level_n(ilev,iv))/(scalar_phi_face_level_n(ilev+1,iv)-scalar_phi_face_level_n(ilev,iv))*&
                                       scalar_delhp_full_level_n(ilev,iv)/scalar_delhp_full_level_np1(ilev,iv)
#endif
!      end do
!      end do
     end do
!$omp end do nowait
!$omp end parallel 

     return
   end subroutine grist_nh_compute_aaa

   subroutine grist_nh_compute_bbb(mesh, nlev, &
                                   scalar_delhp_face_level_np1,&
                                   scalar_aaa_full_level, &
                                   scalar_bbb_face_level)
! io 
     use omp_lib
     type(global_domain),   intent(in)    :: mesh
     integer(i4),           intent(in)    :: nlev
     real(r8), allocatable, intent(in)    :: scalar_delhp_face_level_np1(:,:)
     real(r8), allocatable, intent(in)    :: scalar_aaa_full_level(:,:)
     real(r8), allocatable, intent(inout) :: scalar_bbb_face_level(:,:)
! local
     integer(i4) :: ilev, iv
     integer(i4)    :: ii

!$omp parallel  private(ii,iv,ilev)
!$omp do schedule(static,100) 
     do ii = 1, mesh%nv_halo(1)*(nlev-1),1
        iv=ceiling(ii/real((nlev-1),r8))
        ilev=1+ii-(iv-1)*(nlev-1)
!      do iv   = 1, mesh%nv_halo(1)
!      do ilev = 2, nlev
         scalar_bbb_face_level(ilev,iv) = scalar_delhp_face_level_np1(ilev,iv)-scalar_aaa_full_level(ilev,iv)-scalar_aaa_full_level(ilev-1,iv)
!      end do
!      end do
     end do
!$omp end do nowait
!$omp end parallel 

      scalar_bbb_face_level(1,:)      = fillvalue 
      scalar_bbb_face_level(nlev+1,:) = fillvalue 

     return
   end subroutine grist_nh_compute_bbb

   subroutine grist_nh_compute_rrr(mesh,&
                                   scalar_www_face_level_rk,        &
                                   scalar_et_www_face_level_rk,     &
                                   scalar_et_phi_face_level_rk,     &
                                   scalar_www_face_level_n,         &
                                   scalar_phi_face_level_n,         &
                                   scalar_delhp_face_level_np1,     &
                                   scalar_pressure_full_level_n,    &
                                   scalar_delhp_full_level_n,       &
                                   scalar_potential_temp_full_level_n,   &
                                   scalar_potential_temp_full_level_np1, &
                                   scalar_delhp_full_level_np1,      &
                                   scalar_delp_face_level_rk,        &
                                   scalar_delhp_face_level_rk,       &
                                   scalar_mif_face_level,            &
                                   scalar_rrr_face_level,            &
                                   dtime, nv_full, nlev, nlevp)
! io
     use omp_lib
     type(global_domain),   intent(in)    :: mesh
     real(r4), allocatable, intent(in)    :: scalar_www_face_level_rk(:,:)
     real(r4), allocatable, intent(in)    :: scalar_et_www_face_level_rk(:,:)
     real(r4), allocatable, intent(in)    :: scalar_et_phi_face_level_rk(:,:)
     real(r8), allocatable, intent(in)    :: scalar_www_face_level_n(:,:)
     real(r8), allocatable, intent(in)    :: scalar_phi_face_level_n(:,:)
     real(r8), allocatable, intent(in)    :: scalar_delhp_face_level_np1(:,:)
     real(r8), allocatable, intent(in)    :: scalar_pressure_full_level_n(:,:)
     real(r8), allocatable, intent(in)    :: scalar_delhp_full_level_n(:,:) 
     real(r8), allocatable, intent(in)    :: scalar_potential_temp_full_level_n(:,:)
     real(r8), allocatable, intent(in)    :: scalar_potential_temp_full_level_np1(:,:)
     real(r8), allocatable, intent(in)    :: scalar_delhp_full_level_np1(:,:)
     real(r8), allocatable, intent(in)    :: scalar_delp_face_level_rk(:,:)
     real(r8), allocatable, intent(in)    :: scalar_delhp_face_level_rk(:,:)
     real(r8), allocatable, intent(in)    :: scalar_mif_face_level(:,:)
     real(r8),              intent(inout) :: scalar_rrr_face_level(nlevp,mesh%nv_full)
     real(r8), intent(in)                 :: dtime
     integer(i4),           intent(in)    :: nv_full, nlev, nlevp
! local
!      real(r8)                             :: part_a(nv_full), part_b(nv_full)
     real(r8)                             :: part_a_s, part_b_s
!      real(r8)                             :: delta_pt_a(nv_full), delta_pt_b(nv_full), delta_pt(nv_full)
     real(r8)                             :: delta_pt_a_s, delta_pt_b_s, delta_pt_s
!      real(r8)                             :: delta_big_a(nv_full), delta_big_b(nv_full),delta_big(nv_full)
     real(r8)                             :: delta_big_a_s, delta_big_b_s,delta_big_s
!      real(r8)                             :: phi1_k(nv_full), phi1_kp1(nv_full), phi1_km1(nv_full), delta_phi1(nv_full)
     real(r8)                             :: phi1_k_s, phi1_kp1_s, phi1_km1_s, delta_phi1_s
     real(r8)                             :: gama
     integer(i4)                          :: ilev, iv

      gama = cp/cv
!
! use exploop for default since 202105 (faster), but use imploop for regressing old results
!

#ifndef NHRRR_IMPLOOP
call t_startf("nh_compute_rrr_target")
!$omp target
!$omp parallel  private(iv,ilev,part_a_s,part_b_s,delta_pt_a_s,delta_pt_b_s,delta_pt_s,delta_big_a_s,delta_big_b_s,delta_big_s,phi1_k_s,phi1_kp1_s,phi1_km1_s,delta_phi1_s) 
!$omp do
      do iv = 1, mesh%nv_halo(1)
         do ilev = 2, nlev ! each face level

            part_a_s      = scalar_delhp_face_level_np1(ilev,iv)*&
! #ifdef NHDWPIER
!                               (scalar_www_face_level_n(ilev,iv)-scalar_www_face_level_rk(ilev,iv)-dtime*scalar_et_www_face_level_rk(ilev,iv)-&
! #elif defined (NHDWPIEN)
!                               (-dtime*scalar_et_www_face_level_rk(ilev,iv)-&
! #else
                              (scalar_www_face_level_n(ilev,iv)-dtime*scalar_et_www_face_level_rk(ilev,iv)-&
! #endif
                               gravity*dtime+&
                               scalar_mif_face_level(ilev,iv)*gravity*(1._r8-imkesi)*dtime*scalar_delp_face_level_rk(ilev,iv)/scalar_delhp_face_level_rk(ilev,iv))

#if !defined NEWLGE && !defined LGE1
            delta_pt_a_s  = gama*scalar_pressure_full_level_n(ilev,iv)*(scalar_potential_temp_full_level_np1(ilev,iv)/scalar_potential_temp_full_level_n(ilev,iv))
            delta_pt_b_s  = gama*scalar_pressure_full_level_n(ilev-1,iv)*(scalar_potential_temp_full_level_np1(ilev-1,iv)/scalar_potential_temp_full_level_n(ilev-1,iv))
#endif
! #ifdef LGE1
!             delta_pt_a(iv)  = gama*scalar_pressure_full_level_n(ilev,iv)*(scalar_potential_temp_full_level_np1(ilev,iv)/scalar_potential_temp_full_level_n(ilev,iv)-&
!                                                                           scalar_delhp_full_level_np1(ilev,iv)/scalar_delhp_full_level_n(ilev,iv)+one)
!             delta_pt_b(iv)  = gama*scalar_pressure_full_level_n(ilev-1,iv)*(scalar_potential_temp_full_level_np1(ilev-1,iv)/scalar_potential_temp_full_level_n(ilev-1,iv)-&
!                                                                             scalar_delhp_full_level_np1(ilev-1,iv)/scalar_delhp_full_level_n(ilev-1,iv)+one)
! #endif

            delta_pt_s    = delta_pt_a_s-delta_pt_b_s
#if !defined NEWLGE && !defined LGE1
! old lge beg-----
! #ifdef NHDWPIER
!             phi1_km1(iv)    = scalar_phi_face_level_n(ilev-1,iv)-dtime*scalar_et_phi_face_level_rk(ilev-1,iv)+gravity*dtime*scalar_www_face_level_rk(ilev-1,iv)
!             phi1_k(iv)      = scalar_phi_face_level_n(ilev,iv)  -dtime*scalar_et_phi_face_level_rk(ilev,iv)  +gravity*dtime*scalar_www_face_level_rk(ilev,iv)
! #elif defined (NHDWPIEN)
!             phi1_km1(iv)    = scalar_phi_face_level_n(ilev-1,iv)-dtime*scalar_et_phi_face_level_rk(ilev-1,iv)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev-1,iv)+&
!                                                                                                               gravity*dtime*imbeta*scalar_www_face_level_n(ilev-1,iv)
!             phi1_k(iv)      = scalar_phi_face_level_n(ilev,iv)  -dtime*scalar_et_phi_face_level_rk(ilev,iv)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev,iv)+&
!                                                                                                             gravity*dtime*imbeta*scalar_www_face_level_n(ilev,iv)
! #else
            phi1_km1_s    = scalar_phi_face_level_n(ilev-1,iv)-dtime*scalar_et_phi_face_level_rk(ilev-1,iv)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev-1,iv)
            phi1_k_s      = scalar_phi_face_level_n(ilev,iv)  -dtime*scalar_et_phi_face_level_rk(ilev,iv)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev,iv)
! #endif

            if(ilev.ne.nlev)then
! #ifdef NHDWPIER
!                phi1_kp1(iv) = scalar_phi_face_level_n(ilev+1,iv)-dtime*scalar_et_phi_face_level_rk(ilev+1,iv)+gravity*dtime*scalar_www_face_level_rk(ilev+1,iv)
! #elif defined (NHDWPIEN)
!                phi1_kp1(iv) = scalar_phi_face_level_n(ilev+1,iv)-dtime*scalar_et_phi_face_level_rk(ilev+1,iv)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev+1,iv)+&
!                                                                                                               gravity*dtime*imbeta*scalar_www_face_level_n(ilev+1,iv)
! #else
                  phi1_kp1_s = scalar_phi_face_level_n(ilev+1,iv)-dtime*scalar_et_phi_face_level_rk(ilev+1,iv)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev+1,iv)
! #endif
            else
                  phi1_kp1_s = scalar_phi_face_level_n(ilev+1,iv)
            end if

            delta_phi1_s  = phi1_kp1_s-phi1_k_s
            delta_big_a_s = gama*scalar_pressure_full_level_n(ilev,iv)  *delta_phi1_s/(scalar_phi_face_level_n(ilev+1,iv)-scalar_phi_face_level_n(ilev,iv))

            delta_phi1_s  = phi1_k_s-phi1_km1_s
            delta_big_b_s = gama*scalar_pressure_full_level_n(ilev-1,iv)*delta_phi1_s/(scalar_phi_face_level_n(ilev,iv)  -scalar_phi_face_level_n(ilev-1,iv))
! old lgei end------
#endif
! #ifdef LGE1
!             phi1_km1(iv)    = scalar_phi_face_level_n(ilev-1,iv)-dtime*scalar_et_phi_face_level_rk(ilev-1,iv)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev-1,iv)
!             phi1_k(iv)      = scalar_phi_face_level_n(ilev,iv)  -dtime*scalar_et_phi_face_level_rk(ilev,iv)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev,iv)

!             if(ilev.ne.nlev)then
!                phi1_kp1(iv) = scalar_phi_face_level_n(ilev+1,iv)-dtime*scalar_et_phi_face_level_rk(ilev+1,iv)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev+1,iv)
!             else
!                phi1_kp1(iv) = scalar_phi_face_level_n(ilev+1,iv)
!             end if

!             delta_phi1(iv)  = phi1_kp1(iv)-phi1_k(iv)
!             delta_big_a(iv) = gama*scalar_pressure_full_level_n(ilev,iv)*delta_phi1(iv)/(scalar_phi_face_level_n(ilev+1,iv)-scalar_phi_face_level_n(ilev,iv))*&
!                               scalar_delhp_full_level_n(ilev,iv)/scalar_delhp_full_level_np1(ilev,iv)

!             delta_phi1(iv)  = phi1_k(iv)-phi1_km1(iv)
!             delta_big_b(iv) = gama*scalar_pressure_full_level_n(ilev-1,iv)*delta_phi1(iv)/(scalar_phi_face_level_n(ilev,iv)-scalar_phi_face_level_n(ilev-1,iv))*&
!                               scalar_delhp_full_level_n(ilev-1,iv)/scalar_delhp_full_level_np1(ilev-1,iv)
! #endif

            delta_big_s   = delta_big_a_s-delta_big_b_s

            part_b_s      = scalar_mif_face_level(ilev,iv)*((scalar_pressure_full_level_n(ilev,iv)-scalar_pressure_full_level_n(ilev-1,iv))+delta_pt_s-delta_big_s)

            scalar_rrr_face_level(ilev,iv) = part_a_s+gravity*imkesi*dtime*part_b_s
         end do
      end do
!$omp end do nowait
!$omp end parallel 
!$omp end target
      call t_stopf("nh_compute_rrr_target")
#else

      do ilev = 2, nlev ! each face level

         part_a      = scalar_delhp_face_level_np1(ilev,:)*&
! #ifdef NHDWPIER
!                        (scalar_www_face_level_n(ilev,:)-scalar_www_face_level_rk(ilev,:)-dtime*scalar_et_www_face_level_rk(ilev,:)-&
! #elif defined (NHDWPIEN)
!                        (-dtime*scalar_et_www_face_level_rk(ilev,:)-&
! #else
                       (scalar_www_face_level_n(ilev,:)-dtime*scalar_et_www_face_level_rk(ilev,:)-&
! #endif
                        gravity*dtime+&
                        scalar_mif_face_level(ilev,:)*gravity*(1._r8-imkesi)*dtime*scalar_delp_face_level_rk(ilev,:)/scalar_delhp_face_level_rk(ilev,:))

#if !defined NEWLGE && !defined LGE1
         delta_pt_a  = gama*scalar_pressure_full_level_n(ilev,:)*(scalar_potential_temp_full_level_np1(ilev,:)/scalar_potential_temp_full_level_n(ilev,:))
         delta_pt_b  = gama*scalar_pressure_full_level_n(ilev-1,:)*(scalar_potential_temp_full_level_np1(ilev-1,:)/scalar_potential_temp_full_level_n(ilev-1,:))
#endif
! #ifdef LGE1
!          delta_pt_a  = gama*scalar_pressure_full_level_n(ilev,:)*(scalar_potential_temp_full_level_np1(ilev,:)/scalar_potential_temp_full_level_n(ilev,:)-&
!                                                                   scalar_delhp_full_level_np1(ilev,:)/scalar_delhp_full_level_n(ilev,:)+one)
!          delta_pt_b  = gama*scalar_pressure_full_level_n(ilev-1,:)*(scalar_potential_temp_full_level_np1(ilev-1,:)/scalar_potential_temp_full_level_n(ilev-1,:)-&
!                                                                     scalar_delhp_full_level_np1(ilev-1,:)/scalar_delhp_full_level_n(ilev-1,:)+one)
! #endif

         delta_pt    = delta_pt_a-delta_pt_b
#if !defined NEWLGE && !defined LGE1
! old lge beg-----
! #ifdef NHDWPIER
!          phi1_km1    = scalar_phi_face_level_n(ilev-1,:)-dtime*scalar_et_phi_face_level_rk(ilev-1,:)+gravity*dtime*scalar_www_face_level_rk(ilev-1,:)
!          phi1_k      = scalar_phi_face_level_n(ilev,:)  -dtime*scalar_et_phi_face_level_rk(ilev,:)  +gravity*dtime*scalar_www_face_level_rk(ilev,:)
! #elif defined (NHDWPIEN)
!          phi1_km1    = scalar_phi_face_level_n(ilev-1,:)-dtime*scalar_et_phi_face_level_rk(ilev-1,:)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev-1,:)+&
!                                                                                                      gravity*dtime*imbeta*scalar_www_face_level_n(ilev-1,:)
!          phi1_k      = scalar_phi_face_level_n(ilev,:)  -dtime*scalar_et_phi_face_level_rk(ilev,:)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev,:)+&
!                                                                                                    gravity*dtime*imbeta*scalar_www_face_level_n(ilev,:)
! #else
         phi1_km1    = scalar_phi_face_level_n(ilev-1,:)-dtime*scalar_et_phi_face_level_rk(ilev-1,:)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev-1,:)
         phi1_k      = scalar_phi_face_level_n(ilev,:)  -dtime*scalar_et_phi_face_level_rk(ilev,:)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev,:)
! #endif

         if(ilev.ne.nlev)then
! #ifdef NHDWPIER
!             phi1_kp1 = scalar_phi_face_level_n(ilev+1,:)-dtime*scalar_et_phi_face_level_rk(ilev+1,:)+gravity*dtime*scalar_www_face_level_rk(ilev+1,:)
! #elif defined (NHDWPIEN)
!             phi1_kp1 = scalar_phi_face_level_n(ilev+1,:)-dtime*scalar_et_phi_face_level_rk(ilev+1,:)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev+1,:)+&
!                                                                                                      gravity*dtime*imbeta*scalar_www_face_level_n(ilev+1,:)
! #else
!             phi1_kp1 = scalar_phi_face_level_n(ilev+1,:)-dtime*scalar_et_phi_face_level_rk(ilev+1,:)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev+1,:)
! #endif
         else
            phi1_kp1 = scalar_phi_face_level_n(ilev+1,:)
         end if

         delta_phi1  = phi1_kp1-phi1_k
         delta_big_a = gama*scalar_pressure_full_level_n(ilev,:)  *delta_phi1/(scalar_phi_face_level_n(ilev+1,:)-scalar_phi_face_level_n(ilev,:))

         delta_phi1  = phi1_k-phi1_km1
         delta_big_b = gama*scalar_pressure_full_level_n(ilev-1,:)*delta_phi1/(scalar_phi_face_level_n(ilev,:)  -scalar_phi_face_level_n(ilev-1,:))
! old lgei end------
#endif
! #ifdef LGE1
!          phi1_km1    = scalar_phi_face_level_n(ilev-1,:)-dtime*scalar_et_phi_face_level_rk(ilev-1,:)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev-1,:)
!          phi1_k      = scalar_phi_face_level_n(ilev,:)  -dtime*scalar_et_phi_face_level_rk(ilev,:)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev,:)

!          if(ilev.ne.nlev)then
!             phi1_kp1 = scalar_phi_face_level_n(ilev+1,:)-dtime*scalar_et_phi_face_level_rk(ilev+1,:)+gravity*dtime*(1._r8-imbeta)*scalar_www_face_level_rk(ilev+1,:)
!          else
!             phi1_kp1 = scalar_phi_face_level_n(ilev+1,:)
!          end if

!          delta_phi1  = phi1_kp1-phi1_k
!          delta_big_a = gama*scalar_pressure_full_level_n(ilev,:)*delta_phi1/(scalar_phi_face_level_n(ilev+1,:)-scalar_phi_face_level_n(ilev,:))*&
!                        scalar_delhp_full_level_n(ilev,:)/scalar_delhp_full_level_np1(ilev,:)

!          delta_phi1  = phi1_k-phi1_km1
!          delta_big_b = gama*scalar_pressure_full_level_n(ilev-1,:)*delta_phi1/(scalar_phi_face_level_n(ilev,:)-scalar_phi_face_level_n(ilev-1,:))*&
!                        scalar_delhp_full_level_n(ilev-1,:)/scalar_delhp_full_level_np1(ilev-1,:)
! #endif

         delta_big   = delta_big_a-delta_big_b

         part_b      = scalar_mif_face_level(ilev,:)*((scalar_pressure_full_level_n(ilev,:)-scalar_pressure_full_level_n(ilev-1,:))+delta_pt-delta_big)

         scalar_rrr_face_level(ilev,:) = part_a+gravity*imkesi*dtime*part_b
      end do
#endif
call t_startf("nh_compute_rrr_copy")
!$omp target 
!$omp parallel workshare
      scalar_rrr_face_level(1,:)      = fillvalue
      scalar_rrr_face_level(nlev+1,:) = fillvalue
!$omp end parallel workshare
!$omp end target
      call t_stopf("nh_compute_rrr_copy")
    return
   end subroutine grist_nh_compute_rrr

   subroutine grist_nh_compute_www(mesh, nlev,&
                                   scalar_aaa_full_level, &
                                   scalar_bbb_face_level, &
                                   scalar_rrr_face_level, &
                                   www_top, www_bot  , &
                                   scalar_www_face_level)
! io
     use omp_lib
     type(global_domain),   intent(in)    :: mesh
     integer(i4),           intent(in)    :: nlev
     real(r8),              intent(in)    :: scalar_aaa_full_level(nlev ,mesh%nv_full)
     real(r8),              intent(in)    :: scalar_bbb_face_level(nlev+1,mesh%nv_full)
     real(r8),              intent(in)    :: scalar_rrr_face_level(nlev+1,mesh%nv_full)
     real(r8),              intent(in)    :: www_top(mesh%nv_full), www_bot(mesh%nv_full)
     real(r8), allocatable, intent(inout) :: scalar_www_face_level(:,:)
! local
     integer(i4)                          :: ilev, nsize,iv
     real(r8), allocatable                :: scalar_aaa_face_level(:)
     real(r8), allocatable                :: scalar_bbb_fcut_level(:)
     real(r8), allocatable                :: scalar_ccc_face_level(:)
     real(r8), allocatable                :: scalar_rrr_fcut_level(:)

       nsize = nlev-1

       if(.not.allocated(scalar_aaa_face_level)) allocate(scalar_aaa_face_level(nsize))
       if(.not.allocated(scalar_bbb_fcut_level)) allocate(scalar_bbb_fcut_level(nsize))
       if(.not.allocated(scalar_ccc_face_level)) allocate(scalar_ccc_face_level(nsize))
       if(.not.allocated(scalar_rrr_fcut_level)) allocate(scalar_rrr_fcut_level(nsize))

!$omp parallel  private(iv,scalar_aaa_face_level,scalar_ccc_face_level,scalar_bbb_fcut_level,scalar_rrr_fcut_level) 
!$omp do schedule(dynamic,20)
       do iv = 1, mesh%nv_halo(1)
       scalar_aaa_face_level(1)         = 0._r8
       scalar_aaa_face_level(2:nsize)   = scalar_aaa_full_level(2:nsize,iv)
       scalar_ccc_face_level(1:nsize-1) = scalar_aaa_full_level(2:nsize,iv)
       scalar_ccc_face_level(nsize)     = 0._r8

       scalar_bbb_fcut_level            = scalar_bbb_face_level(2:nsize+1,iv)
       scalar_rrr_fcut_level            = scalar_rrr_face_level(2:nsize+1,iv)

       scalar_rrr_fcut_level(1)         = scalar_rrr_fcut_level(1)-scalar_aaa_full_level(1,iv)*www_top(iv)

       call solve_tridiagonal(scalar_aaa_face_level, scalar_bbb_fcut_level, scalar_ccc_face_level, &
                              scalar_rrr_fcut_level, scalar_www_face_level(2:nsize+1,iv), nsize)

       scalar_www_face_level(1,iv)        = www_top(iv)
       scalar_www_face_level(nlev+1,iv)   = www_bot(iv)

       end do
!$omp end do nowait
!$omp end parallel 

       deallocate(scalar_aaa_face_level)
       deallocate(scalar_bbb_fcut_level)
       deallocate(scalar_ccc_face_level)
       deallocate(scalar_rrr_fcut_level)

      return
   end subroutine grist_nh_compute_www
  end module grist_nh_implicit_module_2d_mixed
