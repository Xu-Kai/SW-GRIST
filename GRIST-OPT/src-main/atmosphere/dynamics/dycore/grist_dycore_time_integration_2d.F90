
!----------------------------------------------------------------------------
! Copyright (c), GRIST-Dev
!
! Unless noted otherwise source code is licensed under the Apache-2.0 license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://https://github.com/grist-dev
!
! Description: This is the major time stepping driver for the 3d model in
!              either, hdc or ndc configuration
!
! Revision history:
!----------------------------------------------------------------------------

  module grist_dycore_time_integration_2d

    use grist_constants,       only: rvap,rdry,ptfactor,cp,p00,zero, i4, r8, gravity, one, rearth
    use grist_nml_module,      only: nlev, mas_adv_flag, pot_adv_flag, ver_adv_flag, nrk, hor_pgf_flag, &
                                     nh_dynamics, gcm_testcase,  nsteps, &
                                     ad_dycore_laplacian_2nd, ad_dycore_laplacian_4th, ad_dycore_laplacian_6th, &
                                     tend_nct_once, working_mode, use_phys, physpkg, &
                                     ptendSubDiabPhys, ptend_wind_rk_on, ptend_heat_rk_on, ptend_dycore_f1_on, ptend_dycore_heat_f1_on, &
                                     ptend_dycore_wind_f1_on, write_stepinfo, write_verbose, use_www_hyperDiffusion, &
                                     restore_hydro, restore_hydro_minsteps,restore_hydro_intsteps, model_timestep, adjphi
#ifdef AMIPW_PHYSICS
    use grist_wrfphys_nml_module, only: step_cu
#endif
    use grist_domain_types,    only: global_domain
    use grist_data_types,      only: scalar_1d_field, scalar_2d_field, exchange_field_list_2d, exchange_field_list_1d
!
! horizontal module
!
    use grist_dycore_primal_flux_operators_2d, only: calc_primal_normal_flux_at_edge
    use grist_dycore_hori_swe_module_2d,       only: calc_tend_nct_at_edge_full_level, calc_grad_kinetic_energy
    use grist_dycore_gcd_recon_module_2d,      only: gradient_operator_2d, &
                                                     gradient_operator_2d_var3, &
                                                     gradient_operator_2d_var4, &
                                                     divergence_operator_2d
    use grist_dycore_ref_atmos
!
! vertical qh
!
    use grist_hpe_hydro_pgf,          only: calc_hpe_hydro               , &
                                            calc_hpe_hpressure_face_level , &
                                            calc_hpe_hpressure_full_level , &
                                            calc_hpe_delhp, &
                                            calc_hpe_get_full_mass

    use grist_hpe_vertical_advection, only: calc_hpe_vert_advection, &
                                            calc_hpe_tend_vert_mass_flux, &
                                            calc_hpe_tend_vert_mass_flux_2d
    use grist_hpe_continuity        , only: calc_hpe_tend_continuity, calc_hpe_tend_continuity_2d
    use grist_hpe_constants,          only: deta_full, deta_face,nlevp,eta_full,eta_full_b
! nh_dynamics
    use grist_nh_driver_module,       only: grist_nh_dynamics_run, ndc_restore_flag
! dycore data
    use grist_dycore_vars_module,     only: dycoreVarCellFace, dycoreVarCellFull, dycoreVarEdgeFull, dycoreVarEdgeFace, dycoreVarSurface
! tracer
    use grist_tracer_transport_vars_module, only: tracerVarCellFace, tracerVarEdgeFull, tracerVarCellFull
!
! physics
! 
    use grist_physics_data_structure, only: ptend_f1, ptend_rk, ptend_f2
! mpi
    use grist_mpi
#ifndef SEQ_GRIST
    use grist_config_partition,       only: debug_data_1d,&
                                            debug_data_2d,&
                                            exchange_data_2d_add , &
                                            exchange_data_2d,&
                                            exchange_data_1d_add,&
                                            exchange_data_1d
    use grist_clocks,                 only: clock_id, clock_begin, clock_end
#endif

#ifdef LAM_DOMAIN
    use grist_datam_glam_data_module, only: update_assignValueArea_1d, update_assignValueArea_2d
#endif

    implicit none

    private

    public :: dycore_time_integration_init,&
              dycore_time_integration_run ,&
              dycore_time_integration_final
!
! local varitable for template io of subroutine
!
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlev_a
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlev_b
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlev_c
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlev_d
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlevp_a
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlevp_b
     real(r8), dimension(:), allocatable :: scalar_template_1d_nlevp_c
     real(r8)                            :: scalar_template_a
     type(scalar_2d_field)               :: scalar_template_2d_ne_b
     type(scalar_2d_field)               :: scalar_pressure_at_pc_face_level_rk  ! nh-diag-var
     type(scalar_2d_field)               :: scalar_www_at_pc_face_level_np1      ! nh-prog-var
     type(scalar_2d_field)               :: scalar_phi_at_pc_face_level_np1      ! nh-prog-var
     type(scalar_2d_field)               :: scalar_hpressure_at_pc_face_level_np1
     type(scalar_2d_field)               :: scalar_delhp_at_pc_face_level_np1
     type(scalar_2d_field)               :: scalar_delhp_at_pc_face_level_rk
     type(scalar_2d_field)               :: scalar_hpressure_at_pc_full_level_np1
     integer(i4) :: clock_rkl, clock_mas, clock_nct, clock_vau, clock_ket, clock_diagom, clock_nhd, clock_ptm, clock_pgf, clock_mainexch, clock_fnlrk

   CONTAINS

!-----------------------------------------------------------------------------
! Name convention:
! for prognostic variables
! q(n+1) = q(n)+F(rk)
! n+1, new state after each sub stage
! n, state at begining of each time step
! rk itermediate state used for explicit tendency computin
!
! for other variable, _n denotes an intermediate state
! np1 denotes next (if needed)
!------------------------------------------------------------------------------

   subroutine dycore_time_integration_run(mesh, dtime, itimestep,idstep,dstep_in_tstep,itstep,tstep_in_mstep)
!
! io
!
     use omp_lib
     use omnicopy_mod
     type(global_domain),  intent(inout) :: mesh
     real(r8)          ,   intent(in)    :: dtime
     integer(i4)       ,   intent(in)    :: itimestep
     integer(i4), optional,intent(in)    :: idstep, dstep_in_tstep, itstep, tstep_in_mstep
!
! variables used within RK step, local rk
!

     type(scalar_2d_field)               :: scalar_www_at_pc_face_level_rk
     type(scalar_2d_field)               :: scalar_phi_at_pc_face_level_rk
     type(scalar_2d_field)               :: scalar_normal_velocity_at_edge_full_level_rk
     type(scalar_2d_field)               :: scalar_delhp_at_pc_full_level_rk
! local np1 vars
     type(scalar_1d_field)               :: scalar_hpressure_at_pc_surface_np1
     type(scalar_2d_field)               :: scalar_normal_velocity_at_edge_full_level_np1
     type(scalar_2d_field)               :: scalar_mass_pt_at_pc_full_level_np1
! local
     real(r8)                            :: rk_substep, fb_step
     real(r8)                            :: cr_num, cr_act, cr_sum, eta, kv 
     real(r8)                            :: div_sum(nlev), tmp1(nlev), tmp2(nlev), tmp3(nlev), tmp4(nlev)
     real(r8)                            :: tend_pidpiv_at_pc_full_level_hori(nlev,mesh%nv_full)
     real(r8)                            :: scalar_mpressure_at_edge_full_level_n(nlev,mesh%ne_full)
     real(r8)                            :: scalar_mpressure_at_pc_full_level_old(nlev,mesh%nv_full)
#ifdef AMIPW_PHYSICS
     real(r8)                            :: scalar_pt_at_pc_full_level_old(nlev,mesh%nv_full)
#endif
     integer(i4)                         :: nsplit                 ! how many small steps within DT
     integer(i4)                         :: rk_number              ! number to divide DT in RK step
     integer(i4)                         :: irk_step
     integer(i4)                         :: ifb_step
     integer(i4)                         :: it, ie, iv, ilev, inb
     integer(i4)                         :: icell1, icell2, index_edge
     integer(i4)                         :: iblock
! local
     real(r8) :: scalar_template_1d_nlevp_d(nlevp)
     real(r8) :: scalar_template_1d_nlevp_e(nlevp)
     real(r8) :: scalar_template_1d_nlevp_f(nlevp)
     real(r8) :: scalar_normal_velocity_n_c(nlev), &
                     tend_normal_velocity_vadv_c(nlev), &
                     tend_normal_velocity_pgf_c(nlev), &
                     tend_normal_velocity_gz_c(nlev), &
                     tend_normal_velocity_ke_c(nlev), &
                     tend_normal_velocity_nct_c(nlev), &
                     tmp1_c(nlev), tmp2_c(nlev)
#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
#include "data_check.inc"
        real(r8), allocatable :: ref_scalar_normal_mass_flux_n_f(:, :)
        real(r8), allocatable :: ref_tend_mass_hori_f(:, :)
        real(r8), allocatable :: ref_tend_normal_velocity_vadv_f(:,:)
        real(r8), allocatable :: ref_scalar_mpressure_at_edge_full_level_n_f(:,:)
        real(r8), allocatable :: ref_tend_mass_pt_hori_f(:, :)
        real(r8), allocatable :: ref_scalar_mpressure_at_edge_full_level_n(:, :)
        real(r8), allocatable :: ref_tend_pidpiv_at_pc_full_level_hori(:, :)
        real(r8), allocatable :: ref_tend_normal_velocity_pgf_f(:, :)
        real(r8), allocatable :: ref_scalar_omega_n_f(:,:)
#endif
      type(exchange_field_list_2d),pointer :: field_head_2d
      type(exchange_field_list_1d),pointer :: field_head_1d
      integer(i4) :: itv
      field_head_2d=>null()
      field_head_1d=>null()
      iblock = mpi_rank()
      !iblock = -1
!================================================
! added for cumulus wrfphys
!================================================

#ifdef AMIPW_PHYSICS
      IF(present(idstep))then
! only at the first dycore and tracer step
         if(idstep.eq.1.and.itstep.eq.1.and.mod(itimestep-1,step_cu).eq.0) dycoreVarCellFull%tend_pt_n%f = zero
      END IF
#endif

!
! prepare for RK dyn integration
!
      if(nh_dynamics)then 
         call t_startf("dycore_time_integration_run_nh_dynamics_copy")
         if(.not.allocated(scalar_www_at_pc_face_level_rk%f))    allocate(scalar_www_at_pc_face_level_rk%f(1:ubound(dycoreVarCellFace%scalar_www_n%f, 1),1:ubound(dycoreVarCellFace%scalar_www_n%f, 2)))
         if(.not.allocated(scalar_phi_at_pc_face_level_rk%f))    allocate(scalar_phi_at_pc_face_level_rk%f(1:ubound(dycoreVarCellFace%scalar_phi_n%f, 1),1:ubound(dycoreVarCellFace%scalar_phi_n%f, 2)))
!$omp target
!$omp parallel workshare         
         scalar_www_at_pc_face_level_rk%f(:,:)      = dycoreVarCellFace%scalar_www_n%f(:,:) 
         scalar_phi_at_pc_face_level_rk%f(:,:)     = dycoreVarCellFace%scalar_phi_n%f(:,:)
         dycoreVarCellFull%scalar_pressure_rk%f(:,:) = dycoreVarCellFull%scalar_pressure_n%f(:,:)  
!$omp end parallel workshare
!$omp end target 
         scalar_www_at_pc_face_level_rk%pos      = dycoreVarCellFace%scalar_www_n%pos
         scalar_phi_at_pc_face_level_rk%pos     = dycoreVarCellFace%scalar_phi_n%pos
         dycoreVarCellFull%scalar_pressure_rk%pos= dycoreVarCellFull%scalar_pressure_n%pos
         call t_stopf("dycore_time_integration_run_nh_dynamics_copy")
      end if
      call t_startf("dycore_time_integration_run_init_copy")
      if(.not.allocated(scalar_normal_velocity_at_edge_full_level_rk%f))    allocate(scalar_normal_velocity_at_edge_full_level_rk%f(1:ubound(dycoreVarEdgeFull%scalar_normal_velocity_n%f, 1),1:ubound(dycoreVarEdgeFull%scalar_normal_velocity_n%f, 2)))
      if(.not.allocated(scalar_delhp_at_pc_full_level_rk%f))    allocate(scalar_delhp_at_pc_full_level_rk%f(1:ubound(dycoreVarCellFull%scalar_delhp_n%f, 1),1:ubound(dycoreVarCellFull%scalar_delhp_n%f, 2)))
      if(.not.allocated(scalar_delhp_at_pc_face_level_rk%f))    allocate(scalar_delhp_at_pc_face_level_rk%f(1:ubound(dycoreVarCellFace%scalar_delhp_n%f, 1),1:ubound(dycoreVarCellFace%scalar_delhp_n%f, 2)))
!$omp target
!$omp parallel workshare
      scalar_normal_velocity_at_edge_full_level_rk%f(:,:) = dycoreVarEdgeFull%scalar_normal_velocity_n%f(:,:)
      scalar_delhp_at_pc_full_level_rk%f(:,:)             = dycoreVarCellFull%scalar_delhp_n%f(:,:)
      scalar_delhp_at_pc_face_level_rk%f(:,:)             = dycoreVarCellFace%scalar_delhp_n%f(:,:)
      scalar_mpressure_at_pc_full_level_old(:,:)        = dycoreVarCellFull%scalar_mpressure_n%f(:,:)
!$omp end parallel workshare
!$omp end target
      scalar_normal_velocity_at_edge_full_level_rk%pos = dycoreVarEdgeFull%scalar_normal_velocity_n%pos
      scalar_delhp_at_pc_full_level_rk%pos             = dycoreVarCellFull%scalar_delhp_n%pos
      scalar_delhp_at_pc_face_level_rk%pos             = dycoreVarCellFace%scalar_delhp_n%pos
      call t_stopf("dycore_time_integration_run_init_copy")
#ifdef AMIPW_PHYSICS
! raw pt
      scalar_pt_at_pc_full_level_old = dycoreVarCellFull%scalar_potential_temp_n%f/(one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:))
#endif
!
! BEGIN THIS RK STEP, use zzz to represent code that should be sleep foreever
!
! we have commented if statement assocaited stencil_exchange_flag,
! stencil_exchange_flag=.true. is default now
!
#ifndef SEQ_GRIST
     call clock_begin(clock_rkl)
#endif

     DO irk_step = 1, nrk

        rk_number  = nrk+1-irk_step
        rk_substep = dtime/rk_number

!-----------------------------------------------------------------------
! compute normal mass flux {DELP*V} and mass flux div {GRAD.(DELP*V)} 
! at each full level
!-----------------------------------------------------------------------
#ifndef SEQ_GRIST
        call clock_begin(clock_mas)
#endif

        mesh%ne = mesh%ne_halo(1)
! inline_opt_here

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
allocate(ref_scalar_normal_mass_flux_n_f, source=dycoreVarEdgeFull%scalar_normal_mass_flux_n%f)
#endif
call t_startf("dycore_time_intergration_target_1")
!$omp target
!$omp parallel  private(ie,ilev)
!$omp do
        do ie = 1, mesh%ne_halo(1)
           do ilev = 1, nlev
              dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev, ie) = 0.5*(scalar_delhp_at_pc_full_level_rk%f(ilev,mesh%edt_v(1,ie))+&
                                                                             scalar_delhp_at_pc_full_level_rk%f(ilev,mesh%edt_v(2,ie)))
           end do
        end do

!$omp end do nowait
!$omp end parallel
!$omp end target 
call t_stopf("dycore_time_intergration_target_1")

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
call t_startf("dycore_time_intergration_1")
!$omp parallel  private(ie,ilev)
!$omp do
      do ie = 1, mesh%ne_halo(1)
        do ilev = 1, nlev
           ref_scalar_normal_mass_flux_n_f(ilev, ie) = 0.5*(scalar_delhp_at_pc_full_level_rk%f(ilev,mesh%edt_v(1,ie))+&
                                                                          scalar_delhp_at_pc_full_level_rk%f(ilev,mesh%edt_v(2,ie)))
        end do
     end do

!$omp end do nowait
!$omp end parallel
call t_stopf("dycore_time_intergration_1")
call data_check(ref_scalar_normal_mass_flux_n_f, dycoreVarEdgeFull%scalar_normal_mass_flux_n%f)
deallocate(ref_scalar_normal_mass_flux_n_f)
#endif
! inline_opt_here
        call t_startf("dycore_time_integration_run_scalar_normal_mass_flux_n_copy")
        !$omp target
        !$omp parallel workshare
        dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(:,:) = dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(:,:)*&
                                                        scalar_normal_velocity_at_edge_full_level_rk%f(:,:)
        !$omp end parallel workshare
        !$omp end target
        call t_stopf("dycore_time_integration_run_scalar_normal_mass_flux_n_copy")
!================================================
! exchange scheme, needed if o3 is used, but o2
! is used so hard coded here!
!        call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_mass_flux_n)
!        call exchange_data_2d(mesh%local_block,field_head_2d)
!        mesh%ne = mesh%ne_halo(1)
!================================================
        ! from ne_halo(1) to nv_halo(1)
        mesh%nv = mesh%nv_halo(1)
#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
allocate(ref_tend_mass_hori_f, source= dycoreVarCellFull%tend_mass_hori%f)
#endif

call t_startf("dycore_time_intergration_target_2")
!$omp target
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev)
!$omp do
        do iv = 1, mesh%nv
           div_sum(1:nlev) = zero
           do inb = 1, mesh%vtx_nnb(iv)
              index_edge  = mesh%vtx_ed(inb,iv)
              do ilev = 1, nlev
                 div_sum(ilev)  =  div_sum(ilev)+dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
              end do
           end do
           dycoreVarCellFull%tend_mass_hori%f(1:nlev,iv) = div_sum(1:nlev)/((rearth**2)*mesh%plg_areag(iv))
        end do
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("dycore_time_intergration_target_2")

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
call t_startf("dycore_time_intergration_2")
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev)
!$omp do
        do iv = 1, mesh%nv
                div_sum(1:nlev) = zero
                do inb = 1, mesh%vtx_nnb(iv)
                index_edge  = mesh%vtx_ed(inb,iv)
                do ilev = 1, nlev
                        div_sum(ilev)  =  div_sum(ilev)+dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
                end do
                end do
                ref_tend_mass_hori_f(1:nlev,iv) = div_sum(1:nlev)/((rearth**2)*mesh%plg_areag(iv))
        end do
!$omp end do nowait
!$omp end parallel
call data_check( dycoreVarCellFull%tend_mass_hori%f, ref_tend_mass_hori_f)
deallocate(ref_tend_mass_hori_f)
call t_stopf("dycore_time_intergration_2")
#endif
! inline_opt_here

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute normal mass flux",mas_adv_flag

        mesh%ne = mesh%ne_compute
        mesh%nv = mesh%nv_compute
!
! compute PS tendency and face level vertical mass flux (m*etadot)
!
        call calc_hpe_tend_continuity_2d(mesh%nv_halo(1), nlev                      , &
                                         dycoreVarCellFull%tend_mass_hori%f         , & ! in
                                         dycoreVarSurface%tend_hpressure_cnty%f     , & ! out
                                         dycoreVarCellFace%scalar_eta_mass_flux_n%f)    ! out
        call t_startf("dycore_time_integration_run_scalar_eta_mass_flux_n_op")
!$omp target
!$omp parallel workshare
        dycoreVarCellFull%scalar_eta_mass_flux_n%f(1:nlev,:) = 0.5_r8*dycoreVarCellFace%scalar_eta_mass_flux_n%f(1:nlev,:)+&
                                                               0.5_r8*dycoreVarCellFace%scalar_eta_mass_flux_n%f(2:nlev+1,:)
!$omp end parallel workshare
!$omp end target
        call t_stopf("dycore_time_integration_run_scalar_eta_mass_flux_n_op")
        if(iblock .eq. 0.and.write_verbose) print*,"finish compute PS and face level V mass flux",nh_dynamics

!================================================
! diag FACE LEVEL VARS
!================================================
#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
allocate(ref_scalar_normal_mass_flux_n_f, source= dycoreVarEdgeFace%scalar_normal_mass_flux_n%f)
#endif
call t_startf("dycore_time_intergration_target_3")

!$omp target
!$omp parallel  private(itv, ilev)
!$omp do
      do itv = 1, ubound(dycoreVarEdgeFace%scalar_normal_mass_flux_n%f, 2)  
        do ilev = 2, nlev
             dycoreVarEdgeFace%scalar_normal_mass_flux_n%f(ilev,itv) = &
                                     0.5*(deta_full(ilev-1)/deta_face(ilev)*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev,itv)+&
                                          deta_full(ilev)  /deta_face(ilev)*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev-1,itv))
        end do
     enddo
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("dycore_time_intergration_target_3")

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)

call t_startf("dycore_time_intergration_3")
!$omp parallel  private(ilev)
!$omp do
        do ilev = 2, nlev
                ref_scalar_normal_mass_flux_n_f(ilev,:) = &
                                        0.5*(deta_full(ilev-1)/deta_face(ilev)*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev,:)+&
                                                deta_full(ilev)  /deta_face(ilev)*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(ilev-1,:))
        end do
!$omp end do nowait
!$omp end parallel
call t_stopf("dycore_time_intergration_3")

call data_check(dycoreVarEdgeFace%scalar_normal_mass_flux_n%f, ref_scalar_normal_mass_flux_n_f)
deallocate(ref_scalar_normal_mass_flux_n_f)

#endif
        call t_startf("dycore_time_integration_run_edgeface_scalar_eta_mass_flux_n_op")
!$omp target
!$omp parallel workshare
        dycoreVarEdgeFace%scalar_normal_mass_flux_n%f(1,:)     = 0.5_r8*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(1,:)
        dycoreVarEdgeFace%scalar_normal_mass_flux_n%f(nlevp,:) = 0.5_r8*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(nlev,:)
!$omp end parallel workshare
!$omp end target
        call t_stopf("dycore_time_integration_run_edgeface_scalar_eta_mass_flux_n_op")
        mesh%nv = mesh%nv_halo(1)  ! only used by nh_dynamics originally, now also used for diagnosing dpidry/dt
        call divergence_operator_2d(mesh,dycoreVarEdgeFace%scalar_normal_mass_flux_n%f, &
                                         dycoreVarCellFace%tend_mass_hori%f, & ! face-level adv
                                         nlev+1)
        mesh%nv = mesh%nv_compute
        call t_startf("dycore_time_integration_run_scalar_hpress_pc_surface_np1_op")
        !!$omp target
        !!$omp parallel workshare
        scalar_hpressure_at_pc_surface_np1%f = dycoreVarSurface%scalar_hpressure_n%f+&
                                               dtime/rk_number*dycoreVarSurface%tend_hpressure_cnty%f
        !!$omp end parallel workshare
        !!$omp end target
        call t_stopf("dycore_time_integration_run_scalar_hpress_pc_surface_np1_op")
#ifndef SEQ_GRIST
        call exchange_data_1d_add(mesh,field_head_1d,scalar_hpressure_at_pc_surface_np1)
        call exchange_data_1d(mesh%local_block,field_head_1d)
#endif
#ifdef LAM_DOMAIN
        call update_assignValueArea_1d(mesh,"hpressure_surface",scalar_hpressure_at_pc_surface_np1)
#endif
!
! renew mass/hpressure state
!
        mesh%nv = mesh%nv_full
        call time_integration_renew_mass_state(mesh, scalar_hpressure_at_pc_surface_np1   ,&
                                                     scalar_hpressure_at_pc_face_level_np1,& ! overwritten each RK step
                                                     dycoreVarCellFull%scalar_delhp_np1   ,& ! overwritten each RK step
                                                     scalar_hpressure_at_pc_full_level_np1,& ! overwritten each RK step
                                                     scalar_delhp_at_pc_face_level_np1  )    ! overwritten each RK step

! get full mass, if dycore, just hpres
        call calc_hpe_get_full_mass(nlev, mesh%nv_halo(1), working_mode            , & ! in
                                          scalar_hpressure_at_pc_face_level_np1%f  , & ! in
                                          scalar_hpressure_at_pc_full_level_np1%f  , & ! in
                                          dycoreVarCellFull%scalar_delhp_np1%f     , & ! in
                                          tracerVarCellFull%scalar_mif_n%f         , & ! in
                                          dycoreVarCellFull%scalar_mpressure_n%f   , & ! out
                                          dycoreVarCellFace%scalar_mpressure_n%f)      ! out
#ifndef SEQ_GRIST
        call clock_end(clock_mas)
#endif

!--------------------------------------------------------------------------
! compute horizontal potential temperature (pt) flux and pt mass flux div
! it computes pt-flux until ne_compute, but we need ne_halo(1), so do excg
! for LAM, because at bdy we use 2nd-order flux, so we can compute until
! ne_halo(1), but only inside this routine
!--------------------------------------------------------------------------

!
! begin pt(m)-flux computing
!
#ifndef SEQ_GRIST
        call clock_begin(clock_ptm)
#endif
        call calc_primal_normal_flux_at_edge(mesh, dycoreVarEdgeFull%scalar_normal_mass_flux_n%f   , &
                                                   dycoreVarEdgeFull%scalar_normal_mass_flux_n%f   , &
                                                   dycoreVarCellFull%scalar_potential_temp_n%f     , &
                                                   dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f, &
                                                   pot_adv_flag(irk_step), rk_substep, nlev)
        call t_startf("dycore_time_integration_run_scalar_template_2d_ne_b_op")
        if (.not. allocated(scalar_template_2d_ne_b%f)) allocate(scalar_template_2d_ne_b%f, mold=dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f)
        !$omp target
        !$omp parallel workshare
        scalar_template_2d_ne_b%f(:,:) = dycoreVarEdgeFull%scalar_normal_pt_mass_flux_n%f(:,:)*&
                                    dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(:,:)
        !$omp end parallel workshare
        !$omp end target
        call t_stopf("dycore_time_integration_run_scalar_template_2d_ne_b_op")
!
! exchange scheme
!
#ifndef SEQ_GRIST
        call exchange_data_2d_add(mesh,field_head_2d,scalar_template_2d_ne_b)
        call exchange_data_2d(mesh%local_block,field_head_2d)
#endif

        mesh%ne = mesh%ne_halo(1)
        mesh%nv = mesh%nv_halo(1)
!inline_opt_here
      !  call divergence_operator_2d(mesh,scalar_template_2d_ne_b%f, &
      !                              dycoreVarCellFull%tend_mass_pt_hori%f, &
      !                              nlev)
#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
allocate(ref_tend_mass_pt_hori_f, source=  dycoreVarCellFull%tend_mass_pt_hori%f)
#endif
call t_startf("dycore_time_intergration_target_4")
!$omp target
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev)
!$omp do
        do iv = 1, mesh%nv
           div_sum(:) = 0._r8
           do inb = 1, mesh%vtx_nnb(iv)
              index_edge  = mesh%vtx_ed(inb,iv)
              do ilev = 1, nlev
                 div_sum(ilev)  = div_sum(ilev)+scalar_template_2d_ne_b%f(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
              end do
           end do
           dycoreVarCellFull%tend_mass_pt_hori%f(:,iv) = div_sum(:)/((rearth**2)*mesh%plg_areag(iv))
        end do
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("dycore_time_intergration_target_4")

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
call t_startf("dycore_time_intergration_4")
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev)
!$omp do
        do iv = 1, mesh%nv
                div_sum(:) = 0._r8
                do inb = 1, mesh%vtx_nnb(iv)
                index_edge  = mesh%vtx_ed(inb,iv)
                do ilev = 1, nlev
                        div_sum(ilev)  = div_sum(ilev)+scalar_template_2d_ne_b%f(ilev,index_edge)*mesh%plg_nr(inb, iv)*rearth*mesh%edp_leng(index_edge)
                end do
                end do
                ref_tend_mass_pt_hori_f(:,iv) = div_sum(:)/((rearth**2)*mesh%plg_areag(iv))
        end do
!$omp end do nowait
!$omp end parallel
call t_stopf("dycore_time_intergration_4")

call data_check( dycoreVarCellFull%tend_mass_pt_hori%f, ref_tend_mass_pt_hori_f)
deallocate(ref_tend_mass_pt_hori_f)
#endif 


!inline_opt_here
        call t_startf("dycore_time_integration_run_tend_mass_pt_hori_fill")
        !$omp target
        !$omp parallel workshare
                dycoreVarCellFull%tend_mass_pt_hori%f(:,:) = -1._r8*dycoreVarCellFull%tend_mass_pt_hori%f(:,:)
        !$omp end parallel workshare
        !$omp end target
        call t_stopf("dycore_time_integration_run_tend_mass_pt_hori_fill")
        if(iblock .eq. 0.and.write_verbose) print*,"finish compute horizontal pt flux"
        mesh%ne = mesh%ne_compute
        mesh%nv = mesh%nv_compute
!
! compute vertical potential temperature mass flux tendency ( partial (m*etadot*theta) )
! 
        call calc_hpe_tend_vert_mass_flux_2d(mesh%nv_full,mesh%nv_halo(1),nlev, & ! in
                                             dycoreVarCellFull%scalar_potential_temp_n%f , & ! in
                                             dycoreVarCellFace%scalar_eta_mass_flux_n%f  , & ! in
                                             ver_adv_flag                               , & ! in
                                             dycoreVarCellFull%tend_mass_pt_vert%f) ! out
!
! mass*pt update
!
        call t_startf("dycore_time_integration_scalar_mass_pt_at_pc_full_level_np1_op1")
        if (.not. allocated(scalar_mass_pt_at_pc_full_level_np1%f)) &
           allocate(scalar_mass_pt_at_pc_full_level_np1%f, mold=dycoreVarCellFull%scalar_mass_pt_n%f)
        !$omp target
        !$omp parallel workshare
        scalar_mass_pt_at_pc_full_level_np1%f(:,:)  = dycoreVarCellFull%scalar_mass_pt_n%f(:,:)+&
                                                 dtime/rk_number*&
                                                (dycoreVarCellFull%tend_mass_pt_hori%f(:,:)+&
                                                 dycoreVarCellFull%tend_mass_pt_vert%f(:,:))
        !$omp end parallel workshare
        !$omp end target
        call t_stopf("dycore_time_integration_scalar_mass_pt_at_pc_full_level_np1_op1")
        if(use_phys.and.ptend_heat_rk_on)then
            call t_startf("dycore_time_integration_scalar_mass_pt_at_pc_full_level_np1_op2")
            !$omp target
            !$omp parallel workshare
            scalar_mass_pt_at_pc_full_level_np1%f(:,:) = scalar_mass_pt_at_pc_full_level_np1%f(:,:)+&
                                                  dtime/rk_number*ptend_rk%tend_mass_pt_at_pc_full_level%f(:,:)
            !$omp end parallel workshare
            !$omp end target
            call t_stopf("dycore_time_integration_scalar_mass_pt_at_pc_full_level_np1_op2")
        end if
        if(ad_dycore_laplacian_2nd)then
            call t_startf("dycore_time_integration_scalar_mass_pt_at_pc_full_level_np1_op3")
            !$omp target
            !$omp parallel workshare
            scalar_mass_pt_at_pc_full_level_np1%f(:,:) = scalar_mass_pt_at_pc_full_level_np1%f(:,:)+&
                                                    dtime/rk_number*dycoreVarCellFull%tend_mass_pt_laplacian_2nd%f(:,:)
            !$omp end parallel workshare
            !$omp end target
            call t_stopf("dycore_time_integration_scalar_mass_pt_at_pc_full_level_np1_op3")
        end if
        call t_startf("dycore_time_integration_scalar_potential_temp_n_op")
        !$omp target
        !$omp parallel workshare
        dycoreVarCellFull%scalar_potential_temp_n%f(:,:) = scalar_mass_pt_at_pc_full_level_np1%f(:,:)/& ! overwritten each RK step
                                                      dycoreVarCellFull%scalar_delhp_np1%f(:,:)
        !$omp end parallel workshare
        !$omp end target
        call t_stopf("dycore_time_integration_scalar_potential_temp_n_op")
#ifndef SEQ_GRIST
        call clock_end(clock_ptm)
#endif

        if(iblock .eq. 0.and.write_verbose) print*,"finish compute vertical pt flux"

!----------------------------------------------------------------------
! compute velocity tendency from Coriolis and Kinetic Energy gradient
!----------------------------------------------------------------------
#ifndef SEQ_GRIST
        call clock_begin(clock_nct)
#endif

        IF((tend_nct_once.and.irk_step.eq.1).or.(.not.tend_nct_once)) &   ! control the frequency of NCT evaluation

        call calc_tend_nct_at_edge_full_level(mesh, scalar_delhp_at_pc_full_level_rk            ,&  ! in
                                                    scalar_normal_velocity_at_edge_full_level_rk,&  ! in
                                                    dycoreVarEdgeFull%scalar_normal_mass_flux_n ,&  ! in
                                                    dycoreVarEdgeFull%tend_normal_velocity_nct  ,&  ! out
                                                    dtime, nlev , irk_step)
#ifndef SEQ_GRIST
        call clock_end(clock_nct)
#endif

        mesh%nv = mesh%nv_halo(1)
        mesh%nt = mesh%nt_halo(1)
!
! KE definition of its gradient
!
#ifndef SEQ_GRIST
        call clock_begin(clock_ket)
#endif
        call calc_grad_kinetic_energy(mesh,scalar_normal_velocity_at_edge_full_level_rk%f,&
                                           dycoreVarEdgeFull%tend_normal_velocity_ke%f  ,&
                                           nlev)
#ifndef SEQ_GRIST
        call clock_end(clock_ket)
#endif

        if(iblock .eq. 0.and.write_verbose) print*,"finish NCT and KE"
#ifndef SEQ_GRIST
        call clock_begin(clock_vau)
#endif

!----------------------------------------------------------
! compute vertical advection of normal velocity at edge
!----------------------------------------------------------

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
allocate(ref_tend_normal_velocity_vadv_f, source=dycoreVarEdgeFull%tend_normal_velocity_vadv%f)
#endif
call t_startf("dycore_time_intergration_target_5")
!$omp target
!$omp parallel private(ie,icell1,icell2,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_a,scalar_template_1d_nlev_b,scalar_template_1d_nlev_c)
!$omp do

        do ie = 1, mesh%ne
!
! center intepolate delhp and mass eta velocity from pc to edge
!
          icell1 = mesh%edt_v(1,ie)
          icell2 = mesh%edt_v(2,ie)
          scalar_template_1d_nlevp_a = &
          0.5_r8*(dycoreVarCellFace%scalar_eta_mass_flux_n%f(:,icell1)+&
                  dycoreVarCellFace%scalar_eta_mass_flux_n%f(:,icell2))
          
          scalar_template_1d_nlev_a  = &
          0.5_r8*(scalar_delhp_at_pc_full_level_rk%f(:,icell1)+&
                  scalar_delhp_at_pc_full_level_rk%f(:,icell2))

          scalar_template_1d_nlev_b(:) = scalar_normal_velocity_at_edge_full_level_rk%f(:,ie)

          call calc_hpe_vert_advection(scalar_template_1d_nlevp_a, &
                                       scalar_template_1d_nlev_a , &
                                       scalar_template_1d_nlev_b , &
                                       scalar_template_1d_nlev_c )

          dycoreVarEdgeFull%tend_normal_velocity_vadv%f(1:nlev,ie) = scalar_template_1d_nlev_c

        end do
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("dycore_time_intergration_target_5")

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
call t_startf("dycore_time_intergration_5")
!$omp parallel private(ie,icell1,icell2,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_a,scalar_template_1d_nlev_b,scalar_template_1d_nlev_c)
!$omp do

        do ie = 1, mesh%ne
!
! center intepolate delhp and mass eta velocity from pc to edge
!
                icell1 = mesh%edt_v(1,ie)
                icell2 = mesh%edt_v(2,ie)
                scalar_template_1d_nlevp_a = &
                0.5_r8*(dycoreVarCellFace%scalar_eta_mass_flux_n%f(:,icell1)+&
                        dycoreVarCellFace%scalar_eta_mass_flux_n%f(:,icell2))
                
                scalar_template_1d_nlev_a  = &
                0.5_r8*(scalar_delhp_at_pc_full_level_rk%f(:,icell1)+&
                        scalar_delhp_at_pc_full_level_rk%f(:,icell2))

                scalar_template_1d_nlev_b(:) = scalar_normal_velocity_at_edge_full_level_rk%f(:,ie)

                call calc_hpe_vert_advection(scalar_template_1d_nlevp_a, &
                                        scalar_template_1d_nlev_a , &
                                        scalar_template_1d_nlev_b , &
                                        scalar_template_1d_nlev_c )

                ref_tend_normal_velocity_vadv_f(1:nlev,ie) = scalar_template_1d_nlev_c

        end do
!$omp end do nowait
!$omp end parallel
call t_stopf("dycore_time_intergration_5")
call data_check(dycoreVarEdgeFull%tend_normal_velocity_vadv%f, ref_tend_normal_velocity_vadv_f)
deallocate(ref_tend_normal_velocity_vadv_f)
#endif


#ifndef SEQ_GRIST
        call clock_end(clock_vau)
#endif

        if(iblock .eq. 0.and.write_verbose) print*,"finish vertical advection of normal wind"

!===========================================================================================
! diagnose dpim/dt (omegaMoist)
! Besides being a diagnostics:
! for HDC this diag is necessary for WRF physics; for NDC this diag produces similar results
! as we use www, g and rhom for diag, but not itentical even in logic;
! 2021-01: to use this for ndc's cu scheme make tropical rainfall looks better real-world 
!===========================================================================================

#ifndef SEQ_GRIST
        call clock_begin(clock_diagom)
#endif
        if(irk_step.eq.nrk)then ! diagnose d(mpressure)/dt for this dycore_step (final rk step)
           mesh%ne = mesh%ne_compute
! pi flux
! inline_opt_here
           !call calc_primal_normal_flux_at_edge(mesh, dycoreVarEdgeFull%scalar_normal_mass_flux_n%f, & ! dum
           !                                           dycoreVarEdgeFull%scalar_normal_mass_flux_n%f, & ! dum
           !                                           scalar_mpressure_at_pc_full_level_old         , & ! halo(1)->ne
           !                                           scalar_mpressure_at_edge_full_level_n         , &
           !                                           2, rk_substep, nlev)
#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
 allocate(ref_scalar_mpressure_at_edge_full_level_n, source=scalar_mpressure_at_edge_full_level_n)
#endif
call t_startf("dycore_time_intergration_target_6")
!$omp target
!$omp parallel  private(ie,ilev,icell1,icell2)
!$omp do
          do ie = 1, mesh%ne
             iCell1      = mesh%edt_v(1, ie)
             iCell2      = mesh%edt_v(2, ie)
             do ilev = 1, nlev
                scalar_mpressure_at_edge_full_level_n(ilev, ie) = 0.5*(scalar_mpressure_at_pc_full_level_old(ilev,iCell1)+&
                                                                       scalar_mpressure_at_pc_full_level_old(ilev,iCell2))
             end do
          end do
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("dycore_time_intergration_target_6")

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
call t_startf("dycore_time_intergration_6")
!$omp parallel  private(ie,ilev,icell1,icell2)
!$omp do
        do ie = 1, mesh%ne
        iCell1      = mesh%edt_v(1, ie)
        iCell2      = mesh%edt_v(2, ie)
        do ilev = 1, nlev
                ref_scalar_mpressure_at_edge_full_level_n(ilev, ie) = 0.5*(scalar_mpressure_at_pc_full_level_old(ilev,iCell1)+&
                                                                scalar_mpressure_at_pc_full_level_old(ilev,iCell2))
        end do
        end do
!$omp end do nowait
!$omp end parallel
call t_stopf("dycore_time_intergration_6")

call data_check(scalar_mpressure_at_edge_full_level_n, ref_scalar_mpressure_at_edge_full_level_n)
deallocate(ref_scalar_mpressure_at_edge_full_level_n)
#endif         
! inline_opt_here

! pi*dpi*V flux (old-time)
           call t_startf("dycore_time_intergration_pi_dpi_V")
           !$omp target
           !$omp parallel workshare
           scalar_mpressure_at_edge_full_level_n(:,:) = scalar_mpressure_at_edge_full_level_n(:,:)*dycoreVarEdgeFull%scalar_normal_mass_flux_n%f(:,:)
           !$omp end parallel workshare
           !$omp end target
           call t_stopf("dycore_time_intergration_pi_dpi_V")
! inline_opt_here
           !call divergence_operator_2d(mesh,scalar_mpressure_at_edge_full_level_n, &
           !                                 tend_pidpiv_at_pc_full_level_hori    , &
           !                                 nlev)
#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
allocate(ref_tend_pidpiv_at_pc_full_level_hori, source=tend_pidpiv_at_pc_full_level_hori)
#endif
call t_startf("dycore_time_intergration_target_7")
!$omp target
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev) 
!$omp do
           do iv = 1, mesh%nv
              div_sum(:) = 0._r8
              do inb = 1, mesh%vtx_nnb(iv)
                 index_edge  = mesh%vtx_ed(inb, iv)
                 do ilev = 1, nlev
                    div_sum(ilev)  = div_sum(ilev)+scalar_mpressure_at_edge_full_level_n(ilev,index_edge)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(index_edge)
                 end do
              end do
              tend_pidpiv_at_pc_full_level_hori(:,iv) = div_sum(:)/((rearth**2)*mesh%plg_areag(iv))
           end do
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("dycore_time_intergration_target_7")

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
call t_startf("dycore_time_intergration_7")
!$omp parallel  private(iv,div_sum,inb,index_edge,ilev) 
!$omp do
        do iv = 1, mesh%nv
                div_sum(:) = 0._r8
                do inb = 1, mesh%vtx_nnb(iv)
                index_edge  = mesh%vtx_ed(inb, iv)
                do ilev = 1, nlev
                div_sum(ilev)  = div_sum(ilev)+scalar_mpressure_at_edge_full_level_n(ilev,index_edge)*mesh%plg_nr(inb,iv)*rearth*mesh%edp_leng(index_edge)
                end do
                end do
                ref_tend_pidpiv_at_pc_full_level_hori(:,iv) = div_sum(:)/((rearth**2)*mesh%plg_areag(iv))
        end do
!$omp end do nowait
!$omp end parallel
call t_stopf("dycore_time_intergration_7")
call data_check(tend_pidpiv_at_pc_full_level_hori, ref_tend_pidpiv_at_pc_full_level_hori)
deallocate(ref_tend_pidpiv_at_pc_full_level_hori)
#endif     

! inline_opt_here
#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
           allocate(ref_scalar_omega_n_f, source=dycoreVarCellFull%scalar_omega_n%f)
#endif
           mesh%nv = mesh%nv_compute
           call t_startf("dycore_time_intergration_7+")
           !$omp target
           !$omp parallel private(iv, ilev)
           !$omp do 
           do iv = 1, mesh%nv
                do ilev = 1, nlev
                   dycoreVarCellFull%scalar_omega_n%f(ilev,iv) = (dycoreVarCellFull%scalar_mpressure_n%f(ilev,iv)-scalar_mpressure_at_pc_full_level_old(ilev,iv))/dtime+&
                                                                  dycoreVarCellFull%scalar_eta_mass_flux_n%f(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)+&
                                                              ((tend_pidpiv_at_pc_full_level_hori(ilev,iv)-dycoreVarCellFull%tend_mass_hori%f(ilev,iv)*&
                                                               scalar_mpressure_at_pc_full_level_old(ilev,iv))/scalar_delhp_at_pc_full_level_rk%f(ilev,iv))
                   if(eta_full_b(ilev).eq.0)then
                      dycoreVarCellFull%scalar_omega_n%f(ilev,iv) = dycoreVarCellFull%scalar_eta_mass_flux_n%f(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)
                   end if
                end do
           end do
           !$omp end do nowait
           !$omp end parallel
           !$omp end target
#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
           do iv = 1, mesh%nv
              do ilev = 1, nlev
                 ref_scalar_omega_n_f(ilev,iv) = (dycoreVarCellFull%scalar_mpressure_n%f(ilev,iv)-scalar_mpressure_at_pc_full_level_old(ilev,iv))/dtime+&
                                                                dycoreVarCellFull%scalar_eta_mass_flux_n%f(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)+&
                                                            ((tend_pidpiv_at_pc_full_level_hori(ilev,iv)-dycoreVarCellFull%tend_mass_hori%f(ilev,iv)*&
                                                             scalar_mpressure_at_pc_full_level_old(ilev,iv))/scalar_delhp_at_pc_full_level_rk%f(ilev,iv))
                 if(eta_full_b(ilev).eq.0)then
                        ref_scalar_omega_n_f(ilev,iv) = dycoreVarCellFull%scalar_eta_mass_flux_n%f(ilev,iv)/tracerVarCellFull%scalar_mif_n%f(ilev,iv)
                 end if
              end do
           end do
           call data_check(dycoreVarCellFull%scalar_omega_n%f, ref_scalar_omega_n_f)
           deallocate(ref_scalar_omega_n_f)
#endif
           call t_stopf("dycore_time_intergration_7+")
        mesh%nv = mesh%nv_halo(1)
        end if
#ifndef SEQ_GRIST
        call clock_end(clock_diagom)
#endif
!
! Until this line, HYDROSTATIC and NONHYDROSTATIC should be the same
!
#ifndef SEQ_GRIST
        call clock_begin(clock_nhd) ! either nhd or hpe
#endif
        IF(nh_dynamics)then

           if(iblock .eq. 0.and.write_verbose) print*,"before nh dynamics"
           call t_startf("dycore_time_intergration_before_nh")
           !$omp target
           !$omp parallel workshare
           dycoreVarCellFace%scalar_delp_n%f(2:nlev,:)  = dycoreVarCellFull%scalar_pressure_rk%f(2:nlev,:)-dycoreVarCellFull%scalar_pressure_rk%f(1:nlev-1,:)
           dycoreVarCellFace%scalar_delp_n%f(1,:)       = dycoreVarCellFull%scalar_pressure_rk%f(1,:)     -scalar_pressure_at_pc_face_level_rk%f(1,:)
           dycoreVarCellFace%scalar_delp_n%f(nlev+1,:)  = scalar_pressure_at_pc_face_level_rk%f(nlev+1,:) -dycoreVarCellFull%scalar_pressure_rk%f(nlev,:)
           !$omp end parallel workshare
           !$omp end target
           call t_stopf("dycore_time_intergration_before_nh")
#ifdef GRIST_DEBUG
           call debug_data_2d(irk_step,1,mesh%e_index,mesh%ne,"o1o",dycoreVarEdgeFace%scalar_normal_mass_flux_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o2o",dycoreVarCellFace%tend_mass_hori%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o3o",dycoreVarCellFull%scalar_eta_mass_flux_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o4o",scalar_phi_at_pc_face_level_rk%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o5o",scalar_www_at_pc_face_level_rk%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o6o",dycoreVarCellFace%scalar_delp_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o7o",dycoreVarCellFace%scalar_delhp_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o8o",dycoreVarCellFull%scalar_pressure_rk%f)
           call debug_data_1d(irk_step,  mesh%v_index,mesh%nv,"o9o",dycoreVarSurface%scalar_geopotential_n%f)!1d
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o10o",dycoreVarCellFull%scalar_pressure_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o11o",dycoreVarCellFace%scalar_phi_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o12o",dycoreVarCellFull%scalar_delhp_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o13o",dycoreVarCellFace%scalar_www_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o14o",dycoreVarCellFull%scalar_mass_pt_n%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o15o",dycoreVarCellFull%scalar_delhp_np1%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o16o",scalar_delhp_at_pc_face_level_np1%f)
           call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"o17o",scalar_mass_pt_at_pc_full_level_np1%f)
#endif
           call grist_nh_dynamics_run(mesh,rk_substep, itimestep, irk_step, idstep, itstep, &
                                         dycoreVarSurface%tend_hpressure_cnty,             & ! time tendency of hps
                                         scalar_hpressure_at_pc_face_level_np1,         & !
                                         dycoreVarEdgeFace%scalar_normal_mass_flux_n,   & ! rk, follow comment, donot fooled by vars' name which is for simplicity
                                         dycoreVarCellFace%tend_mass_hori,              & ! rk
                                         dycoreVarCellFull%scalar_eta_mass_flux_n,      & ! rk
                                         dycoreVarCellFace%scalar_eta_mass_flux_n,      & ! rk
                                         scalar_phi_at_pc_face_level_rk,                & ! rk
                                         scalar_www_at_pc_face_level_rk,                & ! rk
                                         dycoreVarCellFace%scalar_delp_n,               & ! rk
                                         scalar_delhp_at_pc_face_level_rk,              & ! rk
                                         dycoreVarCellFace%scalar_delhp_n,              & ! n
                                         dycoreVarCellFull%scalar_pressure_rk,          & ! rk
                                         dycoreVarSurface%scalar_geopotential_n,        & ! n
                                         dycoreVarCellFull%scalar_pressure_n,           & ! n
                                         dycoreVarCellFace%scalar_phi_n,                & ! n
                                         dycoreVarCellFull%scalar_delhp_n,              & ! n
                                         dycoreVarCellFace%scalar_www_n,                & ! n
                                         dycoreVarCellFull%scalar_mass_pt_n,            & ! n
                                         dycoreVarCellFull%scalar_delhp_np1,            & ! np1
                                         scalar_delhp_at_pc_face_level_np1,             & ! np1
                                         scalar_mass_pt_at_pc_full_level_np1,           & ! np1
                                         tracerVarCellFace%scalar_mif_n,                & ! tracer's n
                                         tracerVarCellFull%scalar_mif_n,                & ! tracer's n
                                         tracerVarCellFull%scalar_tracer_mxrt_n,        & ! tracer's n
                                         scalar_phi_at_pc_face_level_np1,               & ! np1, out
                                         scalar_www_at_pc_face_level_np1,               & ! np1, out
                                         dycoreVarCellFull%scalar_temp_n,               & ! np1, out!
                                         dycoreVarCellFull%scalar_geopotential_n,       & ! np1, out!
                                         dycoreVarCellFace%scalar_geopotential_n,       & ! np1, out
                                         dycoreVarCellFull%scalar_pressure_np1,         & ! np1, out
                                         dycoreVarCellFace%scalar_pressure_n,           & ! np1, out!
                                         dycoreVarCellFull%scalar_delp_np1,             & ! np1, out
                                         dycoreVarCellFull%scalar_alpha_np1)

#ifdef GRIST_DEBUG
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x1x",scalar_phi_at_pc_face_level_np1%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x2x",scalar_www_at_pc_face_level_np1%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x3x",dycoreVarCellFull%scalar_temp_n%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv_halo(1),"x4x",dycoreVarCellFull%scalar_geopotential_n%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x5x",dycoreVarCellFace%scalar_geopotential_n%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv,"x6x",dycoreVarCellFull%scalar_pressure_np1%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv_halo(1),"x7x",dycoreVarCellFace%scalar_pressure_n%f)
          call debug_data_2d(irk_step,1,mesh%v_index,mesh%nv_halo(1),"x8x",dycoreVarCellFull%scalar_delp_np1%f)
#endif
          if(iblock .eq. 0.and.write_verbose) print*,"end nh dynamics"
!
! optionally to obtain a hydrostatic balance state after minsteps for every-intsteps, or based on dynamic control
!
           IF((restore_hydro.and.itimestep.gt.restore_hydro_minsteps.and.mod(itimestep,restore_hydro_intsteps).eq.0).or.& ! mannually control
               ndc_restore_flag.eq.1)then ! dynamic control
! verified: if activate this under ndc every timestep, it produces exactly identical results as use a HDC,
! except some diagnostic variables and www
              call hydrostatic_adjust
              if(iblock.eq.0.and.write_stepinfo) print*,"instant restored at step,", itimestep
           ENDIF

           !IF(adjphi.and.itimestep.lt.int(10800/model_timestep))then
           IF(adjphi.and.itimestep.lt.16)then
              call hydrostatic_adjust
              !if(iblock.eq.0.and.write_stepinfo) print*,"restore hydro during model spinup", itimestep
           ENDIF
        ELSE

           call hydrostatic_adjust
           if(iblock .eq. 0.and.write_verbose) PRINT*,"finish diagnosing geopotential in HDC"

        END IF
#ifndef SEQ_GRIST
        call clock_end(clock_nhd)
#endif

!=================================================================================
! use below subroutines with smallest modification to compute PGF for normal
! velocity; The below part should not change in the case of nhdc or hdc
!=================================================================================
#ifndef SEQ_GRIST
      call clock_begin(clock_pgf)
#endif
      if(iblock .eq. 0.and.write_verbose) print*,"finish geop, grad_delp, grad_pressure gradient"
!
! compute real PGF tendency
!
      select case(hor_pgf_flag)

      case(6)  ! ptb, default with tp1 profile

          call grist_dycore_ref_atmos_create_ptb
          call t_startf("dycore_time_intergration_tend_normal_velocity_gz_fill")
          !$omp target
          !$omp parallel workshare
          dycoreVarEdgeFull%tend_normal_velocity_gz%f(:,:)  = 0._r8
          !$omp end parallel workshare
          !$omp end target
          call t_stopf("dycore_time_intergration_tend_normal_velocity_gz_fill")
#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
allocate(ref_tend_normal_velocity_pgf_f, source=dycoreVarEdgeFull%tend_normal_velocity_pgf%f)
#endif
call t_startf("dycore_time_intergration_target_8")
!$omp target
!$omp parallel  private(ie,icell1,icell2,cr_num,tmp1,tmp2,tmp3,tmp4)    
!$omp do
          do ie = 1, mesh%ne
             icell1 = mesh%edt_v(1,ie)
             icell2 = mesh%edt_v(2,ie)
             cr_num = rearth*mesh%edt_leng(ie)

             tmp1   = scalar_grad_hpres_at_edge_full_level_bar(:,ie)*&
                      0.5_r8*(scalar_alphad_at_pc_full_level_ptb(:,icell1)+scalar_alphad_at_pc_full_level_ptb(:,icell2))

             tmp2   = 0.5_r8*(dycoreVarCellFull%scalar_alpha_np1%f(:,icell1)+dycoreVarCellFull%scalar_alpha_np1%f(:,icell2))*&
                             (scalar_pressure_at_pc_full_level_ptb(:,icell2)-scalar_pressure_at_pc_full_level_ptb(:,icell1))/cr_num

             tmp3   = (scalar_geop_at_pc_full_level_ptb(:,icell2)-scalar_geop_at_pc_full_level_ptb(:,icell1))/cr_num
             tmp4   = 0.5_r8*(scalar_delp_at_pc_full_level_ptb(:,icell1)/dycoreVarCellFull%scalar_delhp_np1%f(:,icell1)+&
                              scalar_delp_at_pc_full_level_ptb(:,icell2)/dycoreVarCellFull%scalar_delhp_np1%f(:,icell2))*&
                             (dycoreVarCellFull%scalar_geopotential_n%f(:,icell2)-dycoreVarCellFull%scalar_geopotential_n%f(:,icell1))/cr_num

             dycoreVarEdgeFull%tend_normal_velocity_pgf%f(:,ie) = -1._r8*(tmp1+tmp2+tmp3+tmp4)*tracerVarEdgeFull%scalar_mif_n%f(:,ie)
          end do
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("dycore_time_intergration_target_8")

#if defined(CHK_DYCORE_TIME_INTERGRATION) || defined(CHK_ALL)
call t_startf("dycore_time_intergration_8")
!$omp parallel private(ie,icell1,icell2,cr_num,tmp1,tmp2,tmp3,tmp4)    
!$omp do
          do ie = 1, mesh%ne
             icell1 = mesh%edt_v(1,ie)
             icell2 = mesh%edt_v(2,ie)
             cr_num = rearth*mesh%edt_leng(ie)

             tmp1   = scalar_grad_hpres_at_edge_full_level_bar(:,ie)*&
                      0.5_r8*(scalar_alphad_at_pc_full_level_ptb(:,icell1)+scalar_alphad_at_pc_full_level_ptb(:,icell2))

             tmp2   = 0.5_r8*(dycoreVarCellFull%scalar_alpha_np1%f(:,icell1)+dycoreVarCellFull%scalar_alpha_np1%f(:,icell2))*&
                             (scalar_pressure_at_pc_full_level_ptb(:,icell2)-scalar_pressure_at_pc_full_level_ptb(:,icell1))/cr_num

             tmp3   = (scalar_geop_at_pc_full_level_ptb(:,icell2)-scalar_geop_at_pc_full_level_ptb(:,icell1))/cr_num
             tmp4   = 0.5_r8*(scalar_delp_at_pc_full_level_ptb(:,icell1)/dycoreVarCellFull%scalar_delhp_np1%f(:,icell1)+&
                              scalar_delp_at_pc_full_level_ptb(:,icell2)/dycoreVarCellFull%scalar_delhp_np1%f(:,icell2))*&
                             (dycoreVarCellFull%scalar_geopotential_n%f(:,icell2)-dycoreVarCellFull%scalar_geopotential_n%f(:,icell1))/cr_num

             ref_tend_normal_velocity_pgf_f(:,ie) = -1._r8*(tmp1+tmp2+tmp3+tmp4)*tracerVarEdgeFull%scalar_mif_n%f(:,ie)
          end do
!$omp end do nowait
!$omp end parallel
call t_stopf("dycore_time_intergration_8")
call data_check(ref_tend_normal_velocity_pgf_f, dycoreVarEdgeFull%tend_normal_velocity_pgf%f)
deallocate(ref_tend_normal_velocity_pgf_f)
#endif 
      case default
          if(iblock .eq. 0) print*, "you must select a PGF term in hor_pgf_flag"
          stop
      end select

      if(iblock .eq. 0.and.write_verbose) print*,"finish PGF"
#ifndef SEQ_GRIST
      call clock_end(clock_pgf)

      call clock_begin(clock_fnlrk)
#endif
!
! normal velocity update
!
         call t_startf("dycore_time_integration_nvu_alloc")
         if (.not. allocated(scalar_normal_velocity_at_edge_full_level_np1%f)) &
            allocate(scalar_normal_velocity_at_edge_full_level_np1%f, mold=dycoreVarEdgeFull%scalar_normal_velocity_n%f)
         call t_stopf("dycore_time_integration_nvu_alloc")
         call t_startf("dycore_time_intergration_normal_velocity_update1")
        !  dycoreVarEdgeFull%scalar_normal_velocity_n%f(:,:) 
         !$omp target parallel do private(scalar_normal_velocity_n_c, tend_normal_velocity_vadv_c, tend_normal_velocity_pgf_c, tend_normal_velocity_gz_c, tend_normal_velocity_ke_c, tend_normal_velocity_nct_c, tmp1_c, tmp2_c) schedule(static, 4)
         do ie = 1, mesh%ne
                call omnicopy(tmp1_c(:), dycoreVarEdgeFull%scalar_normal_velocity_n%f(:,ie))
                call omnicopy(tend_normal_velocity_vadv_c(:), dycoreVarEdgeFull%tend_normal_velocity_vadv%f(:,ie))
                call omnicopy(tend_normal_velocity_pgf_c(:), dycoreVarEdgeFull%tend_normal_velocity_pgf%f(:,ie))
                call omnicopy(tend_normal_velocity_gz_c(:), dycoreVarEdgeFull%tend_normal_velocity_gz%f(:,ie))
                call omnicopy(tend_normal_velocity_ke_c(:), dycoreVarEdgeFull%tend_normal_velocity_ke%f(:,ie))
                call omnicopy(tend_normal_velocity_nct_c(:), dycoreVarEdgeFull%tend_normal_velocity_nct%f(:,ie))
                tmp2_c(:) = tend_normal_velocity_vadv_c + tend_normal_velocity_pgf_c + tend_normal_velocity_gz_c + tend_normal_velocity_ke_c + tend_normal_velocity_nct_c
                tmp1_c(:) = tmp1_c + dtime/rk_number * tmp2_c
                call omnicopy(scalar_normal_velocity_at_edge_full_level_np1%f(:,ie), tmp1_c)
!          scalar_normal_velocity_at_edge_full_level_np1%f(:,ie) = &
!          dycoreVarEdgeFull%scalar_normal_velocity_n%f(:,ie) +dtime/rk_number*&
!        ( dycoreVarEdgeFull%tend_normal_velocity_vadv%f(:,ie) +&
!          dycoreVarEdgeFull%tend_normal_velocity_pgf%f(:,ie)  +&
!          dycoreVarEdgeFull%tend_normal_velocity_gz%f(:,ie)   +&
!          dycoreVarEdgeFull%tend_normal_velocity_ke%f(:,ie)   +&
!          dycoreVarEdgeFull%tend_normal_velocity_nct%f(:,ie))
         end do
         !$omp end target parallel do
         call t_stopf("dycore_time_intergration_normal_velocity_update1")
         call t_startf("dycore_time_intergration_normal_velocity_update2")
         if(use_phys.and.ptend_wind_rk_on)then
             !$omp target
             !$omp parallel workshare
             scalar_normal_velocity_at_edge_full_level_np1%f(:,:) = scalar_normal_velocity_at_edge_full_level_np1%f(:,:)+&
                                                   dtime/rk_number*ptend_rk%tend_normal_velocity_at_edge_full_level%f(:,:)
             !$omp end parallel workshare
             !$omp end target
         end if
         call t_stopf("dycore_time_intergration_normal_velocity_update2")
!
! diffusion terms added here as phys, must be written in this way for bit regression
!
         call t_startf("dycore_time_intergration_ad_dycore_laplacian")
         if(ad_dycore_laplacian_2nd.and. .not.ad_dycore_laplacian_4th)then
             !$omp target
             !$omp parallel workshare
             scalar_normal_velocity_at_edge_full_level_np1%f(:,:) = scalar_normal_velocity_at_edge_full_level_np1%f(:,:)+&
                                                   dtime/rk_number*dycoreVarEdgeFull%tend_hwind_laplacian_2nd%f(:,:)
             !$omp end parallel workshare
             !$omp end target
         else if(ad_dycore_laplacian_4th .and. .not. ad_dycore_laplacian_2nd)then
             !$omp target
             !$omp parallel workshare
             scalar_normal_velocity_at_edge_full_level_np1%f(:,:) = scalar_normal_velocity_at_edge_full_level_np1%f(:,:)+&
                                                   dtime/rk_number*dycoreVarEdgeFull%tend_hwind_laplacian_4th%f(:,:)
             !$omp end parallel workshare
             !$omp end target
         else if(ad_dycore_laplacian_2nd.and.ad_dycore_laplacian_4th)then
             !$omp target
             !$omp parallel workshare
             scalar_normal_velocity_at_edge_full_level_np1%f(:,:) = scalar_normal_velocity_at_edge_full_level_np1%f(:,:)+&
                                                   dtime/rk_number*(dycoreVarEdgeFull%tend_hwind_laplacian_4th%f(:,:)+&
                                                                    dycoreVarEdgeFull%tend_hwind_laplacian_2nd%f(:,:))
             !$omp end parallel workshare
             !$omp end target
         end if
!
! add Lap4th for www equation in nh_dynamics
!
         if(nh_dynamics.and.ad_dycore_laplacian_4th.and.use_www_hyperDiffusion)then
            !$omp target
            !$omp parallel workshare
            scalar_www_at_pc_face_level_np1%f(:,:) = scalar_www_at_pc_face_level_np1%f(:,:)+dtime/rk_number*dycoreVarCellFace%tend_www_laplacian_4th%f(:,:)
            !$omp end parallel workshare
            !$omp end target
         end if

#if (!defined DCMIP21)
         if(nh_dynamics.and.ad_dycore_laplacian_2nd)then
             !$omp target
             !$omp parallel workshare
             scalar_www_at_pc_face_level_np1%f(:,:) = scalar_www_at_pc_face_level_np1%f(:,:)+&
                                                 dtime/rk_number*dycoreVarCellFace%tend_www_laplacian_2nd%f(:,:)
            !$omp end parallel workshare
            !$omp end target
         end if
#endif
!
! diff 6th is only for testing, not used for all cases
!
         if(ad_dycore_laplacian_6th)then
             !$omp target
             !$omp parallel workshare
             scalar_normal_velocity_at_edge_full_level_np1%f(:,:) = scalar_normal_velocity_at_edge_full_level_np1%f(:,:)+&
                                                   dtime/rk_number*dycoreVarEdgeFull%tend_hwind_laplacian_6th%f(:,:)
             !$omp end parallel workshare
             !$omp end target
         end if
         call t_stopf("dycore_time_intergration_ad_dycore_laplacian")
!
! RK-sub-step update
!
#ifdef GRIST_DEBUG
       if(irk_step .eq. 1 .and. itimestep .eq. int(nsteps)) then
       call debug_data_1d(irk_step,      mesh%v_index,mesh%nv,"shaps_n" ,scalar_hpressure_at_pc_surface_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%v_index,mesh%nv,"smpapf_n",scalar_mass_pt_at_pc_full_level_np1%f)
       call debug_data_2d(irk_step,nlevp,mesh%v_index,mesh%nv,"swapf_n" ,scalar_www_at_pc_face_level_np1%f)
       call debug_data_2d(irk_step,nlevp,mesh%v_index,mesh%nv,"sdpapf_n",scalar_phi_at_pc_face_level_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%e_index,mesh%ne,"snvaef_n",scalar_normal_velocity_at_edge_full_level_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%v_index,mesh%nv,"sdhapf_n",dycoreVarCellFull%scalar_delhp_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%v_index,mesh%nv,"spapf_n" ,dycoreVarCellFull%scalar_pressure_np1%f)      
       call debug_data_2d(irk_step,nlevp,mesh%v_index,mesh%nv,"diag_1"  ,scalar_delhp_at_pc_face_level_np1%f)
       call debug_data_2d(irk_step,nlev, mesh%v_index,mesh%nv,"diag_2"  ,dycoreVarCellFull%scalar_potential_temp_n%f) !step 2 error
       call debug_data_2d(irk_step,nlevp,mesh%v_index,mesh%nv,"diag_3"  ,scalar_hpressure_at_pc_face_level_np1%f)
       end if
#endif

#ifdef LAM_DOMAIN
       call update_assignValueArea_2d(mesh,nlev,"normal_velocity",scalar_normal_velocity_at_edge_full_level_np1)
       call update_assignValueArea_2d(mesh,nlev,"mass_pt"        ,scalar_mass_pt_at_pc_full_level_np1)
#endif
!Exchange data
#ifndef SEQ_GRIST
       call clock_begin(clock_mainexch)
       call exchange_data_2d_add(mesh,field_head_2d,scalar_mass_pt_at_pc_full_level_np1)
       call exchange_data_2d_add(mesh,field_head_2d,scalar_normal_velocity_at_edge_full_level_np1)
#endif

       if(nh_dynamics)then
#ifdef LAM_DOMAIN
          call update_assignValueArea_2d(mesh,nlevp,"wwwFace", scalar_www_at_pc_face_level_np1)
          call update_assignValueArea_2d(mesh,nlevp,"phiFace", scalar_phi_at_pc_face_level_np1)
#endif
#ifndef SEQ_GRIST
          call exchange_data_2d_add(mesh,field_head_2d,scalar_www_at_pc_face_level_np1)
          call exchange_data_2d_add(mesh,field_head_2d,scalar_phi_at_pc_face_level_np1)
#endif
       end if

#ifndef SEQ_GRIST
       call exchange_data_2d(mesh%local_block,field_head_2d)
       call clock_end(clock_mainexch)
#endif
       call t_startf("dycore_time_intergration_fnlrk_op")
       !$omp target
       !$omp parallel workshare
       dycoreVarCellFull%scalar_potential_temp_n%f(:,:) = scalar_mass_pt_at_pc_full_level_np1%f(:,:)/& ! Overwritten each RK step
                                                     dycoreVarCellFull%scalar_delhp_np1%f(:,:)
       !$omp end parallel workshare
       !$omp end target
! prog
      if(nh_dynamics)then
         !$omp target
         !$omp parallel workshare
         scalar_www_at_pc_face_level_rk%f(:,:)            = scalar_www_at_pc_face_level_np1%f(:,:)
         scalar_phi_at_pc_face_level_rk%f(:,:)            = scalar_phi_at_pc_face_level_np1%f(:,:)
         dycoreVarCellFull%scalar_pressure_rk%f(:,:)      = dycoreVarCellFull%scalar_pressure_np1%f(:,:)
         scalar_pressure_at_pc_face_level_rk%f(:,:)       = dycoreVarCellFace%scalar_pressure_n%f(:,:)
         !$omp end parallel workshare
         !$omp end target
         scalar_www_at_pc_face_level_rk%pos            = scalar_www_at_pc_face_level_np1%pos
         scalar_phi_at_pc_face_level_rk%pos            = scalar_phi_at_pc_face_level_np1%pos
         dycoreVarCellFull%scalar_pressure_rk%pos      = dycoreVarCellFull%scalar_pressure_np1%pos
         scalar_pressure_at_pc_face_level_rk%pos       = dycoreVarCellFace%scalar_pressure_n%pos
      end if
      !$omp target
      !$omp parallel workshare
      scalar_normal_velocity_at_edge_full_level_rk%f(:,:) = scalar_normal_velocity_at_edge_full_level_np1%f(:,:)
      scalar_delhp_at_pc_full_level_rk%f(:,:)             = dycoreVarCellFull%scalar_delhp_np1%f(:,:)
! diag
      scalar_delhp_at_pc_face_level_rk%f(:,:)             = scalar_delhp_at_pc_face_level_np1%f(:,:)
      !$omp end parallel workshare
      !$omp end target
      scalar_normal_velocity_at_edge_full_level_rk%pos = scalar_normal_velocity_at_edge_full_level_np1%pos
      scalar_delhp_at_pc_full_level_rk%pos             = dycoreVarCellFull%scalar_delhp_np1%pos
      scalar_delhp_at_pc_face_level_rk%pos             = scalar_delhp_at_pc_face_level_np1%pos
      call t_stopf("dycore_time_intergration_fnlrk_op")
      if(iblock .eq. 0 .and. write_stepinfo) print*,"---------- rkfb ",irk_step,"----------"
#ifndef SEQ_GRIST
      call clock_end(clock_fnlrk)
#endif

   END DO
#ifndef SEQ_GRIST
   call clock_end(clock_rkl)
#endif

!----------------------------------------------
! final update for the next dycore step
!----------------------------------------------

      if(nh_dynamics)then
         dycoreVarCellFace%scalar_www_n      = scalar_www_at_pc_face_level_np1
         dycoreVarCellFace%scalar_phi_n      = scalar_phi_at_pc_face_level_np1
         dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_pressure_np1
      end if
      dycoreVarEdgeFull%scalar_normal_velocity_n  = scalar_normal_velocity_at_edge_full_level_np1
      dycoreVarCellFull%scalar_mass_pt_n          = scalar_mass_pt_at_pc_full_level_np1
      dycoreVarSurface%scalar_hpressure_n         = scalar_hpressure_at_pc_surface_np1
      dycoreVarCellFull%scalar_delhp_n            = dycoreVarCellFull%scalar_delhp_np1    ! needed by physics and dtp
      dycoreVarCellFace%scalar_delhp_n            = scalar_delhp_at_pc_face_level_np1

!
! this is only for diagnose PS, not put to dycore_diag because of historry-h1 does not dycore_diag at each step
      dycoreVarSurface%scalar_pressure_n%f = dycoreVarCellFace%scalar_pressure_n%f(nlevp,:)
!
! substract diabatic physical influence, used with F2, so we should substract F2
! here; for idealized test, f2 has the same tendency, while for full-physics
! test, rk has its own total tendency, only the portion that f2 covers should be
! substracted
!
call t_startf("substract_diabatic_physical_influence_target")      
      if(ptendSubDiabPhys)then
!$omp target
!$omp parallel workshare
          dycoreVarCellFull%scalar_mass_pt_n%f(:,:) =  dycoreVarCellFull%scalar_mass_pt_n%f(:,:)-&
                                                 !dtime*ptend_rk%tend_mass_pt_at_pc_full_level%f
                                                 dtime*ptend_f2%tend_mass_pt_at_pc_full_level%f(:,:)
!$omp end parallel workshare
!$omp end target 
#ifndef SEQ_GRIST
          call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
          call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
!$omp target
!$omp parallel workshare
          dycoreVarCellFull%scalar_potential_temp_n%f(:,:) = dycoreVarCellFull%scalar_mass_pt_n%f(:,:)/& ! overwritten each RK step
                                                       dycoreVarCellFull%scalar_delhp_n%f(:,:)
!$omp end parallel workshare
!$omp end target     
        end if
      call t_stopf("substract_diabatic_physical_influence_target")

!================================================
! pt now only affected by adv and diffusion, get 
! this tend and accumulate until this itimestep
! reachs step_cu and all dycore+tracer steps done
! this tend is only for cumulus tiedtke
!================================================

call t_startf("dycore_time_int_last")
#ifdef AMIPW_PHYSICS
      IF(present(idstep))then

      dycoreVarCellFull%tend_pt_n%f = dycoreVarCellFull%tend_pt_n%f+&
                                    (dycoreVarCellFull%scalar_potential_temp_n%f/(one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,:))-&
                                     scalar_pt_at_pc_full_level_old)/dtime
      if(idstep.eq.dstep_in_tstep.and.itstep.eq.tstep_in_mstep.and.mod(itimestep,step_cu).eq.0)then
         dycoreVarCellFull%tend_pt_n%f = dycoreVarCellFull%tend_pt_n%f/(dstep_in_tstep*tstep_in_mstep*step_cu)
      end if
      END IF
#endif

!=====================================================================================
!                          Coupling of F1 physics forcing
! after each dycore step, _n is the updated state, on which model physics reside
! Recent tests suggest that F1 physics for dycore is less useful,
! For HDC, just use ptend_f2_sudden;
! For NDC, use ptend_rk for wind, heat (slow-physics), use ptend_f2_sudden
! for tracer update, and some OS heat and wind (fast-physics)
! These codes are still here for reference purpose, but we added a compipling option.
! So as in tracer transport' F1
! 2021-01-05: Based on many full-physics tests, we can basically ignore F1-tend
! can be replaced by ptend_rk, ptend_f2 and f-s combination
!=====================================================================================

#ifdef USE_PTENDF1
      if(use_phys.and.ptend_dycore_heat_f1_on)then
         dycoreVarCellFull%scalar_mass_pt_n%f           = dycoreVarCellFull%scalar_mass_pt_n%f+&
                                                         ptend_f1%tend_mass_pt_at_pc_full_level%f*dtime
#ifndef SEQ_GRIST
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
         call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
         dycoreVarCellFull%scalar_potential_temp_n%f = dycoreVarCellFull%scalar_mass_pt_n%f/dycoreVarCellFull%scalar_delhp_n%f
      end if

      if(use_phys.and.ptend_dycore_wind_f1_on)then
         dycoreVarEdgeFull%scalar_normal_velocity_n%f = dycoreVarEdgeFull%scalar_normal_velocity_n%f+&
                                                         ptend_f1%tend_normal_velocity_at_edge_full_level%f*dtime
#ifndef SEQ_GRIST
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
         call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
      end if

      if(use_phys.and.ptend_dycore_f1_on)then
         dycoreVarCellFull%scalar_mass_pt_n%f           = dycoreVarCellFull%scalar_mass_pt_n%f+&
                                                         ptend_f1%tend_mass_pt_at_pc_full_level%f*dtime
         dycoreVarEdgeFull%scalar_normal_velocity_n%f = dycoreVarEdgeFull%scalar_normal_velocity_n%f+&
                                                         ptend_f1%tend_normal_velocity_at_edge_full_level%f*dtime
#ifndef SEQ_GRIST
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarCellFull%scalar_mass_pt_n)
         call exchange_data_2d_add(mesh,field_head_2d,dycoreVarEdgeFull%scalar_normal_velocity_n)
         call exchange_data_2d(mesh%local_block,field_head_2d)
#endif
         dycoreVarCellFull%scalar_potential_temp_n%f = dycoreVarCellFull%scalar_mass_pt_n%f/dycoreVarCellFull%scalar_delhp_n%f
      end if
#endif
call t_stopf("dycore_time_int_last")

      return

contains

    subroutine hydrostatic_adjust
!
! the modifiction of this part for mhdc may affect bit reproduce of dry hdc,
! although physically identical
!       

!#define CHK_HYDRO_ADJ
#if defined(CHK_HYDRO_ADJ) || defined(CHK_ALL)
#include <data_check.inc>
        real (r8), allocatable :: facef(:,:), fullf(:,:), np1f(:,:), tempf(:,:)
        integer :: rank
        rank = mpi_rank()
        allocate(facef, source=dycoreVarCellFace%scalar_geopotential_n%f)
        allocate(fullf, source=dycoreVarCellFull%scalar_geopotential_n%f)
        allocate(np1f, source=scalar_phi_at_pc_face_level_np1%f)
        allocate(tempf, source=dycoreVarCellFull%scalar_temp_n%f)
        if (rank == 3) print *, "fulltempf(5, 7420)=", dycoreVarCellFull%scalar_temp_n%f(5,7420), tempf(5,7420)
#endif
        !$omp target
        !$omp parallel workshare
        dycoreVarCellFull%scalar_pressure_np1%f(:,:) = dycoreVarCellFull%scalar_mpressure_n%f(:,:)
        dycoreVarCellFace%scalar_pressure_n%f(:,:)   = dycoreVarCellFace%scalar_mpressure_n%f(:,:)
        !$omp end parallel workshare
        !$omp end target
        call t_startf("hydrostatic_adjust")
!
! diagnose geopotential and evaluate its gradient
!
!$omp target
!$omp parallel  private(iv,scalar_template_1d_nlev_a,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_b,scalar_template_a,scalar_template_1d_nlevp_b,scalar_template_1d_nlev_c)
!$omp do 
        do iv = 1, mesh%nv_halo(1)
          !
          ! recover temperature, and store
          !
           !scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_potential_temp_np1%f(:,iv)*&
           scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_potential_temp_n%f(:,iv)*&
                                      ((dycoreVarCellFull%scalar_pressure_np1%f(:,iv)/p00)**(rdry/cp))/&
                                       (one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv))

           dycoreVarCellFull%scalar_temp_n%f(:,iv) = scalar_template_1d_nlev_a ! output, overwritten each RK step
           scalar_template_1d_nlev_a  = scalar_template_1d_nlev_a*(one+(rvap-rdry)/rdry*(tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv)/&
                                                (one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv)))) ! Tv
          !
          ! new pressure face
          !
           !scalar_template_1d_nlevp_a = scalar_hpressure_at_pc_face_level_n%f(:,iv)
           scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_pressure_n%f(:,iv)
          !
          ! new delp full, not in effect actually
          !
           !scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delhp_np1%f(:,iv)
           scalar_template_1d_nlev_b  = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,iv)-dycoreVarCellFace%scalar_pressure_n%f(1:nlev,iv)
          !
          ! geopotential at surface
          !
           scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)
          ! diagnose geop at face based on pressure, which is related to
          ! the definition of model layer
           call calc_hpe_hydro(scalar_template_1d_nlev_a  ,& ! T or Tv
                               scalar_template_1d_nlevp_a ,& ! hp or p
                               scalar_template_1d_nlev_b  ,& ! delhp or delp
                               scalar_template_a          ,& ! phis
                               scalar_template_1d_nlevp_b ,& ! phi-face
                               scalar_template_1d_nlev_c )   ! phi-full
           dycoreVarCellFace%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlevp_b
           dycoreVarCellFull%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlev_c  ! output
           scalar_phi_at_pc_face_level_np1%f(:,iv)         = dycoreVarCellFace%scalar_geopotential_n%f(:,iv)

        end do
!$omp end do nowait
!$omp end parallel 
!$omp end target
#if defined(CHK_HYDRO_ADJ) || defined(CHK_ALL)
        ! if (rank == 3) print *, "fulltempf(5, 7420)=", dycoreVarCellFull%scalar_temp_n%f(5,7420), tempf(5,7420)
        ! print *, "nvh=", mesh%nv_halo(1)
!$omp parallel  private(iv,scalar_template_1d_nlev_a,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_b,scalar_template_a,scalar_template_1d_nlevp_b,scalar_template_1d_nlev_c)
!$omp do
        do iv = 1, mesh%nv_halo(1)
        !
        ! recover temperature, and store
        !
         !scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_potential_temp_np1%f(:,iv)*&
           scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_potential_temp_n%f(:,iv)*&
                                      ((dycoreVarCellFull%scalar_pressure_np1%f(:,iv)/p00)**(rdry/cp))/&
                                       (one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv))  
  
           tempf(:,iv) = scalar_template_1d_nlev_a ! output, overwritten each RK step
           scalar_template_1d_nlev_a  = scalar_template_1d_nlev_a*(one+(rvap-rdry)/rdry*(tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv)/&
                                                (one+tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv)))) ! Tv
          !
          ! new pressure face
          !
           !scalar_template_1d_nlevp_a = scalar_hpressure_at_pc_face_level_n%f(:,iv)
           scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_pressure_n%f(:,iv)
          !
          ! new delp full, not in effect actually
          !
           !scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delhp_np1%f(:,iv)
           scalar_template_1d_nlev_b  = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,iv)-dycoreVarCellFace%scalar_pressure_n%f(1:nlev,iv)
          !
          ! geopotential at surface
          !
           scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv)
          ! diagnose geop at face based on pressure, which is related to
          ! the definition of model layer
           call calc_hpe_hydro(scalar_template_1d_nlev_a  ,& ! T or Tv
                               scalar_template_1d_nlevp_a ,& ! hp or p
                               scalar_template_1d_nlev_b  ,& ! delhp or delp
                               scalar_template_a          ,& ! phis
                               scalar_template_1d_nlevp_b ,& ! phi-face
                               scalar_template_1d_nlev_c )   ! phi-full
           if (iv == 7420 .and. rank == 3) print *, "tempf(5, 7420)=", tempf(5,7420)
           facef(:,iv) = scalar_template_1d_nlevp_b
           fullf(:,iv) = scalar_template_1d_nlev_c  ! output
           np1f(:,iv)         = facef(:,iv)

        end do
!$omp end do nowait
!$omp end parallel 
        if (rank == 3) print *, "fulltempf(5, 7420)=", dycoreVarCellFull%scalar_temp_n%f(5,7420), tempf(5,7420)
        call data_check(dycoreVarCellFace%scalar_geopotential_n%f, facef)
        call data_check(dycoreVarCellFull%scalar_geopotential_n%f, fullf)
        call data_check(dycoreVarCellFull%scalar_temp_n%f, tempf)
        call data_check(scalar_phi_at_pc_face_level_np1%f, np1f)
        deallocate(facef, fullf, np1f, tempf)
#endif
          call t_stopf("hydrostatic_adjust")
          call t_startf("hydrostatic_adjust_tail")
          ! since in hdc, delp is only used for pgf, and we use (6), so delp should be delhp, to cancel the last term in pgf6, not so !
          ! dycoreVarCellFull%scalar_delp_np1%f   = dycoreVarCellFull%scalar_delhp_np1%f
          ! new
          !$omp target
          !$omp parallel workshare
          dycoreVarCellFull%scalar_delp_np1%f(:,:)   = dycoreVarCellFace%scalar_pressure_n%f(2:nlevp,:)-dycoreVarCellFace%scalar_pressure_n%f(1:nlev,:) 
          ! old dry-hydrostatic model only uses this
          !dycoreVarCellFull%scalar_alpha_np1%f  = rdry*dycoreVarCellFull%scalar_temp_n%f/scalar_hpressure_at_pc_full_level_np1%f
          ! real layer-averaged
          dycoreVarCellFull%scalar_alpha_np1%f(:,:)  = (dycoreVarCellFace%scalar_geopotential_n%f(1:nlev,:)-dycoreVarCellFace%scalar_geopotential_n%f(2:nlevp,:))/&
                                                  dycoreVarCellFull%scalar_delhp_np1%f(:,:)
          !$omp end parallel workshare
          !$omp end target
          call t_stopf("hydrostatic_adjust_tail")
          return
    end subroutine hydrostatic_adjust

   end subroutine dycore_time_integration_run

   subroutine dycore_time_integration_init(mesh)
! dycore-nh
     use grist_nh_driver_module,           only: grist_nh_dynamics_init
     use grist_nh_explicit_tend_module_2d, only: grist_nh_et_init
! io
     type(global_domain),  intent(inout) :: mesh
! local
     integer(i4)      :: iv

     call grist_nh_dynamics_init(mesh)
     call grist_nh_et_init(mesh)

        if(.not.allocated(scalar_template_1d_nlev_a))  allocate(scalar_template_1d_nlev_a(nlev))
        if(.not.allocated(scalar_template_1d_nlev_b))  allocate(scalar_template_1d_nlev_b(nlev))
        if(.not.allocated(scalar_template_1d_nlev_c))  allocate(scalar_template_1d_nlev_c(nlev))
        if(.not.allocated(scalar_template_1d_nlev_d))  allocate(scalar_template_1d_nlev_d(nlev))
        if(.not.allocated(scalar_template_1d_nlevp_a)) allocate(scalar_template_1d_nlevp_a(nlevp))
        if(.not.allocated(scalar_template_1d_nlevp_b)) allocate(scalar_template_1d_nlevp_b(nlevp))
        if(.not.allocated(scalar_template_1d_nlevp_c)) allocate(scalar_template_1d_nlevp_c(nlevp))
        if(.not.allocated(scalar_template_2d_ne_b%f))  allocate(scalar_template_2d_ne_b%f(nlev,mesh%ne))
!
! initially compute hpressure (full&face), delhp (full&face), based on surface hpressure
!
        call time_integration_renew_mass_state(mesh, dycoreVarSurface%scalar_hpressure_n    ,&
                                                     dycoreVarCellFace%scalar_hpressure_n ,&
                                                     dycoreVarCellFull%scalar_delhp_n     ,&
                                                     dycoreVarCellFull%scalar_hpressure_n ,&
                                                     dycoreVarCellFace%scalar_delhp_n    )

        if(trim(working_mode).eq.'dycore')then
           dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_hpressure_n ! init as hpressure
           dycoreVarCellFace%scalar_pressure_n = dycoreVarCellFace%scalar_hpressure_n ! init as hpressure
           scalar_pressure_at_pc_face_level_rk = dycoreVarCellFace%scalar_hpressure_n ! init as hpressure
        end if
! above for dycore mode

        if(trim(working_mode).ne.'dycore')then
           ! already defined in dtp
           !dycoreVarCellFull%scalar_pressure_n = dycoreVarCellFull%scalar_pressure_n
           !dycoreVarCellFace%scalar_pressure_n = dycoreVarCellFace%scalar_pressure_n
           scalar_pressure_at_pc_face_level_rk= dycoreVarCellFace%scalar_pressure_n
        end if

        do iv = 1, mesh%nv
            dycoreVarCellFull%scalar_potential_temp_n%f(:,iv)  =  &
                                                          dycoreVarCellFull%scalar_temp_n%f(:,iv)/&
                                                  ((dycoreVarCellFull%scalar_pressure_n%f(:,iv)/p00)**(rdry/cp))*&
                                                  (one+ptfactor*tracerVarCellFull%scalar_tracer_mxrt_n%f(1,:,iv))

            dycoreVarCellFull%scalar_mass_pt_n%f(:,iv)  =  dycoreVarCellFull%scalar_potential_temp_n%f(:,iv)*&
                                                          dycoreVarCellFull%scalar_delhp_n%f(:,iv)
        end do

        IF(trim(gcm_testcase).ne.'SCHAR'     .and.&
           trim(gcm_testcase).ne.'DCMIP2-1'  .and.&
           trim(gcm_testcase).ne.'real-ERAIM'.and.&
           trim(gcm_testcase).ne.'real-ERAIP'.and.&
           trim(gcm_testcase).ne.'real-GFS'  .and.&
           trim(gcm_testcase).ne.'real-WRFDA'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-BW'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-TC'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-SC-A'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-SC-B'.and.&
           trim(gcm_testcase).ne.'DCMIP2016-SC1' .and.&
           trim(gcm_testcase).ne.'DCMIP2016-SCXX' .and.&
           trim(gcm_testcase).ne.'DCMIP2016-SC')then
        
        do iv = 1, mesh%nv

           scalar_template_1d_nlev_a  = dycoreVarCellFull%scalar_temp_n%f(:,iv)      ! temperature
           scalar_template_1d_nlevp_a = dycoreVarCellFace%scalar_hpressure_n%f(:,iv) ! hpressure face
           scalar_template_1d_nlev_b  = dycoreVarCellFull%scalar_delhp_n%f(:,iv)     ! delhp full
           scalar_template_a          = dycoreVarSurface%scalar_geopotential_n%f(iv) ! geopotential at surface

           call calc_hpe_hydro(scalar_template_1d_nlev_a  ,&
                               scalar_template_1d_nlevp_a ,&
                               scalar_template_1d_nlev_b  ,&
                               scalar_template_a          ,&
                               scalar_template_1d_nlevp_b ,&
                               scalar_template_1d_nlev_c )

           dycoreVarCellFace%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlevp_b
           dycoreVarCellFull%scalar_geopotential_n%f(:,iv) = scalar_template_1d_nlev_c ! output
           dycoreVarCellFace%scalar_phi_n%f(:,iv)          = dycoreVarCellFace%scalar_geopotential_n%f(:,iv)
        end do
        dycoreVarCellFace%scalar_www_n%f   = 0._r8
        dycoreVarCellFace%scalar_phi_n%f(nlev+1,:)  = dycoreVarSurface%scalar_geopotential_n%f ! surface is given in ini condition
        END IF

        mesh%ne = mesh%ne_compute
        call grist_dycore_ref_atmos_init(mesh)
        mesh%ne = mesh%ne_full

!
! init mpressure
!
        call calc_hpe_get_full_mass(nlev, mesh%nv_halo(1), working_mode          , & ! in
                                          dycoreVarCellFace%scalar_hpressure_n%f , & ! in
                                          dycoreVarCellFull%scalar_hpressure_n%f , & ! in
                                          dycoreVarCellFull%scalar_delhp_n%f     , & ! in
                                          tracerVarCellFull%scalar_mif_n%f       , & ! in
                                          dycoreVarCellFull%scalar_mpressure_n%f , & ! out
                                          dycoreVarCellFace%scalar_mpressure_n%f)    ! out
! init array
      scalar_www_at_pc_face_level_np1        = dycoreVarCellFace%scalar_www_n 
      scalar_phi_at_pc_face_level_np1        = dycoreVarCellFace%scalar_www_n
      scalar_delhp_at_pc_face_level_np1      = dycoreVarCellFace%scalar_www_n
      scalar_hpressure_at_pc_face_level_np1  = dycoreVarCellFace%scalar_www_n
      scalar_hpressure_at_pc_full_level_np1  = dycoreVarCellFull%scalar_hpressure_n
#ifndef SEQ_GRIST
      clock_mas      = clock_id('mas')
      clock_nct      = clock_id('nct') 
      clock_nhd      = clock_id('nhd') 
      clock_ptm      = clock_id('ptm')
      clock_pgf      = clock_id('pgf')
      clock_vau      = clock_id('vau')
      clock_ket      = clock_id('ket')
      clock_diagom   = clock_id('diagom')
      clock_mainexch = clock_id('mainexch')
      clock_rkl      = clock_id('rkl')
      clock_fnlrk    = clock_id('fnlrk')
#endif


      return
   end subroutine dycore_time_integration_init

   subroutine dycore_time_integration_final

     use grist_nh_driver_module,             only: grist_nh_dynamics_final
     use grist_nh_explicit_tend_module_2d,   only: grist_nh_et_final

      if(allocated(scalar_template_1d_nlev_a))  deallocate(scalar_template_1d_nlev_a)
      if(allocated(scalar_template_1d_nlev_b))  deallocate(scalar_template_1d_nlev_b)
      if(allocated(scalar_template_1d_nlev_c))  deallocate(scalar_template_1d_nlev_c)
      if(allocated(scalar_template_1d_nlev_d))  deallocate(scalar_template_1d_nlev_d)
      if(allocated(scalar_template_1d_nlevp_a)) deallocate(scalar_template_1d_nlevp_a)
      if(allocated(scalar_template_1d_nlevp_b)) deallocate(scalar_template_1d_nlevp_b)
      if(allocated(scalar_template_1d_nlevp_c)) deallocate(scalar_template_1d_nlevp_c)
      if(allocated(scalar_template_2d_ne_b%f))  deallocate(scalar_template_2d_ne_b%f)

      call grist_nh_dynamics_final
      call grist_nh_et_final
      call grist_dycore_ref_atmos_final

   end subroutine dycore_time_integration_final

!--------------------
!
! PRIVATE BELOW
!
!--------------------

   subroutine time_integration_renew_mass_state(mesh, scalar_hpressure_at_pc_surface    ,&
                                                      scalar_hpressure_at_pc_face_level ,&
                                                      scalar_delhp_at_pc_full_level     ,&
                                                      scalar_hpressure_at_pc_full_level ,&
                                                      scalar_delhp_at_pc_face_level     )
! io
     use omp_lib
     type(global_domain),   intent(in)    :: mesh
     type(scalar_1d_field), intent(in)    :: scalar_hpressure_at_pc_surface
     type(scalar_2d_field), intent(inout) :: scalar_hpressure_at_pc_face_level
     type(scalar_2d_field), intent(inout) :: scalar_delhp_at_pc_full_level
     type(scalar_2d_field), intent(inout) :: scalar_hpressure_at_pc_full_level
     type(scalar_2d_field), intent(inout) :: scalar_delhp_at_pc_face_level
! local
     integer(i4)                          :: iv, ilev
#if defined(CHK_INT_RENEW_MASS_STATE) || defined(CHK_ALL)
#include "data_check.inc"
real(r8), allocatable :: ref_scalar_hpressure_at_pc_face_level_f(:, :)
real(r8), allocatable :: ref_scalar_delhp_at_pc_full_level_f(:, :)
real(r8), allocatable :: ref_scalar_hpressure_at_pc_full_level_f(:, :)
real(r8), allocatable :: ref_scalar_delhp_at_pc_face_level_f(:, :)
allocate(ref_scalar_hpressure_at_pc_face_level_f, source=scalar_hpressure_at_pc_face_level%f)
allocate(ref_scalar_delhp_at_pc_full_level_f, source=scalar_delhp_at_pc_full_level%f)
allocate(ref_scalar_hpressure_at_pc_full_level_f, source=scalar_hpressure_at_pc_full_level%f)
allocate(ref_scalar_delhp_at_pc_face_level_f, source=scalar_delhp_at_pc_face_level%f)
#endif 
!
! renew mass state, compute hpressure at face level, full level, 
! delhp, based on input surface hpressure
!
call t_startf("time_integration_renew_mass_state_target")
!$omp target
!$omp parallel  private(iv,scalar_template_a,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_a,scalar_template_1d_nlev_b,ilev) 
!$omp do
        do iv = 1, mesh%nv
           scalar_template_a   =   scalar_hpressure_at_pc_surface%f(iv)

           call calc_hpe_hpressure_face_level(scalar_template_a, &       ! in
                                             scalar_template_1d_nlevp_a) ! out

           call calc_hpe_delhp(scalar_template_1d_nlevp_a, &  ! in
                               scalar_template_1d_nlev_a)     ! out

           call calc_hpe_hpressure_full_level(scalar_template_a         , & ! in
                                              scalar_template_1d_nlevp_a, & ! in
                                              scalar_template_1d_nlev_a , & ! in
                                              scalar_template_1d_nlev_b)    ! out

           scalar_hpressure_at_pc_face_level%f(:,iv)  = scalar_template_1d_nlevp_a
               scalar_delhp_at_pc_full_level%f(:,iv)  = scalar_template_1d_nlev_a
           scalar_hpressure_at_pc_full_level%f(:,iv)  = scalar_template_1d_nlev_b

           do ilev = 2, nlev
              scalar_delhp_at_pc_face_level%f(ilev,iv)= scalar_hpressure_at_pc_full_level%f(ilev,iv)-&
                                                        scalar_hpressure_at_pc_full_level%f(ilev-1,iv)
           end do

           scalar_delhp_at_pc_face_level%f(1,iv)      = 0.5_r8*scalar_delhp_at_pc_full_level%f(1,iv)
           scalar_delhp_at_pc_face_level%f(nlev+1,iv) = 0.5_r8*scalar_delhp_at_pc_full_level%f(nlev,iv)

        end do
!$omp end do nowait
!$omp end parallel
!$omp end target
call t_stopf("time_integration_renew_mass_state_target")


#if defined(CHK_INT_RENEW_MASS_STATE) || defined(CHK_ALL)
call t_startf("time_integration_renew_mass_state")
!$omp parallel  private(iv,scalar_template_a,scalar_template_1d_nlevp_a,scalar_template_1d_nlev_a,scalar_template_1d_nlev_b,ilev) 
!$omp do
do iv = 1, mesh%nv
        scalar_template_a   =   scalar_hpressure_at_pc_surface%f(iv)

        call calc_hpe_hpressure_face_level(scalar_template_a, &       ! in
                                          scalar_template_1d_nlevp_a) ! out

        call calc_hpe_delhp(scalar_template_1d_nlevp_a, &  ! in
                            scalar_template_1d_nlev_a)     ! out

        call calc_hpe_hpressure_full_level(scalar_template_a         , & ! in
                                           scalar_template_1d_nlevp_a, & ! in
                                           scalar_template_1d_nlev_a , & ! in
                                           scalar_template_1d_nlev_b)    ! out

        ref_scalar_hpressure_at_pc_face_level_f(:,iv)  = scalar_template_1d_nlevp_a
            ref_scalar_delhp_at_pc_full_level_f(:,iv)  = scalar_template_1d_nlev_a
        ref_scalar_hpressure_at_pc_full_level_f(:,iv)  = scalar_template_1d_nlev_b

        do ilev = 2, nlev
           ref_scalar_delhp_at_pc_face_level_f(ilev,iv)= scalar_hpressure_at_pc_full_level%f(ilev,iv)-&
                                                     scalar_hpressure_at_pc_full_level%f(ilev-1,iv)
        end do

        ref_scalar_delhp_at_pc_face_level_f(1,iv)      = 0.5_r8*scalar_delhp_at_pc_full_level%f(1,iv)
        ref_scalar_delhp_at_pc_face_level_f(nlev+1,iv) = 0.5_r8*scalar_delhp_at_pc_full_level%f(nlev,iv)

     end do
!$omp end do nowait
!$omp end parallel
call data_check(ref_scalar_hpressure_at_pc_face_level_f, scalar_hpressure_at_pc_face_level%f)
call data_check(ref_scalar_hpressure_at_pc_full_level_f, scalar_hpressure_at_pc_full_level%f)
call data_check(ref_scalar_delhp_at_pc_full_level_f, scalar_delhp_at_pc_full_level%f)
call data_check(ref_scalar_delhp_at_pc_face_level_f, scalar_delhp_at_pc_face_level%f)
deallocate(ref_scalar_hpressure_at_pc_face_level_f)
deallocate(ref_scalar_hpressure_at_pc_full_level_f)
deallocate(ref_scalar_delhp_at_pc_full_level_f)
deallocate(ref_scalar_delhp_at_pc_face_level_f)
call t_stopf("time_integration_renew_mass_state")

#endif
        return
   end subroutine time_integration_renew_mass_state

  end module grist_dycore_time_integration_2d
