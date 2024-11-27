
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Main Driver for the SWE model
!          1) Set global vars
!          2) Initial conditions
!          3) Time_integration
!          4) Diagnose
!          5) File I/O
! Revision history: 
!----------------------------------------------------------------------------

 module swe_model_driver

   use grist_constants,     only: i4, r8
   use grist_domain_types,  only: global_domain
   use grist_nml_module,    only: swe_timestep, &
                                  nsteps      , &
                                  time_scheme , &
                                  outdir      , &
                                  run_type    , &
                                  h1_restart_freq, &
                                  h1_history_freq, &
                                  test_dual_cell, &
                                  use_tr_mbs
   use grist_geometric_info_icosh,only: init_geometric_info
   use grist_mesh_weight_icosh,   only: calc_trsk_weights_plg , &
                                  init_blocal_prime_cell    , &
                                  init_blocal_prime_cell_4th, &
                                  init_blocal_dual_cell,      &
                                  init_blocal_dual_cell_mbs,  &
                                  init_blocal_dual_cell_4th,  &
                                  init_blocal_dual_cell_4th_mbs, &
                                  init_plg_sub_triangle_info

   use swe_vars_module
   use swe_init_module,     only: swe_initialize, swe_initialize_tr

   use swe_time_integration,only: time_integration_rk       ,&
                                  time_integration_fit      ,&
                                  time_integration_rk3_tr   ,& 
                                  time_integration_fit_tr 

   use swe_diagnose_module, only: diagnose_total_energy     ,&
                                  diagnose_hkinetic_energy  ,&
                                  diagnose_potential_enstropy

   use swe_restart_module,  only: restart_timestep          ,&
                                  swe_read_restart_file     ,&
                                  swe_write_restart_file

   use swe_inout_module,    only: swe_output_file
   use grist_config_partition,  only: debug_data_1d
   use grist_mpi

   implicit none

   private

   public   :: swe_model_init, &
               swe_model_run,  &
               swe_model_final
 
   contains

  subroutine swe_model_init(mesh)
   use swe_inout_module,    only: swe_output_init
   use swe_vars_module,     only: swe_vars_init
   use swe_init_module,     only: swe_initialize, swe_initialize_tr

   type(global_domain), intent(inout) :: mesh

      call swe_vars_init(mesh)
 
      call swe_initialize(mesh,scalar_normal_velocity_at_edge,&
                               scalar_height_at_prime_cell   ,&
                               scalar_topo_at_prime_cell,     &
                               ue_init, ve_init)

      if(test_dual_cell)then
         call swe_initialize_tr(mesh,scalar_normal_velocity_at_tr_edge,&
                                     scalar_height_at_dual_cell       ,&
                                     ue_init, ve_init)
      end if
      call swe_output_init(mesh)

      mesh%nv = mesh%nv_compute
      mesh%ne = mesh%ne_compute
      mesh%nt = mesh%nt_compute
    return
  end subroutine swe_model_init

  subroutine  swe_model_final(mesh)
   use swe_inout_module,    only: swe_output_final
   type(global_domain), intent(inout) :: mesh
     call swe_output_final(mesh)
  end subroutine swe_model_final

!================================================
!          Major stepping of the model
!================================================

  subroutine swe_model_run(mesh)
! io
   type(global_domain), intent(inout) :: mesh
! local
   real(r8)             :: dtime
   real(r8)             :: time_beg
   real(r8)             :: time_end
   real(r8)             :: time_elapse
   real(r8)             :: potential_enstropy
   integer(i4)          :: itimestep,ierr
   integer(i4)          :: it
   integer(i4)          :: ie
   integer(i4)          :: iv
   integer(i4)          :: istep_beg
   integer(i4)          :: istep_end
   logical              :: input_tend
   real(r8)             :: sum_time,max_time

!================================================
!       init global and lifetime variables
!================================================

!     call swe_vars_init(mesh)
!     call calc_trsk_weights_plg(mesh)

!     call init_blocal_prime_cell(mesh)
!     call init_plg_sub_triangle_info(mesh)
!     call init_blocal_prime_cell_4th(mesh)
!     call init_geometric_info(mesh)

     if(use_tr_mbs)then
!        call init_blocal_dual_cell_mbs(mesh)
!        call init_blocal_dual_cell_4th_mbs(mesh)
     else
!        call init_blocal_dual_cell(mesh)
!        call init_blocal_dual_cell_4th(mesh)
     end if

!================================================
!     initial conditions and topography
!================================================

!     call swe_initialize(mesh,scalar_normal_velocity_at_edge,&
!                              scalar_height_at_prime_cell   ,&
!                              scalar_topo_at_prime_cell,     &
!                              ue_init, ve_init)

     if(test_dual_cell)then
!          call swe_initialize_tr(mesh,scalar_normal_velocity_at_tr_edge,&
!                                      scalar_height_at_dual_cell       ,&
!                                      ue_init, ve_init)
     end if

     istep_beg  = 1
     istep_end  = nsteps
     scalar_hb_at_prime_cell%f(:)   = scalar_height_at_prime_cell%f(:)+&
                                      scalar_topo_at_prime_cell%f(:)

     if(trim(run_type).eq.'restart')then ! overwrite data and istep_beg
        call swe_read_restart_file(mesh)
        istep_beg = restart_timestep
     end if

     dtime = swe_timestep

!================================================
!           WRITE INITIAL CONDITIONS
!================================================

     IF(trim(run_type).eq.'init') call swe_output_file(mesh,0)
 
!================================================
!              major time marching
!================================================

     call cpu_time(time_beg)

   if(.not.test_dual_cell)then

   DO itimestep = istep_beg, istep_end

      if(mpi_rank() == 0) print*,"CURRENT TIMESTEP=",ITIMESTEP,":",istep_end
!
! Diagnose State
!
      call diagnose_total_energy(mesh)
      call diagnose_hkinetic_energy(mesh)
      call diagnose_potential_enstropy(mesh)

!
! Time integration
!

    select case(trim(time_scheme))


!     case('rkfb')
!       call time_integration_rkfb(mesh, scalar_normal_velocity_at_edge,     & ! input
!                                        scalar_height_at_prime_cell,        & ! input
!                                        scalar_topo_at_prime_cell,          & ! input
!                                        scalar_normal_velocity_at_edge_next,& ! output
!                                        scalar_height_at_prime_cell_next,   & ! output
!                                        dtime)
! !
! ! time-split version
! !
!     case('splexp1')
!       call time_integration_splexp1(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                          scalar_height_at_prime_cell,        & ! input
!                                          scalar_topo_at_prime_cell,          & ! input
!                                          scalar_normal_velocity_at_edge_next,& ! output
!                                          scalar_height_at_prime_cell_next,   & ! output
!                                          dtime)
!     case('splexp2')
!       call time_integration_splexp2(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                          scalar_height_at_prime_cell,        & ! input
!                                          scalar_topo_at_prime_cell,          & ! input
!                                          scalar_normal_velocity_at_edge_next,& ! output
!                                          scalar_height_at_prime_cell_next,   & ! output
!                                          dtime)
!     case('splexp3')
!       call time_integration_splexp3(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                          scalar_height_at_prime_cell,        & ! input
!                                          scalar_topo_at_prime_cell,          & ! input
!                                          scalar_normal_velocity_at_edge_next,& ! output
!                                          scalar_height_at_prime_cell_next,   & ! output
!                                          dtime)

!     case('splexp5')
!       call time_integration_splexp5(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                          scalar_height_at_prime_cell,        & ! input
!                                          scalar_topo_at_prime_cell,          & ! input
!                                          scalar_normal_velocity_at_edge_next,& ! output
!                                          scalar_height_at_prime_cell_next,   & ! output
!                                          dtime)
! !
! ! the same time integration for the whole system
! ! i.e., no split
! !

!     case('rk4')
!       call time_integration_rk4(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                      scalar_height_at_prime_cell,        & ! input
!                                      scalar_topo_at_prime_cell,          & ! input
!                                      scalar_normal_velocity_at_edge_next,& ! output
!                                      scalar_height_at_prime_cell_next,   & ! output
!                                      tend_normal_velocity_at_edge,       & ! output
!                                      tend_height_at_prime_cell,          & ! output
!                                      itimestep,dtime)
    case('rk3','rk2','fit')
      call time_integration_rk(mesh,scalar_normal_velocity_at_edge,     & ! input
                                    scalar_height_at_prime_cell,        & ! input
                                    scalar_topo_at_prime_cell,          & ! input
                                    scalar_normal_velocity_at_edge_next,& ! output
                                    scalar_height_at_prime_cell_next,   & ! output
                                    tend_normal_velocity_at_edge,       & ! output
                                    tend_height_at_prime_cell,          & ! output
                                    itimestep, dtime)
!     case('pc')
!       call time_integration_pc(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                     scalar_height_at_prime_cell,        & ! input
!                                     scalar_topo_at_prime_cell,          & ! input
!                                     scalar_normal_velocity_at_edge_next,& ! output
!                                     scalar_height_at_prime_cell_next,   & ! output
!                                     tend_normal_velocity_at_edge,       & ! output
!                                     tend_height_at_prime_cell,          & ! output
!                                     itimestep, dtime)
!    case('fit')
!      call time_integration_fit(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                     scalar_height_at_prime_cell,        & ! input
!                                     scalar_topo_at_prime_cell,          & ! input
!                                     scalar_normal_velocity_at_edge_next,& ! output
!                                     scalar_height_at_prime_cell_next,   & ! output
!                                     tend_normal_velocity_at_edge,       & ! output
!                                     tend_height_at_prime_cell,          & ! output
!                                     itimestep,dtime)
!     case('ab3')

!     if(itimestep.le.2) then
!       call time_integration_rk4(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                      scalar_height_at_prime_cell,        & ! input
!                                      scalar_topo_at_prime_cell,          & ! input
!                                      scalar_normal_velocity_at_edge_next,& ! output
!                                      scalar_height_at_prime_cell_next,   & ! output
!                                      tend_normal_velocity_at_edge,       & ! output
!                                      tend_height_at_prime_cell,          & ! output
!                                      itimestep,dtime)
! ! for next step
!        scalar_height_at_prime_cell_nm2%f(:)    = scalar_height_at_prime_cell_nm1%f(:)
!        scalar_normal_velocity_at_edge_nm2%f(:) = scalar_normal_velocity_at_edge_nm1%f(:)
!        scalar_height_at_prime_cell_nm1%f(:)    = scalar_height_at_prime_cell%f(:)
!        scalar_normal_velocity_at_edge_nm1%f(:) = scalar_normal_velocity_at_edge%f(:)

!     else if(itimestep.gt.2.and.itimestep.lt.5)then
!        input_tend = .false.
!        call time_integration_ab3(mesh,scalar_normal_velocity_at_edge_nm2,  & ! n-2
!                                       scalar_normal_velocity_at_edge_nm1,  & ! n-1
!                                       scalar_normal_velocity_at_edge,      & ! n
!                                       scalar_height_at_prime_cell_nm2,     & ! n-2
!                                       scalar_height_at_prime_cell_nm1,     & ! n-1
!                                       scalar_height_at_prime_cell,         & ! n 
!                                       scalar_topo_at_prime_cell,           & ! n 
!                                       scalar_normal_velocity_at_edge_next, & ! out
!                                       scalar_height_at_prime_cell_next,    & ! out
!                                       tend_normal_velocity_at_edge,        & ! out
!                                       tend_height_at_prime_cell,           & ! out
!                                       itimestep, dtime, input_tend)
! ! for next step
!        tend_height_at_prime_cell_nm2%f(:)      = tend_height_at_prime_cell_nm1%f(:)
!        tend_normal_velocity_at_edge_nm2%f(:)   = tend_normal_velocity_at_edge_nm1%f(:)
!        tend_height_at_prime_cell_nm1%f(:)      = tend_height_at_prime_cell%f(:)
!        tend_normal_velocity_at_edge_nm1%f(:)   = tend_normal_velocity_at_edge%f(:)
!        scalar_height_at_prime_cell_nm2%f(:)    = scalar_height_at_prime_cell_nm1%f(:)
!        scalar_normal_velocity_at_edge_nm2%f(:) = scalar_normal_velocity_at_edge_nm1%f(:)
!        scalar_height_at_prime_cell_nm1%f(:)    = scalar_height_at_prime_cell%f(:)
!        scalar_normal_velocity_at_edge_nm1%f(:) = scalar_normal_velocity_at_edge%f(:)

!     else if(itimestep.ge.5)then ! from 5th step to end 
!        input_tend = .true.
!        call time_integration_ab3(mesh,tend_normal_velocity_at_edge_nm2,    & ! n-2
!                                       tend_normal_velocity_at_edge_nm1,    & ! n-1
!                                       scalar_normal_velocity_at_edge,      & ! n
!                                       tend_height_at_prime_cell_nm2,       & ! n-2
!                                       tend_height_at_prime_cell_nm1,       & ! n-1
!                                       scalar_height_at_prime_cell,         & ! n 
!                                       scalar_topo_at_prime_cell,           & ! n 
!                                       scalar_normal_velocity_at_edge_next, & ! out
!                                       scalar_height_at_prime_cell_next,    & ! out
!                                       tend_normal_velocity_at_edge,        & ! out
!                                       tend_height_at_prime_cell,           & ! out
!                                       itimestep, dtime, input_tend)
! ! for next step
!        tend_height_at_prime_cell_nm2%f(:)    = tend_height_at_prime_cell_nm1%f(:)
!        tend_normal_velocity_at_edge_nm2%f(:) = tend_normal_velocity_at_edge_nm1%f(:)
!        tend_height_at_prime_cell_nm1%f(:)    = tend_height_at_prime_cell%f(:)
!        tend_normal_velocity_at_edge_nm1%f(:) = tend_normal_velocity_at_edge%f(:)
!     end if

!     case('ssprk')
!       call time_integration_ssprk(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                        scalar_height_at_prime_cell,        & ! input
!                                        scalar_topo_at_prime_cell,          & ! input
!                                        scalar_normal_velocity_at_edge_next,& ! output
!                                        scalar_height_at_prime_cell_next,   & ! output
!                                        tend_normal_velocity_at_edge,       & ! output
!                                        tend_height_at_prime_cell,          & ! output
!                                        itimestep,dtime)

!     case('abmpc')

!     if(itimestep.eq.1) then
!       call time_integration_rk4(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                      scalar_height_at_prime_cell,        & ! input
!                                      scalar_topo_at_prime_cell,          & ! input
!                                      scalar_normal_velocity_at_edge_next,& ! output
!                                      scalar_height_at_prime_cell_next,   & ! output
!                                      tend_normal_velocity_at_edge,       & ! output
!                                      tend_height_at_prime_cell,          & ! output
!                                      itimestep,dtime)
! ! for next step
!        scalar_height_at_prime_cell_nm1%f(:)    = scalar_height_at_prime_cell%f(:)
!        scalar_normal_velocity_at_edge_nm1%f(:) = scalar_normal_velocity_at_edge%f(:)

!     else if(itimestep.eq.2)then
!        input_tend = .false.
!        call time_integration_abmpc(mesh,scalar_normal_velocity_at_edge_nm1,  & ! n-1
!                                         scalar_normal_velocity_at_edge,      & ! n
!                                         scalar_height_at_prime_cell_nm1,     & ! n-1
!                                         scalar_height_at_prime_cell,         & ! n 
!                                         scalar_topo_at_prime_cell,           & ! n 
!                                         scalar_normal_velocity_at_edge_next, & ! out
!                                         scalar_height_at_prime_cell_next,    & ! out
!                                         tend_normal_velocity_at_edge,        & ! out
!                                         tend_height_at_prime_cell,           & ! out
!                                         itimestep, dtime, input_tend)
! ! for next step
!        tend_height_at_prime_cell_nm1%f(:)      = tend_height_at_prime_cell%f(:)
!        tend_normal_velocity_at_edge_nm1%f(:)   = tend_normal_velocity_at_edge%f(:)
       
!     else if(itimestep.ge.3)then ! from 3th step to end
!        input_tend = .true.
!        call time_integration_abmpc(mesh,tend_normal_velocity_at_edge_nm1,    & ! n-1
!                                         scalar_normal_velocity_at_edge,      & ! n
!                                         tend_height_at_prime_cell_nm1,       & ! n-1
!                                         scalar_height_at_prime_cell,         & ! n
!                                         scalar_topo_at_prime_cell,           & ! n
!                                         scalar_normal_velocity_at_edge_next, & ! out
!                                         scalar_height_at_prime_cell_next,    & ! out
!                                         tend_normal_velocity_at_edge,        & ! out
!                                         tend_height_at_prime_cell,           & ! out
!                                         itimestep, dtime, input_tend)
! ! for next step
!        tend_height_at_prime_cell_nm1%f(:)    = tend_height_at_prime_cell%f(:)
!        tend_normal_velocity_at_edge_nm1%f(:) = tend_normal_velocity_at_edge%f(:)
!     end if

!     case('lftz') ! leapfrog trapezoidal
        
!     if(itimestep.eq.1) then
!       call time_integration_rk4(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                      scalar_height_at_prime_cell,        & ! input
!                                      scalar_topo_at_prime_cell,          & ! input
!                                      scalar_normal_velocity_at_edge_next,& ! output
!                                      scalar_height_at_prime_cell_next,   & ! output
!                                      tend_normal_velocity_at_edge,       & ! output
!                                      tend_height_at_prime_cell,          & ! output
!                                      itimestep,dtime)
! ! for next step
!        scalar_height_at_prime_cell_nm1%f(:)    = scalar_height_at_prime_cell%f(:)
!        scalar_normal_velocity_at_edge_nm1%f(:) = scalar_normal_velocity_at_edge%f(:)

!     else if(itimestep.gt.2)then
!        call time_integration_lftz(mesh,scalar_normal_velocity_at_edge_nm1,  & ! n-1
!                                        scalar_normal_velocity_at_edge,      & ! n
!                                        scalar_height_at_prime_cell_nm1,     & ! n-1
!                                        scalar_height_at_prime_cell,         & ! n 
!                                        scalar_topo_at_prime_cell,           & ! n 
!                                        scalar_normal_velocity_at_edge_next, & ! out
!                                        scalar_height_at_prime_cell_next,    & ! out
!                                        tend_normal_velocity_at_edge,        & ! out
!                                        tend_height_at_prime_cell,           & ! out
!                                        itimestep, dtime)
! ! for next step
!        scalar_height_at_prime_cell_nm1%f(:)    = scalar_height_at_prime_cell%f(:)
!        scalar_normal_velocity_at_edge_nm1%f(:) = scalar_normal_velocity_at_edge%f(:)
!     end if

!     case('lf') ! leapfrog
        
!     if(itimestep.eq.1) then
!       call time_integration_rk4(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                      scalar_height_at_prime_cell,        & ! input
!                                      scalar_topo_at_prime_cell,          & ! input
!                                      scalar_normal_velocity_at_edge_next,& ! output
!                                      scalar_height_at_prime_cell_next,   & ! output
!                                      tend_normal_velocity_at_edge,       & ! output
!                                      tend_height_at_prime_cell,          & ! output
!                                      itimestep,dtime)
! ! for next step
!        scalar_height_at_prime_cell_nm1%f(:)    = scalar_height_at_prime_cell%f(:)
!        scalar_normal_velocity_at_edge_nm1%f(:) = scalar_normal_velocity_at_edge%f(:)

!     else if(itimestep.gt.2)then
!        call time_integration_lf(mesh,scalar_normal_velocity_at_edge_nm1,  & ! n-1
!                                      scalar_normal_velocity_at_edge,      & ! n
!                                      scalar_height_at_prime_cell_nm1,     & ! n-1
!                                      scalar_height_at_prime_cell,         & ! n 
!                                      scalar_topo_at_prime_cell,           & ! n 
!                                      scalar_normal_velocity_at_edge_next, & ! out
!                                      scalar_height_at_prime_cell_next,    & ! out
!                                      tend_normal_velocity_at_edge,        & ! out
!                                      tend_height_at_prime_cell,           & ! out
!                                      itimestep, dtime)
! ! for next step
!        scalar_height_at_prime_cell_nm1%f(:)    = scalar_height_at_prime_cell%f(:)
!        scalar_normal_velocity_at_edge_nm1%f(:) = scalar_normal_velocity_at_edge%f(:)
!     end if

!     case('asselin') ! asselin leapfrog
        
!     if(itimestep.eq.1) then
!       call time_integration_rk4(mesh,scalar_normal_velocity_at_edge,     & ! input
!                                      scalar_height_at_prime_cell,        & ! input
!                                      scalar_topo_at_prime_cell,          & ! input
!                                      scalar_normal_velocity_at_edge_next,& ! output
!                                      scalar_height_at_prime_cell_next,   & ! output
!                                      tend_normal_velocity_at_edge,       & ! output
!                                      tend_height_at_prime_cell,          & ! output
!                                      itimestep,dtime)
! ! for next step
!        scalar_height_at_prime_cell_nm1%f(:)    = scalar_height_at_prime_cell%f(:)
!        scalar_normal_velocity_at_edge_nm1%f(:) = scalar_normal_velocity_at_edge%f(:)

!     else if(itimestep.ge.2)then
!        call time_integration_asselin(mesh,scalar_normal_velocity_at_edge_nm1,  & ! n-1
!                                           scalar_normal_velocity_at_edge,      & ! n
!                                           scalar_height_at_prime_cell_nm1,     & ! n-1
!                                           scalar_height_at_prime_cell,         & ! n 
!                                           scalar_topo_at_prime_cell,           & ! n 
!                                           scalar_normal_velocity_at_edge_next, & ! out
!                                           scalar_height_at_prime_cell_next,    & ! out
!                                           tend_normal_velocity_at_edge,        & ! out
!                                           tend_height_at_prime_cell,           & ! out
!                                           itimestep, dtime)
! ! for next step
!        scalar_height_at_prime_cell_nm1%f(:)    = scalar_height_at_prime_cell%f(:)
!        scalar_normal_velocity_at_edge_nm1%f(:) = scalar_normal_velocity_at_edge%f(:)

!    end if
    case default
      print*," you must select a time scheme"
      stop
    end select
!
! renew for the next step
!
      scalar_height_at_prime_cell    = scalar_height_at_prime_cell_next
      scalar_normal_velocity_at_edge = scalar_normal_velocity_at_edge_next
      scalar_hb_at_prime_cell%f(:)   = scalar_height_at_prime_cell%f(:)+&
                                       scalar_topo_at_prime_cell%f(:)
!
! I/O 
!
      if(itimestep==1 .or. itimestep==nsteps .or. mod(itimestep,h1_history_freq)==0 )then
         call swe_output_file(mesh,itimestep)
      end if

      if(itimestep==1 .or. itimestep==nsteps .or. mod(itimestep,h1_restart_freq)==0 )then
         !call swe_write_restart_file(mesh,(itimestep+1))
      end if
  
      !check output values, debug
      !if(mod(itimestep,1000) .eq. 0)then
      if(itimestep .eq. int(nsteps) .and. .false.)then
      
        call debug_data_1d(itimestep,mesh%v_index,mesh%nv,"PC_AREA",scalar_area_of_prime_cell%f)
        call debug_data_1d(itimestep,mesh%t_index,mesh%nt,"DC_AREA",scalar_area_of_dual_cell%f) 
        call debug_data_1d(itimestep,mesh%e_index,mesh%ne,"HX_EDGE",scalar_leng_of_hx_edge%f) 
        call debug_data_1d(itimestep,mesh%e_index,mesh%ne,"TR_EDGE",scalar_leng_of_tr_edge%f) 
        call debug_data_1d(itimestep,mesh%v_index,mesh%nv,"TOPO",scalar_topo_at_prime_cell%f) 
        call debug_data_1d(itimestep,mesh%t_index,mesh%nt,"H_DC",scalar_height_at_dual_cell%f)

        call debug_data_1d(itimestep,mesh%v_index,mesh%nv,"H",scalar_height_at_prime_cell%f) !H ~
        call debug_data_1d(itimestep,mesh%v_index,mesh%nv,"HB",scalar_hb_at_prime_cell%f) !HB ~
        call debug_data_1d(itimestep,mesh%e_index,mesh%ne,"UE",ue_reconst%f) !UE
        call debug_data_1d(itimestep,mesh%e_index,mesh%ne,"VE",ve_reconst%f) !VE
        call debug_data_1d(itimestep,mesh%v_index,mesh%nv,"UP",up_reconst%f) !UP
        call debug_data_1d(itimestep,mesh%v_index,mesh%nv,"VP",vp_reconst%f) !VP
        call debug_data_1d(itimestep,mesh%v_index,mesh%nv,"DIV",scalar_divergence_at_prime_cell%f) !DIV
        call debug_data_1d(itimestep,mesh%t_index,mesh%nt,"AVOR",scalar_absolute_vorticity_at_dual_cell%f)  !AVOR
        call debug_data_1d(itimestep,mesh%t_index,mesh%nt,"RVOR",scalar_relative_vorticity_at_dual_cell%f)  !RVOR
        call debug_data_1d(itimestep,mesh%t_index,mesh%nt,"PVOR",scalar_potential_vorticity_at_dual_cell%f) !PVOR
      end if
      !stop
  

   END DO

   end if  ! 
   
   if(mpi_rank() == 0) print *,"rk3 is finish, test_dual_cell:",test_dual_cell

   if(test_dual_cell)then
      DO itimestep = istep_beg, istep_end

         print*,"CURRENT TIMESTEP=",ITIMESTEP

         select case(trim(time_scheme))
         case('rk3')
            call time_integration_rk3_tr(mesh,scalar_normal_velocity_at_tr_edge,     & ! input
                                              scalar_height_at_dual_cell,            & ! input
                                              scalar_normal_velocity_at_tr_edge_next,& ! output
                                              scalar_height_at_dual_cell_next,       & ! output
                                              itimestep,dtime)
         case('fit')
            call time_integration_fit_tr(mesh,scalar_normal_velocity_at_tr_edge,     & ! input
                                              scalar_height_at_dual_cell,            & ! input
                                              scalar_normal_velocity_at_tr_edge_next,& ! output
                                              scalar_height_at_dual_cell_next,       & ! output
                                              itimestep,dtime)
         case default
            print*, "you must select a time scheme for tr advection"
            stop
         end select
         scalar_height_at_dual_cell    = scalar_height_at_dual_cell_next
         if(itimestep==1 .or. itimestep==nsteps .or. mod(itimestep,h1_history_freq)==0 )then
            call swe_output_file(mesh,itimestep)
         end if
      END DO
   end if

   call cpu_time(time_end)

   time_elapse = time_end-time_beg

   call reduce(time_elapse, sum_time, 'sum')

   sum_time = sum_time/mpi_size()

   call reduce(time_elapse, max_time, 'max')

   if(mpi_rank() == 0)then
     print *,"average elapsed time:",sum_time,"  secs"
     print *,"max elapsed time:",max_time,"  secs"
     print *,"number of processors:", mpi_size()
     open(1,file=trim(outdir)//"/elapsed_time.txt",status='unknown')
     write(1,*) "max time_elapsed=", max_time,"  secs"
     write(1,*) "average time_elapsed=", sum_time,"  secs"
     write(1,*) "number of processors=", mpi_size()
     close(1)
   end if

   call swe_vars_clean

   return
  end subroutine swe_model_run

  end module swe_model_driver
