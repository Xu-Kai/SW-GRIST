
!----------------------------------------------------------------------------
! Copyright: Copyright (c) 2016-2019
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!    http://www.apache.org/licenses/LICENSE-2.0
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
! Created on 2018
! Version 1.0
! Description: Global-Regional Integrated forecast SysTem-SWM (GRIST-SWM)
!          Major Program driver of GRIST-A, jointly developed by Chinese
!          Academy of Meteorological Sciences and National Supercomputing
!          Center in Wuxi. Unless otherwise mentioned, this model follows
!          the above License.
!
! Revision history:
!----------------------------------------------------------------------------

 PROGRAM GRIST_swe

   use grist_lib
   use grist_mpi
   use grist_destruct_module
   use grist_domain_types,   only: global_domain,block_structure, global_domain_data
   use grist_nml_module,     only: set_global_vars,lbdry_flag
   use grist_grid_file_read, only: grist_read_grid_file
   use grist_config_mesh,    only: config_mesh_from_file
   use grist_mesh_weight_icosh,only: calc_trsk_weights_plg , &
                                   init_blocal_prime_cell    , &
                                   init_blocal_prime_cell_4th, &
                                   init_blocal_dual_cell,      &
                                   init_blocal_dual_cell_mbs,  &
                                   init_blocal_dual_cell_4th,  &
                                   init_blocal_dual_cell_4th_mbs, &
                                   init_plg_sub_triangle_info
   use grist_geometric_info_icosh, only: init_geometric_info
   use grist_config_partition,only: config_sub_domain , &
                                    domain_decompse , &
                                    show_basic_info
   use grist_nml_module,      only: nsteps      , &
                                    time_scheme , &
                                    outdir      , &
                                    run_type    , &
                                    h1_restart_freq, &
                                    !output_freq, &
                                    test_dual_cell, &
                                    use_tr_mbs  , &
                                    stencil_width, &
                                    stencil_exchange_flag
   use swe_model_driver,      only: swe_model_init, swe_model_run, swe_model_final

   implicit none
   
       type(block_structure) ,target         :: local_block
       integer      :: iv, comm, ierr
       logical      :: flag_rank = .true.
       type(global_domain_data) :: mesh

       comm = MPI_COMM_WORLD
       call mpi_init(ierr)
   
       call set_global_vars()
       call grist_read_grid_file(mesh, comm)
       call config_mesh_from_file(mesh)
       call show_basic_info(mesh, flag_rank)
 
!================================================
!       init global and lifetime variables
!================================================
!
! domain decomposision
!
       call domain_decompse(mesh, comm)
       call config_sub_domain(mesh, local_block, comm)
       call destruct_global_mesh_var(mesh)
       local_block%full_domain%local_block => local_block
!
! init mesh info 
!
       call calc_trsk_weights_plg(local_block%full_domain)
       call init_blocal_prime_cell(local_block%full_domain)
       call init_plg_sub_triangle_info(local_block%full_domain)
       call init_blocal_prime_cell_4th(local_block%full_domain)
       call init_geometric_info(local_block%full_domain)

       if(use_tr_mbs)then
         call init_blocal_dual_cell_mbs(local_block%full_domain)
         call init_blocal_dual_cell_4th_mbs(local_block%full_domain)
       else
         call init_blocal_dual_cell(local_block%full_domain)
         call init_blocal_dual_cell_4th(local_block%full_domain)
       end if
!
! init model instance
!
       call swe_model_init(local_block%full_domain)
   
       if(mpi_rank() == 0) then 
          print*,"MPI Pros Num:  ",mpi_size()
          print*,"low order bdry flag:  ",lbdry_flag
       end if
       if(stencil_width >= 4)then
          stencil_exchange_flag = .false.
       if(mpi_rank() == 0) print*,"stencil tpye:  local    halo:",stencil_width
       else
         stencil_exchange_flag = .true.
       if(mpi_rank() == 0) print*,"stencil tpye:  exchange    halo:",stencil_width
       end if
 
       !call tic()
       call swe_model_run(local_block%full_domain)
       !call toc()
       call swe_model_final(local_block%full_domain)
   
       call mpi_finalize(ierr)

 END PROGRAM GRIST_swe
