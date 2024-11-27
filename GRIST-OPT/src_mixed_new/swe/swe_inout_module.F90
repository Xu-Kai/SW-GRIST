
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Input and output file. Currently only for output because IC 
!              is coded
! Revision history: 
!----------------------------------------------------------------------------

  module swe_inout_module

  use grist_lib
  use grist_constants,      only: i4
  use grist_domain_types,   only: global_domain, group_comm, block_structure

  use grist_nml_module,   only: isolate_advection_test,&
                                fname_output,          &
                                time_scheme,           &
                                advection_scheme,      &
                                conserve_scheme,       &
                                swe_timestep,          &
                                testcase,              &
                                initialfield,          &
                                pv_order,              &
                                outdir,                &
                                comm_group_size

  use grist_time_manager,  only: get_current_date

  use grist_recon_module,  only: project_uv, vector_recon_perot_edge2cell_uv
  use grist_util_module,   only: wind_edge_to_cell

  use grist_fileio_list_1d_module_par, only: wrap_output_init_1d,&
                                             wrap_add_field_1d  ,&
                                             wrap_output_1d_group,&
                                             wrap_output_clean_1d

  use swe_vars_module,    only: scalar_normal_velocity_at_edge        , &
                                scalar_height_at_prime_cell           , &
                                scalar_hb_at_prime_cell               , &
                                scalar_area_of_prime_cell             , & 
                                scalar_area_of_dual_cell              , &
                                scalar_leng_of_hx_edge                , &
                                scalar_leng_of_tr_edge                , &
                                lon_nt,lat_nt,lon_ne,lat_ne,lon_nv,lat_nv, &
                                scalar_topo_at_prime_cell             , &
                                ue_init                               , &
                                ve_init                               , &
                                up_init                               , &
                                vp_init                               , &
                                ue_reconst                            , &
                                ve_reconst                            , &
                                up_reconst                            , &
                                vp_reconst                            , &
                                scalar_divergence_at_prime_cell       , &
                                scalar_absolute_vorticity_at_dual_cell, &
                                scalar_relative_vorticity_at_dual_cell, &
                                scalar_potential_vorticity_at_dual_cell,&
                                scalar_height_at_dual_cell

  use swe_diagnose_module, only: diagnose_divergence_vorticity


  implicit none

  private

  public  ::  swe_output_init, &
              swe_output_file, &
              swe_output_final

  contains

   subroutine swe_output_init(mesh)
    ! io
    type(global_domain),   intent(inout), target :: mesh
    ! local
    integer(i4)                                  :: i, ig, np_output
    integer(i4)                                  :: nv_avg, nv_res, vpos, &
                                                    ne_avg, ne_res, epos, &
                                                    nt_avg, nt_res, tpos
    type(set)                                    :: sv_list, se_list, st_list
    integer,               allocatable           :: buff_svi(:), buff_sti(:), buff_sei(:)

    type(block_structure), pointer               :: local_block
    type(group_comm),      pointer               :: comm
    integer                                      :: nprocs, ierr
    integer, allocatable                         :: recvvi(:), recvti(:), recvei(:)
    integer(i4)                                  :: rvtot, rttot, retot

    local_block => mesh%local_block
    comm => mesh%gcomm_write
    comm%all_comm = local_block%comm
    comm%gsize = comm_group_size
    comm%color = mod(mpi_rank(comm%all_comm),comm%gsize)
    call mpi_comm_split(comm%all_comm, comm%color, mpi_rank(comm%all_comm), comm%gcomm,ierr)
    call mpi_comm_rank(comm%gcomm,comm%grank,ierr)

    nprocs = mpi_size(comm%all_comm)
    np_output = ceiling(nprocs/real(comm_group_size, 8))
    ! for vtx
    ! sort v_index & get localv
    do i = 1, size(local_block%output%v)
        call sv_list%insert(local_block%output%v(i))
    enddo
    allocate(buff_svi(sv_list%size()))
    call sv_list%sort(buff_svi)
    call sv_list%final
    allocate(comm%idx_localv(size(buff_svi)))
    do i = 1, size(buff_svi)
       comm%idx_localv(i) = mesh%map_g2l_v%find(buff_svi(i))
    end do
    ! sendv counts & displs
    allocate(comm%counts_sendv(1:nprocs), source=0)
    allocate(comm%displs_sendv(1:nprocs), source=0)
    nv_avg = mesh%nv_all/np_output
    nv_res = mod(mesh%nv_all, np_output)
    vpos = (nv_avg+1)*nv_res
    do i = 1, size(buff_svi)
        if(buff_svi(i) .gt. vpos) then
            ig = (buff_svi(i)-vpos-1)/nv_avg + nv_res
        else
            ig = (buff_svi(i)-1)/(nv_avg+1)
        end if
        comm%counts_sendv(ig*comm_group_size+1) = comm%counts_sendv(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        comm%displs_sendv(i) = comm%displs_sendv(i-1) + comm%counts_sendv(i-1)
    end do
    ! recvv counts & displs
    allocate(comm%counts_recvv(nprocs), source=0)
    allocate(comm%displs_recvv(nprocs), source=0)
    call MPI_Alltoall(comm%counts_sendv, 1, MPI_INTEGER, comm%counts_recvv, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if(comm%color .eq. 0) then
        rvtot = comm%counts_recvv(1)
        do i = 2, nprocs
            comm%displs_recvv(i) = comm%displs_recvv(i-1) + comm%counts_recvv(i-1)
            rvtot = rvtot + comm%counts_recvv(i)
        end do
    else
        rvtot = 1
    end if
    allocate(recvvi(rvtot))
    call MPI_Alltoallv(buff_svi, comm%counts_sendv, comm%displs_sendv, MPI_INTEGER, recvvi, &
                       comm%counts_recvv, comm%displs_recvv, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    !open(101, file="recvv"//i2s(mpi_rank()))
    !do i = 1, size(recvvi)
    !  write(101, *) recvvi(i)
    !end do
    !close(101)
    ! recvv local idx
    allocate(comm%idx_recvv(rvtot))
    if(comm%color .eq. 0) then
        if(comm%grank .le. nv_res) then
            comm%mystart_v = comm%grank*(nv_avg+1) + 1
        else
            comm%mystart_v = vpos + (comm%grank-nv_res)*nv_avg + 1
        end if
        if(comm%grank .lt. nv_res) then
            comm%mycount_v = nv_avg + 1
        else
            comm%mycount_v = nv_avg
        end if
        do i = 1, rvtot
            comm%idx_recvv(i) = recvvi(i)-comm%mystart_v+1
        end do
    end if
    deallocate(buff_svi)
    deallocate(recvvi)

    ! for edg
    ! sort e_index & get locale
    do i = 1, size(local_block%output%e)
        call se_list%insert(local_block%output%e(i))
    enddo
    allocate(buff_sei(se_list%size()))
    call se_list%sort(buff_sei)
    call se_list%final
    allocate(comm%idx_locale(size(buff_sei)))
    do i = 1, size(buff_sei)
       comm%idx_locale(i) = mesh%map_g2l_e%find(buff_sei(i))
    end do
    ! sende counts & displs
    allocate(comm%counts_sende(1:nprocs), source=0)
    allocate(comm%displs_sende(1:nprocs), source=0)
    ne_avg = mesh%ne_all/np_output
    ne_res = mod(mesh%ne_all, np_output)
    epos = (ne_avg+1)*ne_res
    do i = 1, size(buff_sei)
        if(buff_sei(i) .gt. epos) then
            ig = (buff_sei(i)-epos-1)/ne_avg + ne_res
        else
            ig = (buff_sei(i)-1)/(ne_avg+1)
        end if
        comm%counts_sende(ig*comm_group_size+1) = comm%counts_sende(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        comm%displs_sende(i) = comm%displs_sende(i-1) + comm%counts_sende(i-1)
    end do
    ! recve counts & displs
    allocate(comm%counts_recve(nprocs), source=0)
    allocate(comm%displs_recve(nprocs), source=0)
    call MPI_Alltoall(comm%counts_sende, 1, MPI_INTEGER, comm%counts_recve, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if(comm%color .eq. 0) then
        retot = comm%counts_recve(1)
        do i = 2, nprocs
            comm%displs_recve(i) = comm%displs_recve(i-1) + comm%counts_recve(i-1)
            retot = retot + comm%counts_recve(i)
        end do
    else
        retot = 1
    end if
    allocate(recvei(retot))
    call MPI_Alltoallv(buff_sei, comm%counts_sende, comm%displs_sende, MPI_INTEGER, recvei, &
                       comm%counts_recve, comm%displs_recve, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    ! recve local idx
    allocate(comm%idx_recve(retot))
    if(comm%color .eq. 0) then
        if(comm%grank .le. ne_res) then
            comm%mystart_e = comm%grank*(ne_avg+1) + 1
        else
            comm%mystart_e = epos + (comm%grank-ne_res)*ne_avg + 1
        end if
        if(comm%grank .lt. ne_res) then
            comm%mycount_e = ne_avg + 1
        else
            comm%mycount_e = ne_avg
        end if
        do i = 1, retot
            comm%idx_recve(i) = recvei(i)-comm%mystart_e+1
        end do
    end if
    deallocate(buff_sei)
    deallocate(recvei)

    ! for tri
    ! sort t_index & get localt
    do i = 1, size(local_block%output%t)
        call st_list%insert(local_block%output%t(i))
    enddo
    allocate(buff_sti(st_list%size()))
    call st_list%sort(buff_sti)
    call st_list%final
    allocate(comm%idx_localt(size(buff_sti)))
    do i = 1, size(buff_sti)
       comm%idx_localt(i) = mesh%map_g2l_t%find(buff_sti(i))
    end do
    ! sendt counts & displs
    allocate(comm%counts_sendt(1:nprocs), source=0)
    allocate(comm%displs_sendt(1:nprocs), source=0)
    nt_avg = mesh%nt_all/np_output
    nt_res = mod(mesh%nt_all, np_output)
    tpos = (nt_avg+1)*nt_res
    do i = 1, size(buff_sti)
        if(buff_sti(i) .gt. tpos) then
            ig = (buff_sti(i)-tpos-1)/nt_avg + nt_res
        else
            ig = (buff_sti(i)-1)/(nt_avg+1)
        end if
        comm%counts_sendt(ig*comm_group_size+1) = comm%counts_sendt(ig*comm_group_size+1) + 1
    enddo
    do i = 2, nprocs
        comm%displs_sendt(i) = comm%displs_sendt(i-1) + comm%counts_sendt(i-1)
    end do
    ! recvt counts & displs
    allocate(comm%counts_recvt(nprocs), source=0)
    allocate(comm%displs_recvt(nprocs), source=0)
    call MPI_Alltoall(comm%counts_sendt, 1, MPI_INTEGER, comm%counts_recvt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if(comm%color .eq. 0) then
        rttot = comm%counts_recvt(1)
        do i = 2, nprocs
            comm%displs_recvt(i) = comm%displs_recvt(i-1) + comm%counts_recvt(i-1)
            rttot = rttot + comm%counts_recvt(i)
        end do
    else
        rttot = 1
    end if
    allocate(recvti(rttot))
    call MPI_Alltoallv(buff_sti, comm%counts_sendt, comm%displs_sendt, MPI_INTEGER, recvti, &
                       comm%counts_recvt, comm%displs_recvt, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    ! recvt local idx
    allocate(comm%idx_recvt(rttot))
    if(comm%color .eq. 0) then
        if(comm%grank .le. nt_res) then
            comm%mystart_t = comm%grank*(nt_avg+1) + 1
        else
            comm%mystart_t = tpos + (comm%grank-nt_res)*nt_avg + 1
        end if
        if(comm%grank .lt. nt_res) then
            comm%mycount_t = nt_avg + 1
        else
            comm%mycount_t = nt_avg
        end if
        do i = 1, rttot
            comm%idx_recvt(i) = recvti(i)-comm%mystart_t+1
        end do
    end if
    deallocate(buff_sti)
    deallocate(recvti)

    return
  end subroutine swe_output_init

  subroutine swe_output_final(mesh)
    ! io
    type(global_domain), intent(inout), target :: mesh
    ! local
    type(group_comm),      pointer             :: comm
    type(block_structure), pointer             :: local_block

    local_block => mesh%local_block
    comm => mesh%gcomm_write

    deallocate(comm%idx_localv)
    deallocate(comm%idx_localt)
    deallocate(comm%idx_locale)
    deallocate(comm%counts_sendv)
    deallocate(comm%counts_sendt)
    deallocate(comm%counts_sende)
    deallocate(comm%displs_sendv)
    deallocate(comm%displs_sendt)
    deallocate(comm%displs_sende)
    deallocate(comm%counts_recvv)
    deallocate(comm%counts_recvt)
    deallocate(comm%counts_recve)
    deallocate(comm%displs_recvv)
    deallocate(comm%displs_recvt)
    deallocate(comm%displs_recve)
    deallocate(comm%idx_recvv)
    deallocate(comm%idx_recvt)
    deallocate(comm%idx_recve)

    return
  end subroutine swe_output_final

!================================================
! [1] Construct filename according to itimestep
! [2] Calculate Diag Vars using Prog vars
! [3] Do an Output Flow
!================================================

  subroutine swe_output_file(mesh,itimestep)
    use grist_util_module, only: write_string
! io
   type(global_domain), intent(in)    :: mesh
   integer(i4)         , intent(in)    :: itimestep
! local
   character(128)     :: c_glevel
   character(128)     :: c_testcase
   character(128)     :: c_ini
   character(128)     :: c_pv_order
   character(128)     :: c_adv
   character(len=5)   :: day
   character(len=5)   :: sec

!================================================
! [1]  Filename
!================================================

   call write_string(mesh%glevel , c_glevel)
   call write_string(testcase    , c_testcase)
   call write_string(initialfield, c_ini)
   call write_string(pv_order(3) , c_pv_order)
   call write_string(advection_scheme, c_adv)

   call get_current_date(itimestep,swe_timestep,day,sec)

   fname_output="GRIST.SWE.G"//trim(c_glevel)//".TC"//trim(c_testcase)//"."//&
                                         "INI"//trim(c_ini)//"."//&
                           trim(time_scheme)//".PV"//&
                           trim(c_pv_order)//".ADV"//&
                           trim(c_adv)//"."//&
                           trim(conserve_scheme)//"."//&
                           trim(day)//"-"//&
                           trim(sec)//".nc"

!================================================
! [2] Diagnosing Vars
!================================================

    if(itimestep.eq.0) then
       ue_reconst = ue_init
       ve_reconst = ve_init
       up_reconst = up_init
       vp_reconst = vp_init
    else
    ! reconstruct wind for diagnosis
       if(.not.isolate_advection_test) call project_uv(mesh,scalar_normal_velocity_at_edge%f,ue_reconst%f,ve_reconst%f)
    end if

    !if(.not.isolate_advection_test) call wind_edge_to_cell(mesh,ue_reconst%f,up_reconst%f)
    !if(.not.isolate_advection_test) call wind_edge_to_cell(mesh,ve_reconst%f,vp_reconst%f)
    if(.not.isolate_advection_test) call vector_recon_perot_edge2cell_uv(mesh,scalar_normal_velocity_at_edge%f,up_reconst%f,vp_reconst%f)

! div and 3 vor
    if(.not.isolate_advection_test) call diagnose_divergence_vorticity(mesh,scalar_normal_velocity_at_edge        ,&
                                            scalar_height_at_prime_cell           ,&
                                            scalar_divergence_at_prime_cell       ,&
                                            scalar_absolute_vorticity_at_dual_cell,&
                                            scalar_relative_vorticity_at_dual_cell,&
                                            scalar_potential_vorticity_at_dual_cell)
! start output
     call wrap_output_init_1d(mesh)
! add variables
     call wrap_add_field_1d(lon_nv                         ,"lon_nv")
     call wrap_add_field_1d(lat_nv                         ,"lat_nv")
     call wrap_add_field_1d(lon_ne                         ,"lon_ne")
     call wrap_add_field_1d(lat_ne                         ,"lat_ne")
     call wrap_add_field_1d(lon_nt                         ,"lon_nt")
     call wrap_add_field_1d(lat_nt                         ,"lat_nt")
     call wrap_add_field_1d(scalar_area_of_prime_cell      ,"PC_AREA")  ! constant
     call wrap_add_field_1d(scalar_area_of_dual_cell       ,"DC_AREA")  ! constant
     call wrap_add_field_1d(scalar_leng_of_hx_edge         ,"HX_EDGE")  ! constant
     call wrap_add_field_1d(scalar_leng_of_tr_edge         ,"TR_EDGE")  ! constant
     call wrap_add_field_1d(scalar_topo_at_prime_cell      ,"TOPO")
     call wrap_add_field_1d(scalar_height_at_prime_cell    ,"H")
     call wrap_add_field_1d(scalar_height_at_dual_cell     ,"H_DC")

     if(.not.isolate_advection_test)then
        call wrap_add_field_1d(scalar_hb_at_prime_cell                 ,"HB")
        call wrap_add_field_1d(up_reconst                              ,"UP")
        call wrap_add_field_1d(vp_reconst                              ,"VP")
        call wrap_add_field_1d(scalar_divergence_at_prime_cell         ,"DIV")
        call wrap_add_field_1d(ue_reconst                              ,"UE")
        call wrap_add_field_1d(ve_reconst                              ,"VE")
        call wrap_add_field_1d(scalar_absolute_vorticity_at_dual_cell  ,"AVOR")
        call wrap_add_field_1d(scalar_relative_vorticity_at_dual_cell  ,"RVOR")
        call wrap_add_field_1d(scalar_potential_vorticity_at_dual_cell ,"PVOR")
     end if
    
! outputing
    !call wrap_output_1d(mesh, itimestep, swe_timestep, outdir, fname_output)
    call wrap_output_1d_group(mesh, outdir, fname_output)
! clean
    call wrap_output_clean_1d()

    return
  end subroutine swe_output_file

  end module swe_inout_module
