
!----------------------------------------------------------------------------
! Created on 2016
! Version 1.0
! Description: Restart module
! Revision history: 
!----------------------------------------------------------------------------

 module swe_restart_module

  use grist_domain_types, only: global_domain

  use grist_constants,    only: i4

  use grist_nml_module,   only: fname_restart,&
                                time_scheme,  &
                                swe_timestep, &
                                testcase,     &
                                outdir,       &
                                conserve_scheme,  &
                                pv_order,       &
                                advection_scheme

  use grist_fileio_list_1d_module, only: wrap_read_1d,       &
                                   wrap_output_init_1d,&
                                   wrap_add_field_1d,  &
                                   wrap_output_1d,     &
                                   wrap_output_clean_1d

  use grist_time_manager,  only: get_current_date

  use swe_vars_module,    only: tend_normal_velocity_at_edge_nm2   ,&
                                tend_normal_velocity_at_edge_nm1   ,&
                                tend_normal_velocity_at_edge       ,&
                                scalar_normal_velocity_at_edge_nm2 ,&
                                scalar_normal_velocity_at_edge_nm1 ,&
                                scalar_normal_velocity_at_edge     ,&
                                tend_height_at_prime_cell_nm2      ,&   
                                tend_height_at_prime_cell_nm1      ,&
                                tend_height_at_prime_cell          ,&
                                scalar_height_at_prime_cell_nm2    ,&
                                scalar_height_at_prime_cell_nm1    ,&
                                scalar_height_at_prime_cell       
  use grist_util_module, only: write_string
                                          

  implicit none

  private

  public   :: swe_read_restart_file    ,&
              swe_write_restart_file   ,&
              restart_timestep

     integer(i4) :: restart_timestep

  contains

!================================================
!  If Called, read from restart file. 
!  Overwrite data in swe_vars_module
!  Although we write so many data into restart,
!  not all of them are used.
!================================================

  subroutine swe_read_restart_file(mesh)
! io
   type(global_domain), intent(in)    :: mesh
! local
   character(128)     :: c_glevel
   character(128)     :: c_testcase
   character(128)     :: c_adv
   character(128)     :: c_pv
   character(len=5)   :: day
   character(len=5)   :: sec

!================================================
! Construct Filename
!================================================
! read the timestep at which we generate this restart file

   open(1,file=trim(outdir)//'time_restart.txt',status='old')
   read(1,*) restart_timestep
   close(1)

   call write_string(mesh%glevel     , c_glevel)
   call write_string(testcase        , c_testcase)
   call write_string(pv_order(3)     , c_pv)
   call write_string(advection_scheme, c_adv)

   call get_current_date(restart_timestep,swe_timestep,day,sec)

   fname_restart="GRIST.RST.G"//trim(c_glevel)//".TC"//trim(c_testcase)//"."//&
                               trim(time_scheme)//".ADV"//&
                               trim(c_adv)//".PV"//&
                               trim(c_pv)//"."//&
                               trim(conserve_scheme)//"."//&
                               trim(day)//"-"//&
                               trim(sec)//".nc"

   call wrap_read_1d(outdir,fname_restart,'tend_normal_velocity_at_edge_nm2'  ,mesh%ne,tend_normal_velocity_at_edge_nm2)
   call wrap_read_1d(outdir,fname_restart,'tend_normal_velocity_at_edge_nm1'  ,mesh%ne,tend_normal_velocity_at_edge_nm1)
   call wrap_read_1d(outdir,fname_restart,'tend_normal_velocity_at_edge'      ,mesh%ne,tend_normal_velocity_at_edge)
   call wrap_read_1d(outdir,fname_restart,'scalar_normal_velocity_at_edge_nm2',mesh%ne,scalar_normal_velocity_at_edge_nm2)
   call wrap_read_1d(outdir,fname_restart,'scalar_normal_velocity_at_edge_nm1',mesh%ne,scalar_normal_velocity_at_edge_nm1)
   call wrap_read_1d(outdir,fname_restart,'scalar_normal_velocity_at_edge'    ,mesh%ne,scalar_normal_velocity_at_edge)

   call wrap_read_1d(outdir,fname_restart,'tend_height_at_prime_cell_nm2'  ,mesh%nv,tend_height_at_prime_cell_nm2)
   call wrap_read_1d(outdir,fname_restart,'tend_height_at_prime_cell_nm1'  ,mesh%nv,tend_height_at_prime_cell_nm1)
   call wrap_read_1d(outdir,fname_restart,'tend_height_at_prime_cell'      ,mesh%nv,tend_height_at_prime_cell)
   call wrap_read_1d(outdir,fname_restart,'scalar_height_at_prime_cell_nm2',mesh%nv,scalar_height_at_prime_cell_nm2)
   call wrap_read_1d(outdir,fname_restart,'scalar_height_at_prime_cell_nm1',mesh%nv,scalar_height_at_prime_cell_nm1)
   call wrap_read_1d(outdir,fname_restart,'scalar_height_at_prime_cell'    ,mesh%nv,scalar_height_at_prime_cell)

   return

  end subroutine swe_read_restart_file

!================================================
! [1] Construct Filename
! [2] Write Restart File
!================================================

  subroutine swe_write_restart_file(mesh,itimestep)
! io
   type(global_domain), intent(in)    :: mesh
   integer(i4)         , intent(in)    :: itimestep
! local
   character(100)     :: c_glevel
   character(100)     :: c_testcase
   character(100)     :: c_adv
   character(100)     :: c_pv
   character(len=5)   :: day
   character(len=5)   :: sec

!================================================
! [1] Construct Filename
!================================================

   call write_string(mesh%glevel     , c_glevel)
   call write_string(testcase        , c_testcase)
   call write_string(pv_order(3)     , c_pv)
   call write_string(advection_scheme, c_adv)

   call get_current_date(itimestep,swe_timestep,day,sec)

   fname_restart="GRIST.RST.G"//trim(c_glevel)//".TC"//trim(c_testcase)//"."//&
                               trim(time_scheme)//".ADV"//&
                               trim(c_adv)//".PV"//&
                               trim(c_pv)//"."//&
                               trim(conserve_scheme)//"."//&
                               trim(day)//"-"//&
                               trim(sec)//".nc"

!================================================
! [2] Write Restart File
!================================================

    call wrap_output_init_1d(mesh)
! normal wind
    call wrap_add_field_1d(tend_normal_velocity_at_edge_nm2  , "tend_normal_velocity_at_edge_nm2")
    call wrap_add_field_1d(tend_normal_velocity_at_edge_nm1  , "tend_normal_velocity_at_edge_nm1")
    call wrap_add_field_1d(tend_normal_velocity_at_edge      , "tend_normal_velocity_at_edge")
    call wrap_add_field_1d(scalar_normal_velocity_at_edge_nm2, "scalar_normal_velocity_at_edge_nm2")
    call wrap_add_field_1d(scalar_normal_velocity_at_edge_nm1, "scalar_normal_velocity_at_edge_nm1")
    call wrap_add_field_1d(scalar_normal_velocity_at_edge    , "scalar_normal_velocity_at_edge")
! height
    call wrap_add_field_1d(tend_height_at_prime_cell_nm2     , "tend_height_at_prime_cell_nm2")
    call wrap_add_field_1d(tend_height_at_prime_cell_nm1     , "tend_height_at_prime_cell_nm1")
    call wrap_add_field_1d(tend_height_at_prime_cell         , "tend_height_at_prime_cell")
    call wrap_add_field_1d(scalar_height_at_prime_cell_nm2   , "scalar_height_at_prime_cell_nm2")
    call wrap_add_field_1d(scalar_height_at_prime_cell_nm1   , "scalar_height_at_prime_cell_nm1")
    call wrap_add_field_1d(scalar_height_at_prime_cell       , "scalar_height_at_prime_cell")

    call wrap_output_1d(mesh,itimestep,swe_timestep,outdir,fname_restart)
    call wrap_output_clean_1d

! write the timestep at which we generate last restart file
    open(1,file=trim(outdir)//'time_restart.txt',status='unknown')
    write(1,*) itimestep
    close(1)

    return

   end subroutine swe_write_restart_file

 end module swe_restart_module
