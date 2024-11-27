 !======================================================
 !         NAMELIST module for WRF-based physics
 !======================================================

 module grist_wrfphys_nml_module

   use grist_constants,   only: i4, r8
   use grist_mpi
   use grist_nml_module,  only: model_timestep

   implicit none
   
   private

   public  :: set_wrfphys_nml,   &
!
! WRF physics
!
              wrfphys_cu_scheme,  &
              wrfphys_mp_scheme,  &
              wrfphys_bl_scheme,  &
              wrfphys_cf_scheme,  &
              wrfphys_ra_scheme,  &
              wrfphys_rasw_scheme,&
              wrfphys_ralw_scheme,&
              wrfphys_sf_scheme,  &
              wrfphys_lm_scheme,  &
              step_cu,            &
              step_ra,            &
              step_bl,            &
              wphys_has_req,      &
              use_gwdo,           &
              use_cond,           &
              unuse_cu,           &
              force_read_thompson,&
              write_thompson_tables
!
! WRF physics package related
!
   character(100)     :: wrfphys_cu_scheme
   character(100)     :: wrfphys_mp_scheme
   character(100)     :: wrfphys_bl_scheme
   character(100)     :: wrfphys_cf_scheme
   character(100)     :: wrfphys_ra_scheme
   character(100)     :: wrfphys_rasw_scheme
   character(100)     :: wrfphys_ralw_scheme
   character(100)     :: wrfphys_sf_scheme
   character(100)     :: wrfphys_lm_scheme
!
! these steps control those slow physics (PDC using tendency method)
! Fast physics step will be controlled by model_timestep (PDC using operator splitting)
!
   integer(i4)        :: step_cu  = 1  ! how many model_timestep to call cumulus, default every step
   integer(i4)        :: step_ra   ! how many model_timestep to call radiation, default: 3h
   integer(i4)        :: step_bl  = 1  ! not used as bl is currently a fast physics
   integer(i4)        :: wphys_has_req  = 0
   logical            :: use_gwdo
   logical            :: use_cond
   logical            :: unuse_cu = .false.
   logical            :: force_read_thompson = .false.
   logical            :: write_thompson_tables = .false.

  contains

  subroutine set_wrfphys_nml()

!================================================
! global vars have been defined in the header
!================================================

! local
  character(len=300) :: filename
  integer (i4)       :: fileunit

   namelist /wrfphys_para/wrfphys_cu_scheme,  & ! CUmulus
                          wrfphys_mp_scheme,  & ! MicroPhysics
                          wrfphys_bl_scheme,  & ! Boundary Layer
                          wrfphys_cf_scheme,  & ! RAdiation
                          wrfphys_ra_scheme,  & ! RAdiation
                          wrfphys_rasw_scheme,& ! RAdiation
                          wrfphys_ralw_scheme,& ! RAdiation
                          wrfphys_sf_scheme,  & ! Surface Flux
                          wrfphys_lm_scheme,  & ! Land Model
                          step_cu,            &
                          step_ra,            &
                          step_bl,            &
                          wphys_has_req,      &
                          use_gwdo,           &
                          use_cond,           &
                          unuse_cu,           &
                          force_read_thompson,&
                          write_thompson_tables

    filename = "grist_amipw_phys.nml"

    fileunit = 1

    open  (fileunit, status='old',file=filename)
    read  (fileunit, nml=wrfphys_para)
    close (fileunit)

  !  step_ra = 10800/idint(model_timestep)

    if(mpi_rank() .eq. 0)then
       print*,"**********************************************************"
       print*,"     The WRFPhys is used following: ", trim(filename)
       print*,"     step_ra, step_cu, step_bl are:", step_ra, step_cu, step_bl
       print*,"**********************************************************"
    end if

    return
  end subroutine set_wrfphys_nml

  end module grist_wrfphys_nml_module
