module module_ra_rrtmg_dnn
! Chenyu Wang 2024 02 
! MCT Fotran-python
#ifdef RRTMG_DNN
  use grist_nml_module,     only: RA_ml_model
  use, intrinsic :: iso_fortran_env, only : sp => real32
!   use ftorch
  use grist_mpi
  

  implicit none

  ! Generate an object to hold the Torch model
  ! Set precision for reals
!   type(torch_module) :: RA_model
  integer, parameter :: wp = sp

    ! Set up types of input and output data


  contains

!====================================================================
subroutine rrtmg_dnn_lw_init(rthraten,rthratenlw,cldfra,restart,          &
                          ids, ide, jds, jde, kds, kde,                &
                          ims, ime, jms, jme, kms, kme,                &
                          its, ite, jts, jte, kts, kte                 )
!--------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------
   logical , intent(in)           :: restart
   integer , intent(in)           :: ids, ide, jds, jde, kds, kde,  &
                                     ims, ime, jms, jme, kms, kme,  &
                                     its, ite, jts, jte, kts, kte

   real , dimension( ims:ime , kms:kme , jms:jme ) , intent(inout) ::        &
                                                          rthraten, &
                                                        rthratenlw, &
                                                            cldfra
   real :: pi
   integer :: grank, gsize
   integer :: ierr
   integer :: i, j, k, itf, jtf, ktf
  
   ! Initialise the Torch model to be used
   
   ! RA_model = torch_module_load(RA_ml_model)
   ! print *, "batch size is ", 
   i = ime - ims + 1
    call init_ra_weight(i)

    ! deallocate(start)


   !jtf=min0(jte,jde-1)
   !ktf=min0(kte,kde-1)
   !itf=min0(ite,ide-1)

   !if(.not.restart)then
   !  do j=jts,jtf
   !  do k=kts,ktf
   !  do i=its,itf
        rthraten  = 0.
        rthratenlw= 0.
        cldfra    = 0.
   !  enddo
   !  enddo
   !  enddo
   !endif


end subroutine rrtmg_dnn_lw_init

!====================================================================
subroutine rrtmg_dnn_lwrad(p3d,t3d,qv3d                                   &
                ,tsk,emiss,xcoszen,lat,coszr                        &
                ,ids,ide, jds,jde, kds,kde                      &
                ,ims,ime, jms,jme, kms,kme                      &
                ,its,ite, jts,jte, kts,kte                      &
                ,glw, gsw, swdnb)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------
!--------------------------------------------------------------------
    implicit none
!--------------------------------------------------------------------
    integer, intent(in) ::          ids,ide, jds,jde, kds,kde,      &
                                    ims,ime, jms,jme, kms,kme,      &
                                    its,ite, jts,jte, kts,kte      
    real,    dimension(ims:ime, kms:kme, jms:jme), intent(in) ::      &
                                    p3d,t3d,qv3d
    real,    dimension(ims:ime, jms:jme),  intent(in) ::      &
                                    tsk,emiss,xcoszen                    
    real, dimension(ims:ime, jms:jme),                       &
            optional, intent(inout) ::                               &
                                    glw, gsw, swdnb,coszr 
    real,    dimension(ims:ime), intent(in) ::               &
                                        lat  
    real, parameter ::  mpressureFace_max = 103952.4140625,       &
                        mpressureFace_min = 364.34646606,         &
                        temperature_max   = 323.11566162,         &
                        temperature_min   = 164.35217285,         &
                        tracerMxrt_max    = 0.02575767,           &
                        tracerMxrt_min    = 0,                    &
                        tskin_max         = 341.52825928,         &
                        tskin_min         = 187.74468994,         &
                        emiss_max         = 0.9939084,            &
                        emiss_min         = 0.96999866,           &
                        coszr_max         = 0.99770105,           &
                        coszr_min         = -0.99770141,          &
                        fswdsDiag_max     = 1272.85131836,        &
                        fswdsDiag_min     = 0,                    &
                        glw_max           = 525.09002686,         &
                        glw_min           = 47.51123428,          &
                        gsw_max           = 1119.45471191,        &
                        gsw_min           = 0
        
      ! Set up number of dimensions of input tensor and axis order
    integer, parameter :: ptop_res=30
    real, dimension(ims:ime,92)  :: in_data
    real, dimension(ims:ime,3)  :: out_data
    real, dimension(ims:ime, jms:jme) :: glw_old, gsw_old, swdnb_old

    ! Set up the model input and output as Fortran arrays
    integer :: in_shape, out_shape
    integer :: grank
    integer :: i
    real(8)              :: step_beg, step_end, step_elapse

    ! Allocate one-dimensional input/output arrays, based on multiplication of all input/output dimension sizes
    coszr = xcoszen  
    glw_old = glw
    
   in_data(:,1:ptop_res) = (p3d(:, ptop_res:1:-1, 1) - mpressureFace_min)/ (mpressureFace_max - mpressureFace_min)*2 -1
   in_data(:,ptop_res+1:ptop_res*2) = (t3d(:, ptop_res:1:-1, 1) - temperature_min)/(temperature_max - temperature_min)*2 -1
   in_data(:,ptop_res*2+1:ptop_res*3) = (qv3d(:, ptop_res:1:-1, 1) - tracerMxrt_min) /(tracerMxrt_max - tracerMxrt_min)*2 -1
   in_data(:,ptop_res*3+1) = (tsk(:, 1) - tskin_min)/(tskin_max - tskin_min)*2 -1
   ! in_data(:,ptop_res*3+2) = (emiss(:, 1) - emiss_min)/(emiss_max - emiss_min)*2 -1
   in_data(:,ptop_res*3+2) = (coszr(:, 1) - coszr_min)/(coszr_max - coszr_min)*2 -1
   
   call rrtmg_dnn_forward(in_data, out_data, ims, ime)

   glw(:, 1) = (out_data(:, 1) + 1 )/2*(glw_max - glw_min) + glw_min 
   gsw(:, 1) = (out_data(:, 2) + 1 )/2*(gsw_max - gsw_min) + gsw_min
   swdnb(:, 1) = (out_data(:, 3) + 1 )/2*(fswdsDiag_max - fswdsDiag_min) + fswdsDiag_min

   
   
   do i = its, ime
      if ((lat(i).gt. 75.0).or.(lat(i).lt. -75.0)) then
         gsw(i, 1) = 0.0
         swdnb(i, 1) = 0.0
      endif 
   enddo

   where (glw < 0.0)
      glw = 0.0
   end where
   where (gsw < 0.0)
      gsw = 0.0
   end where
   where (swdnb < 0.0)
      swdnb = 0.0
   end where


   
end subroutine rrtmg_dnn_lwrad

subroutine rrtmg_dnn_forward(input, output, ims, ime)
    integer, intent(in) :: ims, ime

    real, dimension(ims:ime,92), intent(in) :: input
    real, dimension(ims:ime,3), intent(out) :: output

    ! real(wp), dimension(ims:ime,92), target :: in_data
    ! real(wp), dimension(ims:ime,3),  target :: out_data

    integer, parameter :: in_dims = 2
    integer, parameter :: out_dims = 2
    integer, parameter :: n_inputs = 1

    integer :: in_layout(in_dims) = [1,2]
    integer :: out_layout(out_dims) = [1, 2]    
    real(wp), allocatable :: in_data(:,:)
    real(wp), allocatable :: out_data(:,:)
    integer :: i, j
   !  type(torch_tensor), dimension(1) :: model_input
   !  type(torch_tensor)  :: model_output

   !  in_data(:,:) = input(:,:)

   !  model_input(1)  = torch_tensor_from_array(in_data, in_layout, torch_kCPU)
   !  model_output = torch_tensor_from_array(out_data, out_layout, torch_kCPU)
   !    ! Infer
   !  call torch_module_forward(RA_model, model_input, n_inputs, model_output)
    
   !  output(:,:) = out_data(:,:)
 
   !  call torch_tensor_delete(model_input(1))
   !  call torch_tensor_delete(model_output)
   allocate(in_data(92, ims:ime))
   allocate(out_data(3, ims:ime))
   do i = ims, ime
      do j = 1, 92
         in_data(j, i) = input(i, j)
      enddo 
   enddo
    i = ime - ims +1
   !  print *, "in data is ", in_data(4, 4), input(1,3)

   call ra_forward(in_data(1, ims), out_data(1, ims), i)
   do i = ims, ime
      do j = 1, 3
         output(i, j) = out_data(j,i)
      enddo 
   enddo

   deallocate(in_data)
   deallocate(out_data)
end subroutine rrtmg_dnn_forward


!====================================================================
subroutine rrtmg_dnn_sw_init(rthraten,rthratenlw,cldfra,restart,          &
                          ids, ide, jds, jde, kds, kde,                &
                          ims, ime, jms, jme, kms, kme,                &
                          its, ite, jts, jte, kts, kte                 )
!--------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------
   logical , intent(in)           :: restart
   integer , intent(in)           :: ids, ide, jds, jde, kds, kde,  &
                                     ims, ime, jms, jme, kms, kme,  &
                                     its, ite, jts, jte, kts, kte

   real , dimension( ims:ime , kms:kme , jms:jme ) , intent(inout) ::        &
                                                          rthraten, &
                                                        rthratenlw, &
                                                            cldfra
   real :: pi

   integer :: i, j, k, itf, jtf, ktf

   !jtf=min0(jte,jde-1)
   !ktf=min0(kte,kde-1)
   !itf=min0(ite,ide-1)

   !if(.not.restart)then
   !  do j=jts,jtf
   !  do k=kts,ktf
   !  do i=its,itf
        rthraten  = 0.
        rthratenlw= 0.
        cldfra    = 0.
   !  enddo
   !  enddo
   !  enddo
   !endif


end subroutine rrtmg_dnn_sw_init

!====================================================================
subroutine rrtmg_dnn_swrad(&
                ids,ide, jds,jde, kds,kde                      &
                ,ims,ime, jms,jme, kms,kme                      &
                ,its,ite, jts,jte, kts,kte)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------
!--------------------------------------------------------------------
    implicit none
!--------------------------------------------------------------------
    integer, intent(in) ::          ids,ide, jds,jde, kds,kde,      &
                                    ims,ime, jms,jme, kms,kme,      &
                                    its,ite, jts,jte, kts,kte                         

end subroutine rrtmg_dnn_swrad

subroutine rrtmg_dnn_lw_final()
   ! call torch_module_delete(RA_model)

end subroutine rrtmg_dnn_lw_final

#endif

end module module_ra_rrtmg_dnn
