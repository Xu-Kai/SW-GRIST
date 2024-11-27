module module_cu_ResNet_hybrid_v2
! Yiming Wang 2023 05 
! Reconstruct ResNet in Fortran
#ifdef Resnet
  use, intrinsic :: iso_fortran_env, only : sp => real32
!   use ftorch
  use grist_constants,          only: r4, i4, r8, one, zero, latvap, cp, one, half,gravity,fillvalue
  use grist_math_module,           only: lininterp
  use grist_hpe_constants,         only: eta_face_a, eta_face_b, eta_full_a, eta_full_b, p0
  use grist_nml_module,            only: nlev, nlevp, nlev_inidata, CU_ml_model
  use grist_datam_static_data_module, only: mldata_bias0, mldata_bias1, mldata_biases, mldata_imin, mldata_imax, &
                                            mldata_omin, mldata_omax, mldata_weight0, mldata_weight1_hybrid, mldata_weights, &
                                            mldata_plev, mldata_ps, mldata_qqq, mldata_inre, mldata_pred
  use grist_dycore_vars_module,     only: dycoreVarCellFace, dycoreVarCellFull, dycoreVarEdgeFull, dycoreVarEdgeFace, dycoreVarSurface
  use grist_tracer_transport_vars_module, only: tracerVarCellFull
!  use grist_wv_saturation,         only: findsp

  implicit none

!   private              ! Make default type private to the module
!   save

  integer,parameter:: nsize=3,in_channels=5,t0_channels=0
  integer,parameter:: nkernels=128,out_channels=2
  integer,parameter:: batch=3720,  nres=10,length=30,lengtho=30
  integer,parameter:: seq_len=1,begintime=4,casename=11
  integer,parameter:: tin_channels=in_channels*seq_len+t0_channels
  real*4,parameter:: ran=0.0d0, alpha=0.0d0
!   type(torch_module) :: CU_model
  integer, parameter :: wp = sp

!   public ResNet_hybrid_tend
!   public ResNet_hybrid_init


  contains
!====================================================================
subroutine ResNet_hybrid_init(rthcuten,rqvcuten,         &
                     restart,                  &
                     ids, ide, jds, jde, kds, kde,                      &
                     ims, ime, jms, jme, kms, kme,                      &
                     its, ite, jts, jte, kts, kte)
!--------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------
   integer , intent(in)           ::  ids, ide, jds, jde, kds, kde, &
                                      ims, ime, jms, jme, kms, kme, &
                                      its, ite, jts, jte, kts, kte

   real,     dimension( ims:ime , kms:kme , jms:jme ) , intent(out) ::  &
                                                              rthcuten, &
                                                              rqvcuten
   integer :: i, j, k, itf, jtf, ktf
   logical , intent(in)           ::  restart
   !jtf=min0(jte,jde-1)
   !ktf=min0(kte,kde-1)
   !itf=min0(ite,ide-1)
   jtf=jme
   ktf=kme
   itf=ime
   
!    CU_model = torch_module_load(CU_ml_model)
    i = ime - ims + 1
   call init_cu_weight(i)

   if(.not.restart)then
     do j=jts,jtf
     do k=kts,ktf
     do i=its,itf
       rthcuten(i,k,j)=0.
       rqvcuten(i,k,j)=0.
     enddo
     enddo
     enddo

     DO j=jts,jtf
     DO k=kts,ktf
     DO i=its,itf
!        rthften(i,k,j)=0.
!        rqvften(i,k,j)=0.
     ENDDO
     ENDDO
     ENDDO

   endif

end subroutine

subroutine ResNet_hybrid_tend(dt,itimestep,stepcu                            &
                ,u3d,v3d,w,th3d,qv3d,qdyn,tdyn,hfx,lh,coszr        &
                ,pcps,p8w,lat,raincv,cu_act_flag,dx             &
                ,ids,ide, jds,jde, kds,kde                      &
                ,ims,ime, jms,jme, kms,kme                      &
                ,its,ite, jts,jte, kts,kte                      &
                ,rthcuten,rqvcuten)
!inputs
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- th3d        3d potential temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- therad3d        3d heating rate due to radiation process ratio (kg/kg)
!-- p8w         3d hydrostatic pressure at full levels (pa)
!-- pcps        3d hydrostatic pressure at half levels (pa)
!s-1)
!-- rthcuten          theta tendency due to 
!                 cumulus scheme precipitation (k/s)
!-- rqvcuten          qv tendency due to 
!                 cumulus scheme precipitation (kg/kg/s)
!-- dz8w        dz between full levels (m)
!-- dt          time step (s)
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
      integer, intent(in) ::            ids,ide, jds,jde, kds,kde,      &
                                        ims,ime, jms,jme, kms,kme,      &
                                        its,ite, jts,jte, kts,kte,      &
                                        itimestep,                      &
                                        stepcu

      real,    intent(in) ::                                            &
                                        dt
      real,    dimension(ims:ime), intent(in) ::            dx


      real,    dimension(ims:ime, jms:jme), intent(out) ::               &
                                        raincv
      real,    dimension(ims:ime), intent(in) ::               &
                                        lat                                       
      real,    dimension(ims:ime, jms:jme), intent(in) ::               &
                                        hfx,                             &
                                        lh,                              &
                                        coszr

      logical, dimension(ims:ime,jms:jme), intent(inout) ::             &
                                        cu_act_flag

      real,    dimension(ims:ime, kms:kme, jms:jme), intent(in) ::      &
                                        qdyn,                           &
                                        pcps,                           &
                                        tdyn,                            &
                                        qv3d,                           &
                                        p8w,                          &
                                        th3d,                            &
                                        u3d,                            &
                                        v3d,                            &
                                        w
!      real,    dimension(kms:kme,ims:ime),           intent(in) ::      &
!                                        w

!--------------------------- optional vars ----------------------------
!--------------------------- optional vars ----------------------------

      real, dimension(ims:ime, kms:kme, jms:jme),                       &
               optional, intent(inout) ::                               &
                                        rqvcuten,                       &
                                        rthcuten
!--------------------------- local vars ------------------------------
      real      ::                                      &
                                        delt,                           &
                                        rdelt

      real     , dimension(its:ite) ::                  &
                                        rcs,                            &
                                        rn,                             &
                                        evap,                           &
                                        heatflux
      integer  , dimension(its:ite) ::  slimsk


      real     , dimension(its:ite, kts:kte+1) ::         &
                                        prsi,           &
                                        ghti,           &
                                          zl,           &
                                          zi
      integer, dimension(its:ite) ::                                    &
                                        kbot,                           &
                                        ktop

      integer ::                                                        &
                                        i,                              &
                                        im,                             &
                                        j,                              &
                                        k,                              &
                                        km,                             &
                                        kp,                             &
                                        kx,                             &
                                        kx1,                            &
                                        ilev
      real     , dimension(its:ite, 1:lengtho) ::  q1, q2
      real     , dimension(1:30) ::  q1w,q2w  

      integer                         :: nlev_eff,nlev_eff_t,nlev_eff_b
      real                            :: pface(29+1), hpface(29+1)
      real                            :: pfull(29),   hpfull(29)
      real                            :: dpfull(30),  dhpfull(29)
      real                            :: hpmfull(nlev),hpmface(nlevp)
      integer(i4)                         ::   qlev
!-------other local variables----
      integer                      :: zz
!-----------------------------------------------------------------------
!
!
!***  check to see if this is a convection timestep
!

!-----------------------------------------------------------------------
      do j=jts,jte
         do i=its,ite
           cu_act_flag(i,j)=.true.
         enddo
      enddo
      im=ite-its+1
      kx=kte-kts+1
      kx1=kx+1
      delt=dt*stepcu
      rdelt=1./delt
      raincv=0

!-------------  j loop (outer)
!--------------------------------------------------
   do j=jts,jte

      
!    call ResNet(q1(:,:),q2(:,:),u3d(:,:,j),v3d(:,:,j),w(:,:,j), &
!                th3d(:,:,j),qv3d(:,:,j),lh(:,1),hfx(:,1),coszr(:,1), &
!                pcps(:,:,j),qdyn(:,:,j),tdyn(:,:,j),kms,kme, its, ite)

!        print*, 'q1:',q1(i,1:19)
!        q1=q1-q2
        q1=q1/86400   !K/s
        q2=q2/86400*cp/latvap !kg/kg/s
    do i = its, ite
        if ( lat(i).lt. -70.0) then
            q1(i,:)=zero
            q2(i,:)=zero
!            q2=q2/10
        endif

        do ilev=11,30
           dpfull(ilev)=p8w(i,31-ilev,j)-p8w(i,32-ilev,j)
           raincv(i,1)=raincv(i,1)+q2(i, ilev)/gravity*dpfull(ilev)
        enddo


!        raincv(i,1)=raincv(i,1)/(abs(pfull(1)-pfull(nlev_eff)))
          
        raincv(i,1)=amax1(0.0,raincv(i,1)) !kg/m^2/s ->mm/s

!        rthcuten=0.
!        rqvcuten=0.
        rthcuten(i,1:18,1)=q1(i, 30:13:-1)
        rqvcuten(i,1:18,1)=-q2(i, 30:13:-1)
!        rthcuten(i,1:30,1)=q1(30:1:-1)
!        rqvcuten(i,1:30,1)=q2(30:1:-1)

      enddo
   enddo

end subroutine 

subroutine ResNet(q1, q2, u3d, v3d, w, &
                  th3d, qv3d, lh,hfx,coszr, &
                  pcps,qdyn,tdyn, kms,kme,its, ite)
    implicit none
    integer , intent(in)           ::   its, ite
    integer , intent(in)           ::   kms, kme


    real,    dimension(its:ite, kms:kme), intent(in) ::      &
                                        pcps,                            &
                                        qv3d,                           &
                                        qdyn,                          &
                                        th3d,                            &
                                        u3d,                            &
                                        v3d,                            &
                                        w,                              &
                                        tdyn

    real,     dimension(its:ite, 1:lengtho) , intent(out) ::  &
                                                              q1, &
                                                              q2
    real,      dimension(its:ite), intent(in)   :: &
                       lh, &
                       hfx, &
                       coszr
                                                             ! rthften,rqvften
    integer:: iname, ires
    integer:: iseq
    integer:: i, j, k
    integer,parameter :: ptop_res=30, ptop_reso=30

    real(r4)::  inputs( length ,  in_channels*seq_len+t0_channels)
    real(r4)::  outputs( lengtho,  out_channels)
    real(r4)::  predictions(lengtho , out_channels)
    real(r4)::  sinput(its:ite,length, in_channels*seq_len+t0_channels), stemp(length, nkernels )
    real(r4)::  soutput(its:ite,lengtho, out_channels), spred(lengtho, out_channels)
    real(r4)::  output_gap(nkernels),output_dense(lengtho*out_channels)
    real(r4)::  inputs_max( length, in_channels*seq_len+t0_channels)
    real(r4)::  inputs_min( length, in_channels*seq_len+t0_channels)
    real(r4)::  outputs_max( lengtho, out_channels), outputs_min(lengtho,out_channels)
    real(r4)       :: weights0(nsize,tin_channels,nkernels), bias0(nkernels)
    real(r4)       :: weightsa(nsize,nkernels,nkernels),biasa(nkernels),weightsb(nsize,nkernels,nkernels),biasb(nkernels)
    real(r4)       :: weights1(nsize,nkernels,out_channels), bias1(out_channels)
    real::  qs,es,ilev,qtemp,qfix,tfix
    integer :: status

    weights0=mldata_weight0%f
    bias0=mldata_bias0%f
    weights1=mldata_weight1_hybrid%f
    bias1=mldata_bias1%f

    inputs_max=mldata_imax%f
    inputs_min=mldata_imin%f
    outputs_min=mldata_omin%f
    outputs_max=mldata_omax%f
    sinput(:,1:ptop_res,1)=qv3d(:,ptop_res:1:-1)
    sinput(:,1:ptop_res,2)=th3d(:,ptop_res:1:-1)
!    sinput(1:ptop_res,3)=qdyn(ptop_res:1:-1)*86400/cp*latvap
!    sinput(1:ptop_res,4)=-tdyn(ptop_res:1:-1)*86400
    sinput(:,1:ptop_res,3)=u3d(:,ptop_res:1:-1)
    sinput(:,1:ptop_res,4)=v3d(:,ptop_res:1:-1)

!    sinput(1:ptop_res,5)=w(:)
!    sinput(1:ptop_res,5)=w(ptop_res:1:-1)
    sinput(:,1:ptop_res,5)=pcps(:,ptop_res:1:-1)
!    do i=1,ptop_res
!        sinput(i,5)=pcps(1)
!       sinput(i,8)=lh
!       sinput(i,9)=coszr
!    enddo

!    sinput=(sinput+1)*(inputs_max- inputs_min)/2 +inputs_min

    do i = its, ite
        sinput(i,:,:)=(sinput(i,:,:)-inputs_min(:,:))/(inputs_max(:,:)-inputs_min(:,:))*2.d0-1.d0
    end do

    if (.true.) then
    do i=1, in_channels*seq_len+t0_channels
    do j=1,length
    do k=its,ite
        if (sinput(k,j,i)>=1.d0+ran)then
            sinput(k,j,i)=1.d0+ran
        elseif (sinput(k,j,i)<=-1.d0-ran)then
            sinput(k,j,i)=-1.d0-ran
        endif
    enddo
    enddo
    enddo
    endif

    if (.true.)then
    do i=1,seq_len*in_channels+t0_channels
        do j=1,length
        do k=its,ite
        if (inputs_max(j,i)<=1e-29)then
            sinput(k,j,i)=0.d0
        endif
        enddo
        enddo
    enddo
    endif


   !!!!!
   !forward
   !!!!!
   call cu_dnn_forward(sinput, soutput, its, ite)


        ! ResNet complete
    if (.true.)then
    do i=1, out_channels
        do j=1,lengtho
        do k = its,ite
        if (soutput(k,j,i)>=1.d0+ran)then
            soutput(k,j,i)=1.d0+ran
        elseif (soutput(k,j,i)<=-1.d0-ran)then
            soutput(k,j,i)=-1.d0-ran
        endif
        enddo
        enddo
    enddo
    endif

    do i = its, ite
        soutput(i,:,:)=(soutput(i,:,:)+1.d0)*(outputs_max(:,:)-outputs_min(:,:))/2.d0+outputs_min(:,:)
    end do

    if (.true.)then
    do i=1, out_channels
        do j=1,lengtho
        do k=its,ite
        if (abs(soutput(k,j,i))<=1e-29)then
            soutput(k,j,i)=0.d0
        endif
        enddo
        enddo
    enddo
    endif


    q1=0.d0
    q2=0.d0

    q1(:,1:ptop_reso)=soutput(:,1:ptop_reso,1)!*3-sinput(1:ptop_res,4)*2
    q2(:,1:ptop_reso)=soutput(:,1:ptop_reso,2)!*3-sinput(1:ptop_res,3)*2
!    do ilev=1,ptop_reso
!       call findsp(qv3d(ptop_reso-ilev+1),th3d(ptop_reso-ilev+1),pcps(ptop_reso-ilev+1),.false.,es,qs,status)
!       if (qs .eq. fillvalue ) then
!          return
!       end if
!       print*,'qs:',qs
!       qtemp=qv3d(ptop_reso-ilev+1)-q2(ilev)/86400*cp/latvap
!       if (qtemp>qs) then
!           qfix=qtemp-qs
!           print*,'qfix:',qfix
!           tfix=-qfix*latvap/cp*86400
!           print*,'tfix:',tfix
!           q1(ilev)=q1(ilev)-tfix
!           q2(ilev)=q2(ilev)+qfix*latvap/cp*86400
!       endif
!    enddo



!    q2(1:ptop_reso)= sinput(1:ptop_res,3)
!    q1(1:ptop_reso)= sinput(1:ptop_res,4)


end subroutine

subroutine cu_dnn_forward(input, output, ims, ime)
    integer, intent(in) :: ims, ime

    real(r4), dimension(ims:ime, length, in_channels*seq_len+t0_channels), intent(in) :: input
    real(r4), dimension(ims:ime,lengtho, out_channels), intent(out) :: output

    ! real(wp), dimension(ims:ime, length, in_channels*seq_len+t0_channels), target :: in_data
    ! real(wp), dimension(ims:ime,lengtho, out_channels),  target :: out_data

    real(wp), allocatable :: in_data(:,:,:)
    real(wp), allocatable :: out_data(:,:,:)
    integer, parameter :: in_dims = 3
    integer, parameter :: out_dims = 3
    integer, parameter :: n_inputs = 1

    integer :: in_layout(in_dims) = [1,2,3]
    integer :: out_layout(out_dims) = [1,2,3]    

    ! type(torch_tensor), dimension(1) :: model_input
    ! type(torch_tensor)  :: model_output
    integer :: i, j, k
    ! in_data(:,:,:) = input(:,:,:)

    allocate(in_data( length, in_channels*seq_len+t0_channels, ims:ime))
    allocate(out_data(lengtho, out_channels, ims:ime))

    ! print *, "can shu", length, in_channels*seq_len+t0_channels, lengtho, out_channels
    do i = ims, ime
        do k = 1, length
             do j = 1, in_channels*seq_len+t0_channels
                in_data(k, j, i) =  input(i, k, j)
            enddo
        enddo 
    enddo 
    ! print *, "in is ", in_data(5, 1, ims)

    ! model_input(1)  = torch_tensor_from_array(in_data, in_layout, torch_kCPU)
    ! model_output = torch_tensor_from_array(out_data, out_layout, torch_kCPU)
      ! Infer
    i = ime - ims + 1
    call  cu_resnet_forward(in_data(1, 1, ims), out_data(1, 1, ims), i)
    ! call torch_module_forward(CU_model, model_input, n_inputs, model_output)
    
    ! output(:,:,:) = out_data(:,:,:)
 
    do i = ims, ime
        do k = 1, lengtho
             do j = 1, out_channels
                output(i, k, j) =  out_data(k, j, i)
            enddo
        enddo 
    enddo 
    ! print *, "out is ", out_data(5, 1, ims)

    deallocate(in_data)
    deallocate(out_data)
    ! call torch_tensor_delete(model_input(1))
    ! call torch_tensor_delete(model_output)

end subroutine cu_dnn_forward

#endif
end module module_cu_ResNet_hybrid_v2

