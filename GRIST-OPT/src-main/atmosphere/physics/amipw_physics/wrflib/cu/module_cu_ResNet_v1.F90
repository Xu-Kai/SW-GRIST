module module_cu_ResNet_v1
! Yiming Wang 2023 05 
! Reconstruct ResNet in Fortran
#ifdef Resnet

  use grist_constants,          only: r4, i4, r8, one, zero, latvap, cp, one, half,gravity
  use grist_math_module,           only: lininterp
  use grist_hpe_constants,         only: eta_face_a, eta_face_b, eta_full_a, eta_full_b, p0
  use grist_nml_module,            only: nlev, nlevp, nlev_inidata
  use grist_datam_static_data_module, only: mldata_bias0, mldata_bias1, mldata_biases, mldata_imin, mldata_imax, &
                                            mldata_omin, mldata_omax, mldata_weight0, mldata_weight1, mldata_weights, &
                                            mldata_plev, mldata_ps, mldata_qqq, mldata_inre, mldata_pred
  use grist_dycore_vars_module,     only: dycoreVarCellFace, dycoreVarCellFull, dycoreVarEdgeFull, dycoreVarEdgeFace, dycoreVarSurface
  use grist_tracer_transport_vars_module, only: tracerVarCellFull

  implicit none

  private              ! Make default type private to the module
  save

  integer,parameter:: nsize=3,in_channels=6,t0_channels=0
  integer,parameter:: nkernels=128,out_channels=2
  integer,parameter:: batch=3720,  nres=10,length=30,lengtho=29
  integer,parameter:: seq_len=1,begintime=4,casename=11
  integer,parameter:: tin_channels=in_channels*seq_len+t0_channels
  real*4,parameter:: ran=0.0d0, alpha=0.0d0

  public ResNet_tend
  public ResNet_init


  contains
!====================================================================
subroutine ResNet_init(rthcuten,rqvcuten,         &
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

subroutine ResNet_tend(dt,itimestep,stepcu                            &
                ,u3d,v3d,w,th3d,qv3d,therad3d        &
                ,dz8w,pcps,p8w,raincv,cu_act_flag,dx             &
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

      logical, dimension(ims:ime,jms:jme), intent(inout) ::             &
                                        cu_act_flag

      real,    dimension(ims:ime, kms:kme, jms:jme), intent(in) ::      &
                                        dz8w,                           &
                                        pcps,                           &
                                        p8w,                            &
                                        qv3d,                           &
!                                        therad3d,                          &
                                        th3d,                            &
                                        u3d,                            &
                                        v3d,                            &
                                        w
!      real,    dimension(kms:kme,ims:ime),           intent(in) ::      &
!                                        w
      real,    dimension(kms:kme,ims:ime), intent(in) ::                &
                                        therad3d

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
      real     , dimension(1:lengtho) ::  q1, q2
      real     , dimension(1:30) ::  q1w,q2w  

      integer                         :: nlev_eff,nlev_eff_t,nlev_eff_b
      real                            :: pface(29+1), hpface(29+1)
      real                            :: pfull(29),   hpfull(29)
      real                            :: dpfull(29),  dhpfull(29)
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

      do i=its,ite 
        call ResNet(q1(:),q2(:),u3d(i,:,j),v3d(i,:,j),w(i,:,j), &
                    th3d(i,:,j),qv3d(i,:,j),therad3d(:,i), &
                    pcps(i,:,j),kms,kme)
!        print*, 'q1:',q1(i,1:19)
        q1=q1/86400   !K/s
        q2=q2/86400*cp/latvap !kg/kg/s
       
!
! nominal pressure at ERA pressure levels
!
        pfull(1:29)=mldata_plev%f
!        print*, 'pfull :', pfull
!
! check effective full level number
!
        nlev_eff = 99999
        do ilev = 1, 29
           if(pfull(ilev).ge.mldata_ps%f(i))then
              nlev_eff = ilev-1
              exit
           end if
        end do
        if(nlev_eff.eq.99999) nlev_eff = 29

        pface(2:nlev_eff)       = (pfull(1:nlev_eff-1)+pfull(2:nlev_eff))*half
        pface(1)                = zero
        pface(nlev_eff+1)       = mldata_ps%f(i) ! if ps<pface(nlev_eff), let it be, good assumption! no problem till now...
        dpfull(1:nlev_eff)      = pface(2:nlev_eff+1)-pface(1:nlev_eff)
        hpface(1)               = pface(1)
!        print*,'hpface:',hpface
!        print*,'dpfull:',dpfull
!        print*,'qqq:',mldata_qqq%f(:,i)
!
! evaluate dry mass at ERA pressure level and surface
!
        do ilev = 2, nlev_eff+1
           hpface(ilev) = hpface(ilev-1)+dpfull(ilev-1)*(one-mldata_qqq%f(ilev-1,i))
        end do
!        print*,'qqq:',mldata_qqq%f(:,i)
        hpfull(1:nlev_eff) = half*(hpface(1:nlev_eff)+hpface(2:nlev_eff+1))

        hpmface(1:nlevp)  = eta_face_a(1:nlevp)*p0+eta_face_b(1:nlevp)*hpface(nlev_eff+1) ! Pa
        hpmfull(1:nlev)   = eta_full_a(1:nlev) *p0+eta_full_b(1:nlev) *hpface(nlev_eff+1)
!        print*,'nlev_eff',nlev_eff
!        print*,'q1w:', q1
!        print*,'q2w:', q2
!        print*,'dpfull:',dpfull
!        print*,'hpmfull:',hpmfull
!        print*,'rthcuten:',rthcuten(i,:,1)
!        print*,'rqvcuten:',rqvcuten(i,:,1)
!        print*,'ii',i
        if(nlev_eff.ge.29) nlev_eff = 29

        do ilev=10,nlev_eff
           raincv(i,1)=raincv(i,1)+q2(ilev)/gravity*dpfull(ilev)
        enddo

!        raincv(i,1)=raincv(i,1)/(abs(pfull(1)-pfull(nlev_eff)))
          
        raincv(i,1)=amax1(0.0,raincv(i,1)) !kg/m^2/s ->mm/s

!        rthcuten=0.
!        rqvcuten=0.

        nlev_eff_b = 99999
        do ilev = 1, 30
           if(dycoreVarCellFull%scalar_hpressure_n%f(ilev,i).ge.hpfull(nlev_eff))then
              nlev_eff_b = ilev
              exit
           end if
        end do
        if(nlev_eff_b.eq.99999) nlev_eff_b = 30

        nlev_eff_t = 99999
        do ilev = 1, 30
           if(dycoreVarCellFull%scalar_hpressure_n%f(ilev,i).ge.hpfull(1))then
              nlev_eff_t = ilev
              exit
           end if
        end do
        if(nlev_eff_t.eq.99999) nlev_eff_t = 1



!
! interpolate from ERA pressure level to GRIST model level based on dry air mass
!
!        call lininterp (q1w(1:19,i), hpfull(1:19), 1, 19, &
!                        rthcuten(i,6:30,1), hpmfull(6:30), nlev-5)
!        call lininterp (q2w(1:19,i), hpfull(1:19), 1, 19, &
!                        rqvcuten(i,6:30,1), hpmfull(6:30), nlev-5)
        call lininterp (q1(1:nlev_eff), hpfull(1:nlev_eff), 1, nlev_eff, &
                        q1w(nlev_eff_t:nlev_eff_b), dycoreVarCellFull%scalar_hpressure_n%f(nlev_eff_t:nlev_eff_b,i), nlev_eff_b-nlev_eff_t+1)
        call lininterp (q2(1:nlev_eff), hpfull(1:nlev_eff), 1, nlev_eff, &
                        q2w(nlev_eff_t:nlev_eff_b), dycoreVarCellFull%scalar_hpressure_n%f(nlev_eff_t:nlev_eff_b,i), nlev_eff_b-nlev_eff_t+1)

!        q1w(1:nlev_eff_t)=0
!        q1w(nlev_eff_b:30)=0
!        q2w(1:nlev_eff_t)=0
!        q2w(nlev_eff_b:30)=0a
        q2w(1:10)=0

        rthcuten(i,1:30,1)=q1w(30:1:-1)
        rqvcuten(i,1:30,1)=q2w(30:1:-1)
!         rthcuten(i,12:30,1)=q1(19:1:-1)
!         rqvcuten(i,12:30,1)=q2(19:1:-1)
!        rthcuten(i,1:30,1)=rthcuten(i,1:30,1)-therad3d(i,1:30,1)
 

!        do ilev=1,25
!           raincv(i,1)=raincv(i,1)-rqvcuten(i,ilev,1)/gravity*1000
!        enddo

!        raincv(i,1)=raincv(i,1)/(abs(pfull(1)-pfull(nlev_eff)))*1000

!        raincv(i,1)=amax1(0.0,raincv(i,1))

!        rthcuten(i,6:24,1)=q1(1:19)
!        rqvcuten(i,6:24,1)=q2(1:19)
!        rthcuten=rthcuten*1e-10
!        rqvcuten=rqvcuten*1e-10

!        print*,'rthcuten_resnet', rthcuten(i,:,:)
!        print*, 'rqvcuten_resnet', rqvcuten(i,:,:)
!        print*,'i=',i
!        rthcuten=rthcuten*1e-10*0
!        rqvcuten=rqvcuten*1e-10*0

!        print*,'q1w:',q1w(1:19,i)
      enddo
   enddo

end subroutine 

subroutine ResNet(q1, q2, u3d, v3d, w, &
                  th3d, qv3d, therad3d, &
                  pcps, kms,kme)
    implicit none

    integer , intent(in)           ::   kms, kme


    real,    dimension(kms:kme), intent(in) ::      &
                                        pcps,                            &
                                        qv3d,                           &
                                        therad3d,                          &
                                        th3d,                            &
                                        u3d,                            &
                                        v3d,                            &
                                        w

    real,     dimension(1:lengtho) , intent(out) ::  &
                                                              q1, &
                                                              q2
                                                             ! rthften,rqvften
    integer:: iname, ires
    integer:: iseq
    integer:: i, j
    integer,parameter :: ptop_res=30, ptop_reso=29

    real(r4)::  inputs( length ,  in_channels*seq_len+t0_channels)
    real(r4)::  outputs( lengtho,  out_channels)
    real(r4)::  predictions(lengtho , out_channels)
    real(r4)::  sinput(length, in_channels*seq_len+t0_channels), stemp(length, nkernels )
    real(r4)::  soutput(lengtho, out_channels), spred(lengtho, out_channels)
    real(r4)::  output_gap(nkernels),output_dense(lengtho*out_channels)
    real(r4)::  inputs_max( length, in_channels*seq_len+t0_channels)
    real(r4)::  inputs_min( length, in_channels*seq_len+t0_channels)
    real(r4)::  outputs_max( lengtho, out_channels), outputs_min(lengtho,out_channels)
    real(r4)       :: weights0(nsize,tin_channels,nkernels), bias0(nkernels)
    real(r4)       :: weightsa(nsize,nkernels,nkernels),biasa(nkernels),weightsb(nsize,nkernels,nkernels),biasb(nkernels)
!    real(r4)       :: weights(nsize,nkernels,nkernels,nres*2),bias(nkernels,nres*2)
    real(r4)       :: weights1(nkernels,out_channels*lengtho), bias1(out_channels*lengtho)

    weights0=mldata_weight0%f
    bias0=mldata_bias0%f
    weights1=mldata_weight1%f
    bias1=mldata_bias1%f
!    bias=mldata_biases%f
!    print*,'bias1', mldata_bias1%f
!    do i=1,nres*2  
!       weights(1,:,:,i)=mldata_weights%f(3*i-2,:,:)
!       weights(2,:,:,i)=mldata_weights%f(3*i-1,:,:)
!       print*,'weights=',weights
!       print*,sizeof(mldata_weights)
!       weights(3,:,:,i)=mldata_weights%f(3*i,:,:)
!    enddo

!    print*,'weitghts0',weights0
!    print*, 'weights:', weights(:,:,:,:)

    inputs_max=mldata_imax%f
    inputs_min=mldata_imin%f
    outputs_min=mldata_omin%f
    outputs_max=mldata_omax%f

!    inputs_max(:,6)=inputs_max(:,6)+2000
!    inputs_min(:,6)=inputs_min(:,6)-2000


    sinput(1:ptop_res,1)=qv3d(ptop_res:1:-1)
    sinput(1:ptop_res,2)=th3d(ptop_res:1:-1)
    sinput(1:ptop_res,3)=u3d(ptop_res:1:-1)
    sinput(1:ptop_res,4)=v3d(ptop_res:1:-1)
!    sinput(1:ptop_res,5)=w(:)
    sinput(1:ptop_res,5)=w(ptop_res:1:-1)
    sinput(1:ptop_res,6)=pcps(ptop_res:1:-1)
!    sinput(1:ptop_res,7)=therad3d(ptop_res:1:-1)
!    sinput(1:ptop_res,7)=therad3d(6:ptop_res+5)
!    print*, 'inputs_ori: ', sinput(:,:)

!    sinput=mldata_inre%f(29,:,:)

!    sinput=(sinput+1)*(inputs_max- inputs_min)/2 +inputs_min

!    print*, 'sinput_therad3d: ', sinput(:,7)
!    print*, 'input_max_the:', inputs_max(:,:)
!    print*, 'input_min_the:', inputs_min(:,:)
!     print*, 'inputs: ', sinput(:,:)

!    print*, 'inputs_python: ', mldata_inre%f(29,:,:)

    sinput=(sinput-inputs_min)/(inputs_max-inputs_min)*2.d0-1.d0
!    print*,'inputs_max:',inputs_max(:,:)
!    print*,'inputs_min:',inputs_min(:,:)
!    print*, 'sinput: ', sinput(:,:)

    if (.true.) then
    do i=1, in_channels*seq_len+t0_channels
    do j=1,length
        if (sinput(j,i)>=1.d0+ran)then
            sinput(j,i)=1.d0+ran
        elseif (sinput(j,i)<=-1.d0-ran)then
            sinput(j,i)=-1.d0-ran
        endif
        enddo
    enddo
    endif

    if (.true.)then
    do i=1,seq_len*in_channels+t0_channels
        do j=1,length
        if (inputs_max(j,i)<=1e-29)then
            sinput(j,i)=0.d0
        endif
        enddo
    enddo
    endif

!    print*,'initial_sinput:',sinput(:,:)

    call conv1d(stemp, sinput, weights0,bias0,length,seq_len*in_channels+t0_channels, nkernels)

!    print*,'conv1d:', stemp(:,1)

    do ires =1,nres
!       print*,'nres:',ires
       biasa=mldata_biases%f(:,2*ires-1)
       biasb=mldata_biases%f(:,2*ires)
       weightsa(1,:,:)=mldata_weights%f(3*(2*ires-1)-2,:,:)
       weightsa(2,:,:)=mldata_weights%f(3*(2*ires-1)-1,:,:)
       weightsa(3,:,:)=mldata_weights%f(3*(2*ires-1),:,:)
       weightsb(1,:,:)=mldata_weights%f(3*(2*ires)-2,:,:)
       weightsb(2,:,:)=mldata_weights%f(3*(2*ires)-1,:,:)
       weightsb(3,:,:)=mldata_weights%f(3*(2*ires),:,:)
!    print*,'bias1', mldata_bias1%f
!    do i=1,nres*2  
!       weights(1,:,:,i)=mldata_weights%f(3*i-2,:,:)
!       weights(2,:,:,i)=mldata_weights%f(3*i-1,:,:)
!       print*,'weights=',weights
!       print*,sizeof(mldata_weights)
!       weights(3,:,:,i)=mldata_weights%f(3*i,:,:)
!    enddo

!       print*,'weights1:',weights(1,1,:,2*(ires-1)+1)
!       print*,'weights2:',weights(1,1,:,2*(ires-1)+2)
!       print*,'bias1',bias(:,2*(ires-1)+1)
!       print*,'bias2',biasb(:)
       call ResUnit(stemp,  weightsa(:,:,:), biasa(:),weightsb(:,:,:),biasb(:), length, nkernels )
!       print*, 'resunit:',stemp(:,1)
!    do i=1,nres*2
!       weights(1,:,:,i)=mldata_weights%f(3*i-2,:,:)
!       weights(2,:,:,i)=mldata_weights%f(3*i-1,:,:)
!       print*,'weights=',weights
!       print*,sizeof(mldata_weights)
!       weights(3,:,:,i)=mldata_weights%f(3*i,:,:)
!    enddo

!    do ires =1,nres
!       print*,'nres:',ires
!       print*,'weights1:',weights(1,1,:,2*(ires-1)+1)
!       print*,'weights2:',weights(1,1,:,2*(ires-1)+2)
!       print*,'bias1',bias(:,2*(ires-1)+1)
!       print*,'bias2',bias(:,2*(ires-1)+2)
!       call ResUnit(stemp,  weights(:,:,:,2*(ires-1)+1), bias(:,2*(ires-1)+1),weights(:,:,:,2*(ires-1)+2),bias(:,2*(ires-1)+2), length, nkernels )
!       print*, 'resunit:',stemp(:,1)
!    enddo


    enddo

!    print*, 'resunit:',stemp(:,:)

    call relu(stemp, length, nkernels)

!    print*, 'relu:',stemp(:,:)

    output_gap=stemp(1,:)
    call global_average_pooling_1d(stemp, output_gap, length, nkernels)
!    print*,'pooling:',output_gap
    call dense(output_gap,output_dense,weights1,bias1,"tanh")
!    print*, 'dense:',output_dense
    soutput=reshape(output_dense,(/lengtho,2/),order=(/2,1/))

!    print*,'output_fortran:',soutput
!    soutput=mldata_pred%f(29,:,:)
!    print*,'output_pyhton:', soutput

        ! ResNet complete
    if (.true.)then
    do i=1, out_channels
        do j=1,lengtho
        if (soutput(j,i)>=1.d0+ran)then
            soutput(j,i)=1.d0+ran
        elseif (soutput(j,i)<=-1.d0-ran)then
            soutput(j,i)=-1.d0-ran
        endif
        enddo
    enddo
    endif

!    print*,'remove_fillvalue:',soutput
!    print*,'outputs_max:', outputs_max
!    print*,'outputs_min:', outputs_min
    ! reverse of normalization
    soutput=(soutput+1.d0)*(outputs_max-outputs_min)/2.d0+outputs_min
!    print*,'res_norm:',soutput
    if (.true.)then
    do i=1, out_channels
        do j=1,lengtho
        if (abs(soutput(j,i))<=1e-29)then
            soutput(j,i)=0.d0
        endif
        enddo
    enddo
    endif

!    print*,'re_norm:',soutput

    q1=0.d0
    q2=0.d0

    q1(1:ptop_reso)=soutput(1:ptop_reso,1)
    q2(1:ptop_reso)=soutput(1:ptop_reso,2)
 
!    print*,'q1:',q1
!    print*,'q2:',q2

end subroutine

subroutine conv1d(outputs,inputs_in,kernels,bias, length, nkernels1, nkernels2)
    implicit none
    integer,intent(in):: length,nkernels1,nkernels2
    integer,parameter::size=3,padup=1,paddown=1
    integer:: i,j,k
    real(r4),intent(in):: kernels(size,nkernels1,nkernels2)
    real(r4),intent(in):: bias(nkernels2)
    real(r4),intent(in):: inputs_in(length,nkernels1)
    real(r4)::inputs(padup+length+paddown,nkernels1)
    real(r4),intent(out):: outputs(length,nkernels2)

    inputs=0.d0
    outputs=0.d0
    inputs(1+padup:length+padup,:)=inputs_in
    !print*,'inputs+pad in conv', inputs(:,1)
    ! inputs_size=(:,26,40) kernel_size=(3,40,128)
    do k=1 , nkernels2
        do i=1 ,nkernels1
            do j=1 , length
                outputs(j , k)=outputs(j , k)+ sum(inputs (j: j+size-1, i)* kernels(:, i , k))
            end do
        end do
    end do
    !print*,'outputs bf +bias:', outputs(:,1) 
          do i=1,length
               outputs(i, :) = outputs(i, :)+bias
          enddo

end subroutine

subroutine relu(inputs ,length, nkernels)
    implicit none
    integer,intent(in):: length, nkernels
    integer::k,j
    !real(r4),parameter:: alpha=1.d0
    real(r4),intent(inout):: inputs(length,nkernels)
    do k=1,nkernels
        do j=1, length
            if (inputs(j,k)>0.d0)then
                inputs(j,k)= inputs(j, k)
            else
                inputs(j,k)= alpha*inputs(j,k)
            endif
        enddo
    enddo
end subroutine

subroutine  ResUnit(inputs, kernels1,bias1,kernels2, bias2, length ,nkernels)
    implicit none
    integer,intent(in):: length,nkernels
    integer,parameter::size=3
    real*4,intent(in):: kernels1(size,nkernels,nkernels)
    real*4,intent(in):: kernels2(size,nkernels,nkernels)
    real*4,intent(in):: bias1(nkernels)
    real*4,intent(in):: bias2(nkernels)
    real*4,intent(inout):: inputs(length, nkernels)
    real*4:: outputs(length, nkernels), inputsm(length, nkernels)

    outputs=0.d0
    inputsm=inputs

    call relu(inputs,length,nkernels)
    call conv1d(outputs,inputs,kernels1,bias1,length,nkernels,nkernels)
    inputs=outputs
    outputs=0.d0
    call relu(inputs,length,nkernels)
    call conv1d(outputs,inputs,kernels2,bias2,length,nkernels,nkernels)
    inputs=outputs
    inputs=inputs+inputsm

end subroutine

subroutine global_average_pooling_1d(input_array, output_array, input_size, channels)
    real(4), intent(in) :: input_array(:,:)
    real(4), intent(out) :: output_array(:)
    integer, intent(in) :: input_size, channels

    integer :: i, j
    do i = 1, channels
        output_array(i) = 0.0
        do j = 1, input_size
            output_array(i) = output_array(i) + input_array(j, i)
        end do
        output_array(i) = output_array(i) / real(input_size, 8)
    end do
end subroutine

subroutine dense(input_array, output_array, weight_matrix, bias_array, activation)
    real(4), intent(in) :: input_array(:)
    real(4), intent(out) :: output_array(:)
    real(4), intent(in) :: weight_matrix(:,:)
    real(4), intent(in) :: bias_array(:)
    character(len=*), intent(in) :: activation

    integer :: input_size, output_size
    integer :: i, j

    input_size = size(input_array) !128
    output_size = size(output_array) !38

    ! 计算输出数组的值
    do i = 1, output_size
        output_array(i) = 0.0
        do j = 1, input_size
            output_array(i) = output_array(i) + input_array(j)*weight_matrix(j, i)
        end do
        output_array(i) = output_array(i) + bias_array(i)

        ! 应用激活函数
        if (trim(activation) == "relu") then
            if (output_array(i) < 0.0) output_array(i) = 0.0
        elseif (trim(activation) == "sigmoid") then
            output_array(i) = 1.0 / (1.0 + exp(-output_array(i)))
        elseif (trim(activation) == "tanh") then
            output_array(i)=tanh(output_array(i))
        ! 添加其他激活函数的条件分支
        end if
    end do
end subroutine
#endif

end module module_cu_ResNet_v1

