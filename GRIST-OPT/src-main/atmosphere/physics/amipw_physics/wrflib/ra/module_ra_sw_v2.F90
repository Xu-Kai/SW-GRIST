!wrf:model_layer:physics
!
module module_ra_sw_v2

contains

!------------------------------------------------------------------
   subroutine swrad(dt,rthraten,gsw,xlat,xlong,albedo,            &
                    rho_phy,t3d,qv3d,qc3d,qr3d,                   &
                    qi3d,qs3d,qg3d,p3d,pi3d,dz8w,gmt,             &
                    r,cp,g,julday,                                &
                    xtime,declin,solcon,                          &
                    p_qv,p_qc,p_qr,p_qi,p_qs,p_qg,                &
                    p_first_scalar,                               &
                    radfrq,icloud,degrad,warm_rain,               &
                    ids,ide, jds,jde, kds,kde,                    & 
                    ims,ime, jms,jme, kms,kme,                    &
                    its,ite, jts,jte, kts,kte                     ) 
!------------------------------------------------------------------
   implicit none
!------------------------------------------------------------------
   integer,    intent(in   ) ::        ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte, &
                                       icloud,p_qv,p_qc,p_qr,     &
                                       p_qi,p_qs,p_qg,            &
                                       p_first_scalar
   logical,    intent(in   ) ::        warm_rain

   real, intent(in    )      ::        radfrq,degrad,             &
                                       xtime,declin,solcon
!
   real, dimension( ims:ime, kms:kme, jms:jme ),                  &
         intent(in    ) ::                                   p3d, &
                                                            pi3d, &
                                                             t3d, &
                                                            qv3d, &
                                                            qc3d, &
                                                            qr3d, &
                                                            qi3d, &
                                                            qs3d, &
                                                            qg3d, &
                                                         rho_phy, &
                                                            dz8w

   real, dimension( ims:ime, kms:kme, jms:jme ),                  &
         intent(inout)  ::                              rthraten
!
   real, dimension( ims:ime, jms:jme ),                           &
         intent(in   )  ::                                  xlat, &
                                                           xlong, &
                                                          albedo
!
   real, dimension( ims:ime, jms:jme ),                           &
         intent(inout)  ::                                   gsw
!
   real, intent(in   )   ::                        gmt,r,cp,g,dt
!
   integer, intent(in  ) ::                               julday  
 
! local vars
 
   real, dimension( kts:kte ) ::                                  &
                                                          tten1d, &
                                                          rho01d, &
                                                             p1d, &
                                                              dz, &
                                                             t1d, &
                                                            qv1d, &
                                                            qc1d, &
                                                            qr1d, &
                                                            qi1d, &
                                                            qs1d, &
                                                            qg1d
!
   real::      xlat0,xlong0,alb0,gsw0
!
   integer :: i,j,k,nk

!------------------------------------------------------------------
   j_loop: do j=jts,jte
   i_loop: do i=its,ite

! reverse vars 
         do k=kts,kte
            qv1d(k)=0.
            qc1d(k)=0.
            qr1d(k)=0.
            qi1d(k)=0.
            qs1d(k)=0.
            qg1d(k)=0.
         enddo

         do k=kts,kte
            nk=kme-1-k+kms
            tten1d(k)=0.

            t1d(k)=t3d(i,nk,j)
            p1d(k)=p3d(i,nk,j)
            rho01d(k)=rho_phy(i,nk,j)
            dz(k)=dz8w(i,nk,j)
         enddo

         if (p_qv .ge. p_first_scalar) then
            do k=kts,kte
               nk=kme-1-k+kms
               qv1d(k)=qv3d(i,nk,j)
               qv1d(k)=max(0.,qv1d(k))
            enddo
         endif

         if (p_qc .ge. p_first_scalar) then
            do k=kts,kte
               nk=kme-1-k+kms
               qc1d(k)=qc3d(i,nk,j)
               qc1d(k)=max(0.,qc1d(k))
            enddo
         endif

         if (p_qr .ge. p_first_scalar) then
            do k=kts,kte
               nk=kme-1-k+kms
               qr1d(k)=qr3d(i,nk,j)
               qr1d(k)=max(0.,qr1d(k))
            enddo
         endif

!
         if (p_qi .ge. p_first_scalar) then
            do k=kts,kte          
               nk=kme-1-k+kms
               qi1d(k)=qi3d(i,nk,j)
               qi1d(k)=max(0.,qi1d(k))
            enddo
         else
            if (.not. warm_rain) then
               do k=kts,kte          
               if(t1d(k) .lt. 273.15) then
                  qi1d(k)=qc1d(k)
                  qc1d(k)=0.
                  qs1d(k)=qr1d(k)
                  qr1d(k)=0.
               endif
               enddo
            endif
         endif

         if (p_qs .ge. p_first_scalar) then
            do k=kts,kte          
               nk=kme-1-k+kms
               qs1d(k)=qs3d(i,nk,j)
               qs1d(k)=max(0.,qs1d(k))
            enddo
         endif

         if (p_qg .ge. p_first_scalar) then
            do k=kts,kte          
               nk=kme-1-k+kms
               qg1d(k)=qg3d(i,nk,j)
               qg1d(k)=max(0.,qg1d(k))
            enddo
         endif

         xlat0=xlat(i,j)
         xlong0=xlong(i,j)
         alb0=albedo(i,j)

         call swpara(tten1d,gsw0,xlat0,xlong0,alb0,              &
                     t1d,qv1d,qc1d,qr1d,qi1d,qs1d,qg1d,p1d,      &
                     xtime,gmt,rho01d,dz,                        &
                     r,cp,g,declin,solcon,                       &
                     radfrq,icloud,degrad,                       &
                     kts,kte                                     )

         gsw(i,j)=gsw0
         do k=kts,kte          
            nk=kme-1-k+kms
            rthraten(i,k,j)=rthraten(i,k,j)+tten1d(nk)/pi3d(i,k,j)
         enddo
!
   enddo i_loop
   enddo j_loop                                          

   end subroutine swrad

!------------------------------------------------------------------
   subroutine swpara(tten,gsw,xlat,xlong,albedo,                  &
                     t,qv,qc,qr,qi,qs,qg,p,                       &
                     xtime, gmt, rho0, dz,                        &
                     r,cp,g,declin,solcon,                        &
                     radfrq,icloud,degrad,                        &
                     kts,kte                                      )
!------------------------------------------------------------------
!     to calculate short-wave absorption and scattering in clear
!     air and reflection and absorption in cloud layers (stephens,
!     1984)
!     changes:
!       reduce effects of ice clouds and precip on liquid water path
!       add effect of graupel
!------------------------------------------------------------------

  integer, intent(in ) ::                 kts,kte
!
  real, dimension( kts:kte ), intent(in   )  ::                   &
                                                            rho0, &
                                                               t, &
                                                               p, &
                                                              dz, &
                                                              qv, &
                                                              qc, &
                                                              qr, &
                                                              qi, &
                                                              qs, &
                                                              qg

   real, dimension( kts:kte ), intent(inout)::              tten
!
   real, intent(in  )   ::               xtime,gmt,r,cp,g,declin, &
                                        solcon,xlat,xlong,albedo, &
                                                  radfrq, degrad
!
   real, intent(inout)  ::                                   gsw
!
! local vars
!
   real, dimension( kts:kte+1 ) ::                         sdown

   real, dimension( kts:kte )   ::                          xlwp, &
                                                            xatp, &
                                                            xwvp, &
                                                              ro
!
   real, dimension( 4, 5 ) ::                             albtab, &
                                                          abstab

   real, dimension( 4    ) ::                             xmuval

!------------------------------------------------------------------

      data albtab/0.,0.,0.,0., &
           69.,58.,40.,15.,    &
           90.,80.,70.,60.,    &
           94.,90.,82.,78.,    &
           96.,92.,85.,80./

      data abstab/0.,0.,0.,0., &
           0.,2.5,4.,5.,       &
           0.,2.6,7.,10.,      &
           0.,3.3,10.,14.,     &
           0.,3.7,10.,15./

      data xmuval/0.,0.2,0.5,1.0/

      gsw=0.0
 
      soltop=solcon
      xt24=amod(xtime+radfrq*0.5,1440.)
      tloctm=gmt+xt24/60.+xlong/15.
      hrang=15.*(tloctm-12.)*degrad
      xxlat=xlat*degrad
      csza=sin(xxlat)*sin(declin)+cos(xxlat)*cos(declin)*cos(hrang)

!     return if night
      if(csza.le.1.e-9)goto 7
!
      do k=kts, kte

! p in the unit of 10mb
         ro(k)=p(k)/(r*t(k))
         xwvp(k)=ro(k)*qv(k)*dz(k)*1000.
! kg/m**2
          xatp(k)=ro(k)*dz(k)
      enddo
!
!     g/m**2
!     reduce weight of liquid and ice in short-wave scheme
!     add graupel effect (assumed same as rain)
!
      if (icloud.eq.0)then
         do k=kts, kte
            xlwp(k)=0.
         enddo
      else
         do k=kts, kte
            xlwp(k)=ro(k)*1000.*dz(k)*(qc(k)+0.1*qi(k)+0.05* &
                    qr(k)+0.02*qs(k)+0.05*qg(k))
         enddo
      endif
!
      xmu=csza
      sdown(1)=soltop*xmu
!     set ww (g/m**2) liquid water path integrated down
!     set uv (g/m**2) water vapor path integrated down
      ww=0.
      uv=0.
      oldalb=0.
      oldabc=0.
      totabs=0.
!     contributions due to clear air and cloud
      dsca=0.
      dabs=0.
      dscld=0.
!
      do 200 k=kts,kte
         ww=ww+xlwp(k)
         uv=uv+xwvp(k)
!     wgm is ww/cos(theta) (g/m**2)
!     ugcm is uv/cos(theta) (g/cm**2)
         wgm=ww/xmu
         ugcm=uv*0.0001/xmu
!
         oldabs=totabs
!     water vapor absorption as in lacis and hansen (1974)
         totabs=2.9*ugcm/((1.+141.5*ugcm)**0.635+5.925*ugcm)
!     approximate rayleigh + aerosol scattering
         xsca=1.e-5*xatp(k)/xmu
!     layer vapor absorption done first
         xabs=(totabs-oldabs)*(sdown(1)-dscld-dsca)/sdown(k)
         if(xabs.lt.0.)xabs=0.
!
         alw=alog10(wgm+1.)
         if(alw.gt.3.999)alw=3.999
!
         do ii=1,3
            if(xmu.gt.xmuval(ii))then
              iil=ii
              iu=ii+1
              xi=(xmu-xmuval(ii))/(xmuval(ii+1)-xmuval(ii))+float(iil)
            endif
         enddo
!
         jjl=ifix(alw)+1
         ju=jjl+1
         yj=alw+1.
!     cloud albedo
         alba=(albtab(iu,ju)*(xi-iil)*(yj-jjl)   &
              +albtab(iil,ju)*(iu-xi)*(yj-jjl)   &
              +albtab(iu,jjl)*(xi-iil)*(ju-yj)   &
              +albtab(iil,jjl)*(iu-xi)*(ju-yj))  &
             /((iu-iil)*(ju-jjl))
!     cloud absorption
         absc=(abstab(iu,ju)*(xi-iil)*(yj-jjl)   &
              +abstab(iil,ju)*(iu-xi)*(yj-jjl)   &
              +abstab(iu,jjl)*(xi-iil)*(ju-yj)   &
              +abstab(iil,jjl)*(iu-xi)*(ju-yj))  &
             /((iu-iil)*(ju-jjl))
!     layer albedo and absorption
         xalb=(alba-oldalb)*(sdown(1)-dsca-dabs)/sdown(k)
         xabsc=(absc-oldabc)*(sdown(1)-dsca-dabs)/sdown(k)
         if(xalb.lt.0.)xalb=0.
         if(xabsc.lt.0.)xabsc=0.
         dscld=dscld+(xalb+xabsc)*sdown(k)*0.01
         dsca=dsca+xsca*sdown(k)
         dabs=dabs+xabs*sdown(k)
         oldalb=alba
         oldabc=absc
!     layer transmissivity
         trans0=100.-xalb-xabsc-xabs*100.-xsca*100.
         if(trans0.lt.1.)then
           ff=99./(xalb+xabsc+xabs*100.+xsca*100.)
           xalb=xalb*ff
           xabsc=xabsc*ff
           xabs=xabs*ff
           xsca=xsca*ff
           trans0=1.
         endif
         sdown(k+1)=amax1(1.e-9,sdown(k)*trans0*0.01)
         tten(k)=sdown(k)*(xabsc+xabs*100.)*0.01/( &
                 ro(k)*cp*dz(k))
  200   continue
!
        gsw=(1.-albedo)*sdown(kte+1)

    7 continue
!
   end subroutine swpara

!====================================================================
   subroutine swinit(rthraten,rthratensw,restart,                   &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte                   )
!--------------------------------------------------------------------
   implicit none
!--------------------------------------------------------------------
   logical , intent(in)           :: restart
   integer , intent(in)           :: ids, ide, jds, jde, kds, kde,  &
                                     ims, ime, jms, jme, kms, kme,  &
                                     its, ite, jts, jte, kts, kte

   real , dimension( ims:ime , kms:kme , jms:jme ) , intent(inout) ::        &
                                                          rthraten, &
                                                        rthratensw
   integer :: i, j, k, itf, jtf, ktf

   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)

   !if(.not.restart)then
   !  do j=jts,jtf
   !  do k=kts,ktf
   !  do i=its,itf
        rthraten   = 0.
        rthratensw = 0.
   !  enddo
   !  enddo
   !  enddo
   !endif

   end subroutine swinit


end module module_ra_sw_v2
