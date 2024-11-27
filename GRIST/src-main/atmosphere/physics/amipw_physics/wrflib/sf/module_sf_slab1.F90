!wrf:model_layer:physics
!
module module_sf_slab1

   !---specify constants and layers for soil model
   !---soil diffusion constant set (m^2/s)

   real, parameter :: difsl=5.e-7

   !---factor to make soil step more conservative

   real , parameter :: soilfac=1.25

contains

!----------------------------------------------------------------
   subroutine slab(t3d,qv3d,p3d,flhc,flqc,                      &
                   psfc,xland,tmn,hfx,qfx,lh,tsk,qsfc,chklowq,  &
                   gsw,glw,capg,thc,snowc,emiss,mavail,         &
                   deltsm,rovcp,xlv,dtmin,ifsnow,               &
                   svp1,svp2,svp3,svpt0,ep2,                    &
                   karman,eomeg,stbolt,                         &
                   tslb,zs,dzs,num_soil_layers,radiation,       &
                   ids,ide, jds,jde, kds,kde,                   &
                   ims,ime, jms,jme, kms,kme,                   &
                   its,ite, jts,jte, kts,kte                    )
!----------------------------------------------------------------
    implicit none
!----------------------------------------------------------------
!                                                                        
!     subroutine slab calculates the ground temperature tendency 
!     according to the residual of the surface energy budget           
!     (blackadar, 1978b).                                              
!                                                                      
!     changes:                                                         
!          for soil sub-timesteps update surface hfx and qfx as tg     
!          changes to prevent possible instability for long model      
!          steps (dt > ~200 sec).                                      
!                                                                      
!          put snow cover check on soil sub-timesteps                  
!                                                                      
!          make upper limit on soil sub-step length more conservative  
!                                                                      
!----------------------------------------------------------------          
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- p3d         3d pressure (pa)
!-- flhc        exchange coefficient for heat (m/s)
!-- flqc        exchange coefficient for moisture (m/s)
!-- psfc        surface pressure (pa)
!-- xland       land mask (1 for land, 2 for water)
!-- tmn         soil temperature at lower boundary (k)
!-- hfx         upward heat flux at the surface (w/m^2)
!-- qfx         upward moisture flux at the surface (kg/m^2/s)
!-- lh          latent heat flux at the surface (w/m^2)
!-- tsk         surface temperature (k)
!-- gsw         downward short wave flux at ground surface (w/m^2)      
!-- glw         downward long wave flux at ground surface (w/m^2)
!-- capg        heat capacity for soil (j/k/m^3)
!-- thc         thermal inertia (cal/cm/k/s^0.5)
!-- snowc       flag indicating snow coverage (1 for snow cover)
!-- emiss       surface emissivity (between 0 and 1)
!-- deltsm      time step (second)
!-- rovcp       r/cp
!-- xlv         latent heat of melting (j/kg)
!-- dtmin       time step (minute)
!-- ifsnow      ifsnow=1 for snow-cover effects
!-- svp1        constant for saturation vapor pressure (kpa)
!-- svp2        constant for saturation vapor pressure (dimensionless)
!-- svp3        constant for saturation vapor pressure (k)
!-- svpt0       constant for saturation vapor pressure (k)
!-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
!-- ep2         constant for specific humidity calculation 
!               (r_d/r_v) (dimensionless)
!-- karman      von karman constant
!-- eomeg       angular velocity of earth's rotation (rad/s)
!-- stbolt      stefan-boltzmann constant (w/m^2/k^4)
!-- tslb        soil temperature in 5-layer model
!-- zs          depths of centers of soil layers
!-- dzs         thicknesses of soil layers
!-- num_soil_layers   the number of soil layers
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
!----------------------------------------------------------------
   integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,  &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte

   integer, intent(in)       ::     num_soil_layers
   logical, intent(in)       ::     radiation

   integer,  intent(in   )   ::     ifsnow

!
   real,     intent(in   )   ::     dtmin,xlv,rovcp,deltsm

   real,     intent(in )     ::     svp1,svp2,svp3,svpt0
   real,     intent(in )     ::     ep2,karman,eomeg,stbolt

   real,     dimension( ims:ime , 1:num_soil_layers, jms:jme ), &
             intent(inout)   :: tslb

   real,     dimension(1:num_soil_layers), intent(in)::zs,dzs

   real,    dimension( ims:ime, kms:kme, jms:jme )            , &
            intent(in   )    ::                           qv3d, &
                                                           p3d, &
                                                           t3d
!
   real,    dimension( ims:ime, jms:jme )                     , &
            intent(in   )    ::                          snowc, &
                                                         xland, &
                                                         emiss, &
                                                        mavail, &
                                                           tmn, &
                                                           gsw, &
                                                           glw, &
                                                           thc

!chklowq is declared as memory size
!
   real,    dimension( ims:ime, jms:jme )                     , &
            intent(inout)    ::                            hfx, &
                                                           qfx, &
                                                            lh, &
                                                          capg, &
                                                           tsk, &
                                                          qsfc, &
                                                       chklowq

   real,     dimension( ims:ime, jms:jme )                    , &
             intent(in   )               ::               psfc
!
   real,    dimension( ims:ime, jms:jme ), intent(inout) ::     &
                                                          flhc, &
                                                          flqc

! local vars

   real,     dimension( its:ite ) ::                      qv1d, &
                                                           p1d, &
                                                           t1d
   integer ::  i,j

   do j=jts,jte

      do i=its,ite
         t1d(i) =t3d(i,1,j)
         qv1d(i)=qv3d(i,1,j)
         p1d(i) =p3d(i,1,j)
      enddo

      call slab1d(j,t1d,qv1d,p1d,flhc(ims,j),flqc(ims,j),       &
           psfc(its,j),xland(ims,j),tmn(ims,j),hfx(ims,j),      &
           qfx(ims,j),tsk(ims,j),qsfc(ims,j),chklowq(ims,j),    &
           lh(ims,j),gsw(ims,j),glw(ims,j),                     &
           capg(ims,j),thc(ims,j),snowc(ims,j),emiss(ims,j),    &
           mavail(ims,j),deltsm,rovcp,xlv,dtmin,ifsnow,         &
           svp1,svp2,svp3,svpt0,ep2,karman,eomeg,stbolt,        &
           tslb(ims,1,j),zs,dzs,num_soil_layers,radiation,      &
           ids,ide, jds,jde, kds,kde,                           &
           ims,ime, jms,jme, kms,kme,                           &
           its,ite, jts,jte, kts,kte                            )

   enddo

   end subroutine slab

!----------------------------------------------------------------
   subroutine slab1d(j,t1d,qv1d,p1d,flhc,flqc,                  &
                   psfcpa,xland,tmn,hfx,qfx,tsk,qsfc,chklowq,   &
                   lh,gsw,glw,capg,thc,snowc,emiss,mavail,      &
                   deltsm,rovcp,xlv,dtmin,ifsnow,               &
                   svp1,svp2,svp3,svpt0,ep2,                    &
                   karman,eomeg,stbolt,                         &
                   tslb2d,zs,dzs,num_soil_layers,radiation,     &
                   ids,ide, jds,jde, kds,kde,                   &
                   ims,ime, jms,jme, kms,kme,                   &
                   its,ite, jts,jte, kts,kte                    )
!----------------------------------------------------------------
    implicit none
!----------------------------------------------------------------
!                                                                        
!     subroutine slab calculates the ground temperature tendency 
!     according to the residual of the surface energy budget           
!     (blackadar, 1978b).                                              
!                                                                      
!     changes:                                                         
!          for soil sub-timesteps update surface hfx and qfx as tg     
!          changes to prevent possible instability for long model      
!          steps (dt > ~200 sec).                                      
!                                                                      
!          put snow cover check on soil sub-timesteps                  
!                                                                      
!          make upper limit on soil sub-step length more conservative  
!                                                                      
!----------------------------------------------------------------          

   integer,  intent(in   )   ::     ids,ide, jds,jde, kds,kde,  &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte,j 

   integer , intent(in)      ::     num_soil_layers
   logical,  intent(in   )   ::     radiation

   integer,  intent(in   )   ::     ifsnow
!
   real,     intent(in   )   ::     dtmin,xlv,rovcp,deltsm

   real,     intent(in )     ::     svp1,svp2,svp3,svpt0
   real,     intent(in )     ::     ep2,karman,eomeg,stbolt

   real,     dimension( ims:ime , 1:num_soil_layers ),          &
             intent(inout)   :: tslb2d

   real,     dimension(1:num_soil_layers), intent(in)::zs,dzs

!
   real,    dimension( ims:ime )                              , &
            intent(inout)    ::                            hfx, &
                                                           qfx, &
                                                            lh, &
                                                          capg, &
                                                           tsk, &
                                                          qsfc, &
                                                       chklowq
!
   real,    dimension( ims:ime )                              , &
            intent(in   )    ::                          snowc, &
                                                         xland, &
                                                         emiss, &
                                                        mavail, &
                                                           tmn, &
                                                           gsw, &
                                                           glw, &
                                                           thc
!
   real,    dimension( its:ite )                              , &
            intent(in   )    ::                           qv1d, &
                                                           p1d, &
                                                           t1d
!
   real,     dimension( its:ite )                             , &
             intent(in   )               ::             psfcpa

!
   real,    dimension( ims:ime ), intent(inout) ::              &
                                                          flhc, &
                                                          flqc
! local vars

   real,    dimension( its:ite )          ::              psfc

   real,    dimension( its:ite )          ::                    &
                                                           thx, &
                                                            qx, &
                                                          scr3 

   real,    dimension( its:ite )          ::            dthgdt, &
                                                           tg0, &
                                                         thtmn, &
                                                          xld1, &
                                                         tscvn, &
                                                          oltg, &
                                                        upflux, &
                                                            hm, &
                                                          rnet, &
                                                         xinet, &
                                                            qs, &
                                                         dtsdt
!
   real, dimension( its:ite, num_soil_layers )        :: flux
!
   integer :: i,k,nsoil,itsoil,l,nk,radswtch
   real    :: ps,ps1,xldcol,tskx,rnsoil,rhog1,rhog2,rhog3,lamdag
   real    :: thg,esg,qsg,hfxt,qfxt,cs,csw,lamg(4),thcon,pl
 
!----------------------------------------------------------------------          
!-----determine if any points in column are land (rather than ocean)             
!       points.  if not, skip down to the print statements since ocean           
!       surface temperatures are not allowed to change.                          
!                                                                                
! from sfcrad   
!----------------------------------------------------------------------
   data csw/4.183e6/
   data lamg/1.407e-8, -1.455e-5, 6.290e-3, 0.16857/

   do i=its,ite
! in cmb
      psfc(i)=psfcpa(i)/1000.
   enddo


      do i=its,ite
! pl cmb
         pl=p1d(i)/1000.
         scr3(i)=t1d(i)
         thcon=(100./pl)**rovcp
         thx(i)=scr3(i)*thcon
         qx(i)=0.
      enddo

!     if(idry.eq.1) goto 81
      do i=its,ite
         qx(i)=qv1d(i)
      enddo
   81 continue

!
!-----the slab thermal capacity capg(i) are dependent on:
!     thc(i) - soil thermal inertial, only.
!
      do i=its,ite
         capg(i)=3.298e6*thc(i)
         if(num_soil_layers .gt. 1)then

! capg represents soil heat capacity (j/k/m^3) when difsl=5.e-7 (m^2/s)
! to give a correct thermal inertia (=capg*difsl^0.5)

            capg(i)=5.9114e7*thc(i)
         endif
      enddo
!        
      xldcol=2.0                                                                 
      do 10 i=its,ite
        xldcol=amin1(xldcol,xland(i))                                          
   10 continue                                                                   
!                                                                                
      if(xldcol.gt.1.5)goto 90                                                   
!                                                                                
!                                                                                
!-----convert slab temperature to potential temperature and                      
!     set xld1(i) = 0. for ocean points:                                         
!                                                                                
!                                                                                
      do 20 i=its,ite
        if((xland(i)-1.5).ge.0)then                                            
          xld1(i)=0.                                                             
        else                                                                     
          xld1(i)=1.                                                             
        endif                                                                    
   20 continue                                                                   
!                                                                                
!-----convert 'tsk(thetag)' to 'tg' for 'iup' calculation ....                   
!       if we are using the blackadar multi-level (high-resolution)              
!       pbl model                                                                
!                                                                                
      do 50 i=its,ite
        if(xld1(i).lt.0.5)goto 50                                                

! ps cmb
        ps=psfc(i)

! tsk is temperature at gound sfc
!       tg0(i)=tsk(i)*(ps*0.01)**rovcp                                         
        tg0(i)=tsk(i)
   50 continue                                                                   
!                                                                                
!-----compute the surface energy budget:                                         
!                                                                                
!     if(isoil.eq.1)nsoil=1                                                      
      if(num_soil_layers .gt. 1)nsoil=1                                                      


      if (radiation) then
        radswtch=1
      else
        radswtch=0
      endif

      do 70 i=its,ite
        if(xld1(i).lt.0.5)goto 70
        oltg(i)=tsk(i)*(100./psfc(i))**rovcp
        upflux(i)=radswtch*stbolt*tg0(i)**4                            
        xinet(i)=emiss(i)*(glw(i)-upflux(i))    
        rnet(i)=gsw(i)+xinet(i)                                                
        hm(i)=1.18*eomeg*(tg0(i)-tmn(i))                                       
!       moisture flux calculated here (overwrites sfc layer value for land)
                esg=svp1*exp(svp2*(tg0(i)-svpt0)/(tg0(i)-svp3))
                qsg=ep2*esg/(psfc(i)-esg)
                qfx(i)=flqc(i)*(qsg-qx(i))
                lh(i)=qfx(i)*xlv
        qs(i)=hfx(i)+qfx(i)*xlv                                
!       if(isoil.eq.0)then                                                       
        if(num_soil_layers .eq. 1)then                                                       
          dthgdt(i)=(rnet(i)-qs(i))/capg(i)-hm(i)                              
        else
          dthgdt(i)=0.                                                           
        endif                                                                    
   70 continue                                                                   
!     if(isoil.eq.1)then                                                         
      if(num_soil_layers .gt. 1)then                                                         
        nsoil=1+ifix(soilfac*4*difsl/dzs(1)*deltsm/dzs(1))   
        rnsoil=1./float(nsoil)                                                   
!                                                                                
!     soil sub-timestep                                                          
!                                                                                
        do itsoil=1,nsoil                                                        
          do i=its,ite
             do l=1,num_soil_layers-1
              if(xld1(i).lt.0.5)goto 75                                          
              if(l.eq.1.and.itsoil.gt.1)then                                     
                ps1=(psfc(i)*0.01)**rovcp    

! for rk scheme a and b are the same
                ps=psfc(i)
                thg=tslb2d(i,1)/ps1                                              
                esg=svp1*exp(svp2*(tslb2d(i,1)-svpt0)/(tslb2d(i,1) & 
                    -svp3))                                                      
                qsg=ep2*esg/(ps-esg)                                             
!     update fluxes for new ground temperature                                   
                hfxt=flhc(i)*(thg-thx(i))                                     
                qfxt=flqc(i)*(qsg-qx(i))
                qs(i)=hfxt+qfxt*xlv                                
!     sum hfx and qfx over soil timesteps                                        
                hfx(i)=hfx(i)+hfxt                                           
                qfx(i)=qfx(i)+qfxt                                           
              endif                                                              
              flux(i,1)=rnet(i)-qs(i)                                            
              flux(i,l+1)=-difsl*capg(i)*(tslb2d(i,l+1)-tslb2d(i,l))/( & 
                          zs(l+1)-zs(l))                                         
              dtsdt(i)=-(flux(i,l+1)-flux(i,l))/(dzs(l)*capg(i))               
              tslb2d(i,l)=tslb2d(i,l)+dtsdt(i)*deltsm*rnsoil                     
              if(ifsnow.eq.1.and.l.eq.1)then                              
                if((snowc(i).gt.0..and.tslb2d(i,1).gt.273.16))then             
                  tslb2d(i,1)=273.16                                             
                endif                                                            
              endif                                                              
              if(l.eq.1)dthgdt(i)=dthgdt(i)+rnsoil*dtsdt(i)                      
              if(itsoil.eq.nsoil.and.l.eq.1)then                                 
!     average hfx and qfx over soil timesteps for output to pbl                  
                hfx(i)=hfx(i)*rnsoil                                         
                qfx(i)=qfx(i)*rnsoil                                         
                lh(i)=qfx(i)*xlv
              endif                                                              
   75         continue                                                           
            enddo                                                                
          enddo                                                                  
        enddo                                                                    
      endif                                                                      
!                                                                                
      do 80 i=its,ite
        if(xld1(i).lt.0.5) goto 80                                                
        tskx=tg0(i)+deltsm*dthgdt(i)                                             

! tsk is temperature
!       tsk(i)=tskx*(100./ps1)**rovcp                                          
        tsk(i)=tskx
   80 continue                                                                   

!                                                                                
!-----modify the the ground temperature if the snow cover effects are            
!     considered: limit the ground temperature under 0 c.                        
!                                                                                
      if(ifsnow.eq.0)goto 90                                              
      do 85 i=its,ite
        if(xld1(i).lt.0.5)goto 85                                                
!       ps1=(psfc(i)*0.01)**rovcp             
!       tscvn(i)=tsk(i)*ps1                                            
        tscvn(i)=tsk(i)
        if((snowc(i).gt.0..and.tscvn(i).gt.273.16))then                        
          tscvn(i)=273.16                                                        
        else                                                                     
          tscvn(i)=tscvn(i)                                                      
        endif                                                                    
!       tsk(i)=tscvn(i)/ps1                                                    
        tsk(i)=tscvn(i)
   85 continue                                                                   
!                                                                                
   90 continue                                                                   
      do i=its,ite
! qsfc and chklowq needed by eta pbl
        qsfc(i)=qx(i)+qfx(i)/flqc(i)
        chklowq(i)=mavail(i)
      enddo
!                                                                                
  140 continue                                                                   

   end subroutine slab1d

!================================================================
   subroutine slabinit(tsk,tmn,                                 &
                       tslb,zs,dzs,num_soil_layers,             &
                       restart,                                 &
                       ids,ide, jds,jde, kds,kde,               &
                       ims,ime, jms,jme, kms,kme,               &
                       its,ite, jts,jte, kts,kte                )
!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
   logical , intent(in)      ::      restart
   integer, intent(in   )    ::      ids,ide, jds,jde, kds,kde, &
                                     ims,ime, jms,jme, kms,kme, &
                                     its,ite, jts,jte, kts,kte

   integer, intent(in   )    ::      num_soil_layers
!   
   real,     dimension( ims:ime , 1:num_soil_layers , jms:jme ), intent(inout) :: tslb

   real,     dimension(1:num_soil_layers), intent(in)  ::  zs,dzs

   real,    dimension( ims:ime, jms:jme )                     , &
            intent(in)    ::                            tsk, &
                                                           tmn
!  locar var

   integer                   ::      l,j,i,itf,jtf
!----------------------------------------------------------------
 
   itf=min0(ite,ide-1)
   jtf=min0(jte,jde-1)

   end subroutine slabinit

!-------------------------------------------------------------------          

end module module_sf_slab1
