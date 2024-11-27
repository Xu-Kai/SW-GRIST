!wrf:model_layer:physics
!

module module_mp_lin_v2

!   use     module_wrf_error
   use grist_handle_error, only: endrun

   integer :: wrf_err_message = 6
!
   real    , parameter, private ::       rh = 1.0
!  real    , parameter, private ::   episp0 = 0.622*611.21
   real    , parameter, private ::     xnor = 8.0e6
   real    , parameter, private ::     xnos = 3.0e6

!  lin
!  real    , parameter, private ::     xnog = 4.0e4
!  real    , parameter, private ::     rhograul = 917.

! hobbs
  real     , parameter, private ::     xnog = 4.0e6
  real     , parameter, private ::     rhograul = 400.

!
  real     , parameter, private ::                              &
             qi0 = 1.0e-3, ql0 = 7.0e-4, qs0 = 6.0e-4,          &
             xmi50 = 4.8e-10, xmi40 = 2.46e-10,                 &
             constb = 0.8, constd = 0.25,                       &
             o6 = 1./6.,  cdrag = 0.6,                          &
             avisc = 1.49628e-6, adiffwv = 8.7602e-5,           &
             axka = 1.4132e3, di50 = 1.0e-4, xmi = 4.19e-13,    &
             cw = 4.187e3, vf1s = 0.78, vf2s = 0.31,            &
             xni0 = 1.0e-2, xmnin = 1.05e-18, bni = 0.5,        &
             ci = 2.093e3
contains

!-------------------------------------------------------------------
!  lin et al., 1983, jam, 1065-1092, and
!  rutledge and hobbs, 1984, jas, 2949-2972
!-------------------------------------------------------------------
  subroutine lin_et_al(th, qv, ql, qr, qi, qs, qg,                 &
                       qrold, qsold, qgold,                        &
                       rho, pii, p, rainnc, rainncv,               &
                       dt_in, z, ht, dz8w,                         &
                       grav, cp, rair, rvapor,                     &
                       xls, xlv, xlf, rhowater, rhosnow,           &
                       ep2,svp1,svp2,svp3,svpt0,                   &
                       p_qi,p_qs,p_qg,                             &
                       p_first_scalar,                             &
                       ids,ide, jds,jde, kds,kde,                  & ! domain dims
                       ims,ime, jms,jme, kms,kme,                  & ! memory dims
                       its,ite, jts,jte, kts,kte                   ) ! tile   dims
!-------------------------------------------------------------------
  implicit none
!-------------------------------------------------------------------
!
! shuhua 12/17/99
!
  integer,      intent(in   )    ::   ids,ide, jds,jde, kds,kde , &
                                      ims,ime, jms,jme, kms,kme , &
                                      its,ite, jts,jte, kts,kte , &
                                      p_qi,p_qs,p_qg,             &
                                      p_first_scalar

  real, dimension( ims:ime , kms:kme , jms:jme ),                 &
        intent(inout) ::                                          &
                                                              th, &
                                                              qv, &
                                                              ql, &
                                                              qr, &
                                                              qi, &
                                                              qs, &
                                                              qg
!
  real, dimension( ims:ime , kms:kme , jms:jme ),                 &
        intent(in   ) ::                                          &
                                                           qrold, &
                                                           qsold, &
                                                           qgold, &
                                                             rho, &
                                                             pii, &
                                                               p, &
                                                            dz8w

  real, dimension( ims:ime , kms:kme , jms:jme ),                 &
        intent(in   ) ::                                       z

  real, dimension( ims:ime , jms:jme ),                           &
        intent(inout) ::                                  rainnc, &
                                                         rainncv


  real , dimension( ims:ime , jms:jme ) , intent(in) ::       ht

  real, intent(in   ) ::                                   dt_in, &
                                                            grav, &
                                                            rair, &
                                                          rvapor, &
                                                              cp, &
                                                             xls, &
                                                             xlv, &
                                                             xlf, &
                                                        rhowater, &
                                                         rhosnow

  real, intent(in   ) :: ep2,svp1,svp2,svp3,svpt0

! local var

  integer             ::                            min_q, max_q

  real, dimension( its:ite , jts:jte )                            &
                               ::        rain, snow, graupel,ice

  real, dimension( kts:kte )   ::                  qvz, qlz, qrz, &
                                                   qiz, qsz, qgz, &
                                                  qrzold, qszold, &
                                                     qgzold, thz, &
                                                     tothz, rhoz, &
                                                   orhoz, sqrhoz, &
                                                        prez, zz, &
                                                        dzw
!
  real    ::         dt, pptrain, pptsnow, pptgraul, rhoe_s,      &
                     gindex, pptice

  integer ::               i,j,k
!
   dt=dt_in
   rhoe_s=1.29
!
   if (p_qi .lt. p_first_scalar .or. p_qs .lt. p_first_scalar) then
      call endrun ( 'module_mp_lin: improper use of lin et al scheme; no ice phase. please chose another one.')
   endif

   gindex=1.0

   if (p_qg .lt. p_first_scalar) gindex=0.

   j_loop:  do j = jts, jte
   i_loop:  do i = its, ite
!
!- write data from 3-d to 1-d
!
   do k = kts, kte
      qvz(k)=qv(i,k,j)
      qlz(k)=ql(i,k,j)
      qrz(k)=qr(i,k,j)
      qrzold(k)=qrold(i,k,j)
      thz(k)=th(i,k,j)
!
      rhoz(k)=rho(i,k,j)
      orhoz(k)=1./rhoz(k)
      prez(k)=p(i,k,j)
      sqrhoz(k)=sqrt(rhoe_s*orhoz(k))
      tothz(k)=pii(i,k,j)
      zz(k)=z(i,k,j)
      dzw(k)=dz8w(i,k,j)
   end do

   if (p_qi .ge. p_first_scalar) then
      do k = kts, kte
         qiz(k)=qi(i,k,j)
      enddo
   else
      do k = kts, kte
         qiz(k)=0.
      enddo
   endif

   if (p_qs .ge. p_first_scalar) then
      do k = kts, kte
         qsz(k)=qs(i,k,j)
         qszold(k)=qsold(i,k,j)
      enddo
   else
      do k = kts, kte
         qsz(k)=0.
         qszold(k)=0.
      enddo
   endif

   if (p_qg .ge. p_first_scalar) then
      do k = kts, kte
         qgz(k)=qg(i,k,j)
         qgzold(k)=qgold(i,k,j)
      enddo
   else
      do k = kts, kte
         qgz(k)=0.
         qgzold(k)=0.
      enddo
   endif
!
   pptrain=0.
   pptsnow=0.
   pptgraul=0.
   pptice=0.
   call clphy1d(    dt, qvz, qlz, qrz, qiz, qsz, qgz,         &
                    qrzold, qszold, qgzold,                   &
                    thz, tothz, rhoz, orhoz, sqrhoz,          &
                    prez, zz, dzw, ht(i,j),                   &
                    pptrain, pptsnow, pptgraul, pptice,       &
                    grav, cp, rair, rvapor, gindex,           &
                    xls, xlv, xlf, rhowater, rhosnow,         &
                    ep2,svp1,svp2,svp3,svpt0,                 &
                    kts, kte, i, j                            )

!
! precipitation from cloud microphysics -- only for one time step
!
! unit is transferred from m to mm

!
   rain(i,j)=pptrain
   snow(i,j)=pptsnow
   graupel(i,j)=pptgraul
   ice(i,j)=pptice
!
   rainncv(i,j)= pptrain + pptsnow + pptgraul + pptice
   rainnc(i,j)=rainnc(i,j) + pptrain + pptsnow + pptgraul + pptice

!
!- update data from 1-d back to 3-d
!
!
   do k = kts, kte
      qv(i,k,j)=qvz(k)
      ql(i,k,j)=qlz(k)
      qr(i,k,j)=qrz(k)
      th(i,k,j)=thz(k)
   end do
!
   if (p_qi .ge. p_first_scalar) then
      do k = kts, kte
         qi(i,k,j)=qiz(k)
      enddo
   endif

   if (p_qs .ge. p_first_scalar) then
      do k = kts, kte
         qs(i,k,j)=qsz(k)
      enddo
   endif

   if (p_qg .ge. p_first_scalar) then
      do k = kts, kte
         qg(i,k,j)=qgz(k)
      enddo
   endif
!
   enddo i_loop
   enddo j_loop

   end subroutine lin_et_al


!-----------------------------------------------------------------------
   subroutine clphy1d(dt, qvz, qlz, qrz, qiz, qsz, qgz,                &
                      qrzold, qszold, qgzold,                          &
                      thz, tothz, rho, orho, sqrho,                    &
                      prez, zz, dzw, zsfc, pptrain, pptsnow, pptgraul, &
                      pptice, grav, cp, rair, rvapor, gindex,          &
                      xls, xlv, xlf, rhowater, rhosnow,                &
                      ep2,svp1,svp2,svp3,svpt0,                        &
                      kts, kte, i, j                                   )
!-----------------------------------------------------------------------
    implicit none
!-----------------------------------------------------------------------
!  this program handles the vertical 1-d cloud micphysics
!-----------------------------------------------------------------------
! avisc: constant in empirical formular for dynamic viscosity of air
!         =1.49628e-6 [kg/m/s] = 1.49628e-5 [g/cm/s]
! adiffwv: constant in empirical formular for diffusivity of water
!          vapor in air
!          = 8.7602e-5 [kgm/s3] = 8.7602 [gcm/s3]
! axka: constant in empirical formular for thermal conductivity of air
!       = 1.4132e3 [m2/s2/k] = 1.4132e7 [cm2/s2/k]
! qi0: mixing ratio threshold for cloud ice aggregation [kg/kg]
! xmi50: mass of a 50 micron ice crystal
!        = 4.8e-10 [kg] =4.8e-7 [g]
! xmi40: mass of a 40 micron ice crystal
!        = 2.46e-10 [kg] = 2.46e-7 [g]
! di50: diameter of a 50 micro (radius) ice crystal
!       =1.0e-4 [m]
! xmi: mass of one cloud ice crystal
!      =4.19e-13 [kg] = 4.19e-10 [g]
! oxmi=1.0/xmi
!
! xni0=1.0e-2 [m-3] the value given in lin et al. is wrong.(see
!                   hsie et al.(1980) and rutledge and hobbs(1983) )
! bni=0.5 [k-1]
! xmnin: mass of a natural ice nucleus
!    = 1.05e-18 [kg] = 1.05e-15 [g] this values is suggested by
!                    hsie et al. (1980)
!    = 1.0e-12 [kg] suggested by rutlegde and hobbs (1983)
! rhowater: density of water=1.0 g/cm3=1000.0 kg/m3
! consta: constant in empirical formular for terminal
!         velocity of raindrop
!         =2115.0 [cm**(1-b)/s] = 2115.0*0.01**(1-b) [m**(1-b)/s]
! constb: constant in empirical formular for terminal
!         velocity of raindrop
!         =0.8
! xnor: intercept parameter of the raindrop size distribution
!       = 0.08 cm-4 = 8.0e6 m-4
! araut: time sacle for autoconversion of cloud water to raindrops
!       =1.0e-3 [s-1]
! ql0: mixing ratio threshold for cloud watercoalescence [kg/kg]
! vf1r: ventilation factors for rain =0.78
! vf2r: ventilation factors for rain =0.31
! rhosnow: density of snow=0.1 g/cm3=100.0 kg/m3
! constc: constant in empirical formular for terminal
!         velocity of snow
!         =152.93 [cm**(1-d)/s] = 152.93*0.01**(1-d) [m**(1-d)/s]
! constd: constant in empirical formular for terminal
!         velocity of snow
!         =0.25
! xnos: intercept parameter of the snowflake size distribution
! vf1s: ventilation factors for snow =0.78
! vf2s: ventilation factors for snow =0.31
!
!----------------------------------------------------------------------

  integer, intent(in   )               :: kts, kte, i, j

  real,    dimension( kts:kte ),                                      &
           intent(inout)               :: qvz, qlz, qrz, qiz, qsz,    &
                                          qgz, thz

  real,    dimension( kts:kte ),                                      &
           intent(in   )               :: tothz, rho, orho, sqrho,    &
                                          prez, zz, qrzold, qszold,   &
                                          qgzold, dzw

  real,    intent(in   )               :: dt, grav, cp, rair, rvapor, &
                                          xls, xlv, xlf, rhowater,    &
                                          rhosnow,ep2,svp1,svp2,svp3,svpt0

  real,    intent(inout)               :: pptrain, pptsnow, pptgraul, pptice

  real,    intent(in   )               :: zsfc

! local vars

   real                                :: obp4, bp3, bp5, bp6, odp4,  &
                                          dp3, dp5, dp5o2


! temperary vars

   real                                :: tmp, tmp0, tmp1, tmp2,tmp3,  &
                                          tmp4,delta2,delta3, delta4,  &
                                          tmpa,tmpb,tmpc,tmpd,alpha1,  &
                                          qic, abi,abr, abg, odtberg,  &
                                          vti50,eiw,eri,esi,esr, esw,  &
                                          erw,delrs,term0,term1,araut, &
                                          constg2, vf1r, vf2r,alpha2,  &
                                          ap, bp, egw, egi, egs, egr,  &
                                          constg, gdelta4, g1sdelt4,   &
                                          factor, tmp_r, tmp_s,tmp_g,  &
                                          qlpqi, rsat, a1, a2, xnin

  integer                              :: k
!
  real, dimension( kts:kte )    ::  oprez, tem, temcc, theiz, qswz,    &
                                    qsiz, qvoqswz, qvoqsiz, qvzodt,    &
                                    qlzodt, qizodt, qszodt, qrzodt,    &
                                    qgzodt

  real, dimension( kts:kte )    :: psnow, psaut, psfw,  psfi,  praci,  &
                                   piacr, psaci, psacw, psdep, pssub,  &
                                   pracs, psacr, psmlt, psmltevp,      &
                                   prain, praut, pracw, prevp, pvapor, &
                                   pclw,  pladj, pcli,  pimlt, pihom,  &
                                   pidw,  piadj, pgraupel, pgaut,      &
                                   pgfr,  pgacw, pgaci, pgacr, pgacs,  &
                                   pgacip,pgacrp,pgacsp,pgwet, pdry,   &
                                   pgsub, pgdep, pgmlt, pgmltevp,      &
                                   qschg, qgchg
!

  real, dimension( kts:kte )    :: qvsbar, rs0, viscmu, visc, diffwv,  &
                                   schmidt, xka

  real, dimension( kts:kte )    :: vtr, vts, vtg,                      &
                                   vtrold, vtsold, vtgold, vtiold,     &
                                   xlambdar, xlambdas, xlambdag,       &
                                   olambdar, olambdas, olambdag

  real                          :: episp0k, dtb, odtb, pi, pio4,       &
                                   pio6, oxlf, xlvocp, xlfocp, consta, &
                                   constc, ocdrag, gambp4, gamdp4,     &
                                   gam4pt5, cpor, oxmi, gambp3, gamdp3,&
                                   gambp6, gam3pt5, gam2pt75, gambp5o2,&
                                   gamdp5o2, cwoxlf, ocp, xni50, es
!
  real                          :: qvmin=1.e-20
  real                          :: gindex
  real                          :: temc1,save1,save2,xni50mx

! for terminal velocity flux

  integer                       :: min_q, max_q
  real                          :: t_del_tv, del_tv, flux, fluxin, fluxout ,tmpqrz
  logical                       :: notlast
!

   dtb=dt
   odtb=1./dtb
   pi=acos(-1.)
   pio4=acos(-1.)/4.
   pio6=acos(-1.)/6.
   ocp=1./cp
   oxlf=1./xlf
   xlvocp=xlv/cp
   xlfocp=xlf/cp
   consta=2115.0*0.01**(1-constb)
   constc=152.93*0.01**(1-constd)
   ocdrag=1./cdrag
!  episp0k=rh*episp0
   episp0k=rh*ep2*1000.*svp1
!
   gambp4=ggamma(constb+4.)
   gamdp4=ggamma(constd+4.)
   gam4pt5=ggamma(4.5)
   cpor=cp/rair
   oxmi=1.0/xmi
   gambp3=ggamma(constb+3.)
   gamdp3=ggamma(constd+3.)
   gambp6=ggamma(constb+6)
   gam3pt5=ggamma(3.5)
   gam2pt75=ggamma(2.75)
   gambp5o2=ggamma((constb+5.)/2.)
   gamdp5o2=ggamma((constd+5.)/2.)
   cwoxlf=cw/xlf
   delta2=0.
   delta3=0.
   delta4=0.
!
!-----------------------------------------------------------------------
!     oprez       1./prez ( prez : pressure)
!     qsw         saturated mixing ratio on water surface
!     qsi         saturated mixing ratio on ice surface
!     episp0k     rh*e*saturated pressure at 273.15 k
!     qvoqsw      qv/qsw
!     qvoqsi      qv/qsi
!     qvzodt      qv/dt
!     qlzodt      ql/dt
!     qizodt      qi/dt
!     qszodt      qs/dt
!     qrzodt      qr/dt
!     qgzodt      qg/dt
!
!     temcc       temperature in dregee c
!

      obp4=1.0/(constb+4.0)
      bp3=constb+3.0
      bp5=constb+5.0
      bp6=constb+6.0
      odp4=1.0/(constd+4.0)
      dp3=constd+3.0
      dp5=constd+5.0
      dp5o2=0.5*(constd+5.0)
!
      do k=kts,kte
         oprez(k)=1./prez(k)
      enddo

      do k=kts,kte
         qlz(k)=amax1( 0.0,qlz(k) )
         qiz(k)=amax1( 0.0,qiz(k) )
         qvz(k)=amax1( qvmin,qvz(k) )
         qsz(k)=amax1( 0.0,qsz(k) )
         qrz(k)=amax1( 0.0,qrz(k) )
         qgz(k)=amax1( 0.0,qgz(k) )
!
         tem(k)=thz(k)*tothz(k)
         temcc(k)=tem(k)-273.15
!
!        qswz(k)=episp0k*oprez(k)* &
!               exp( svp2*temcc(k)/(tem(k)-svp3) )
         es=1000.*svp1*exp( svp2*temcc(k)/(tem(k)-svp3) )
         qswz(k)=ep2*es/(prez(k)-es)
         if (tem(k) .lt. 273.15 ) then
!           qsiz(k)=episp0k*oprez(k)* &
!                  exp( 21.8745584*(tem(k)-273.16)/(tem(k)-7.66) )
            es=1000.*svp1*exp( 21.8745584*(tem(k)-273.16)/(tem(k)-7.66) )
            qsiz(k)=ep2*es/(prez(k)-es)
            if (temcc(k) .lt. -40.0) qswz(k)=qsiz(k)
         else
            qsiz(k)=qswz(k)
         endif
!
         qvoqswz(k)=qvz(k)/qswz(k)
         qvoqsiz(k)=qvz(k)/qsiz(k)
         qvzodt(k)=amax1( 0.0,odtb*qvz(k) )
         qlzodt(k)=amax1( 0.0,odtb*qlz(k) )
         qizodt(k)=amax1( 0.0,odtb*qiz(k) )
         qszodt(k)=amax1( 0.0,odtb*qsz(k) )
         qrzodt(k)=amax1( 0.0,odtb*qrz(k) )
         qgzodt(k)=amax1( 0.0,odtb*qgz(k) )

         theiz(k)=thz(k)+(xlvocp*qvz(k)-xlfocp*qiz(k))/tothz(k)
      enddo


!
!
!-----------------------------------------------------------------------
! in this simple stable cloud parameterization scheme, only five
! forms of water substance (water vapor, cloud water, cloud ice,
! rain and snow are considered. the prognostic variables are total
! water (qp),cloud water (ql), and cloud ice (qi). rain and snow are
! diagnosed following nagata and ogura, 1991, mwr, 1309-1337. eq (a7).
! the micro physics are based on (1) hsie et al.,1980, jam, 950-977 ;
! (2) lin et al., 1983, jam, 1065-1092 ; (3) rutledge and hobbs, 1983,
! jas, 1185-1206 ; (4) rutledge and hobbs, 1984, jas, 2949-2972.
!-----------------------------------------------------------------------
!
! rhowater: density of water=1.0 g/cm3=1000.0 kg/m3
! rhosnow: density of snow=0.1 g/cm3=100.0 kg/m3
! xnor: intercept parameter of the raindrop size distribution
!       = 0.08 cm-4 = 8.0e6 m-4
! xnos: intercept parameter of the snowflake size distribution
!       = 0.03 cm-4 = 3.0e6 m-4
! xnog: intercept parameter of the graupel size distribution
!       = 4.0e-4 cm-4 = 4.0e4 m-4
! consta: constant in empirical formular for terminal
!         velocity of raindrop
!         =2115.0 [cm**(1-b)/s] = 2115.0*0.01**(1-b) [m**(1-b)/s]
! constb: constant in empirical formular for terminal
!         velocity of raindrop
!         =0.8
! constc: constant in empirical formular for terminal
!         velocity of snow
!         =152.93 [cm**(1-d)/s] = 152.93*0.01**(1-d) [m**(1-d)/s]
! constd: constant in empirical formular for terminal
!         velocity of snow
!         =0.25
! avisc: constant in empirical formular for dynamic viscosity of air
!         =1.49628e-6 [kg/m/s] = 1.49628e-5 [g/cm/s]
! adiffwv: constant in empirical formular for diffusivity of water
!          vapor in air
!          = 8.7602e-5 [kgm/s3] = 8.7602 [gcm/s3]
! axka: constant in empirical formular for thermal conductivity of air
!       = 1.4132e3 [m2/s2/k] = 1.4132e7 [cm2/s2/k]
! qi0: mixing ratio threshold for cloud ice aggregation [kg/kg]
!      = 1.0e-3 g/g = 1.0e-3 kg/gk
! ql0: mixing ratio threshold for cloud watercoalescence [kg/kg]
!      = 2.0e-3 g/g = 2.0e-3 kg/gk
! qs0: mixing ratio threshold for snow aggregation
!      = 6.0e-4 g/g = 6.0e-4 kg/gk
! xmi50: mass of a 50 micron ice crystal
!        = 4.8e-10 [kg] =4.8e-7 [g]
! xmi40: mass of a 40 micron ice crystal
!        = 2.46e-10 [kg] = 2.46e-7 [g]
! di50: diameter of a 50 micro (radius) ice crystal
!       =1.0e-4 [m]
! xmi: mass of one cloud ice crystal
!      =4.19e-13 [kg] = 4.19e-10 [g]
! oxmi=1.0/xmi
!


! if gindex=1.0 include graupel
! if gindex=0. no graupel
!
!
      do k=kts,kte
         psnow(k)=0.0
         psaut(k)=0.0
         psfw(k)=0.0
         psfi(k)=0.0
         praci(k)=0.0
         piacr(k)=0.0
         psaci(k)=0.0
         psacw(k)=0.0
         psdep(k)=0.0
         pssub(k)=0.0
         pracs(k)=0.0
         psacr(k)=0.0
         psmlt(k)=0.0
         psmltevp(k)=0.0
!
         prain(k)=0.0
         praut(k)=0.0
         pracw(k)=0.0
         prevp(k)=0.0
!
         pvapor(k)=0.0
!
         pclw(k)=0.0
         pladj(k)=0.0
!
         pcli(k)=0.0
         pimlt(k)=0.0
         pihom(k)=0.0
         pidw(k)=0.0
         piadj(k)=0.0
      enddo

!
!!! graupel
!
      do k=kts,kte
         pgraupel(k)=0.0
         pgaut(k)=0.0
         pgfr(k)=0.0
         pgacw(k)=0.0
         pgaci(k)=0.0
         pgacr(k)=0.0
         pgacs(k)=0.0
         pgacip(k)=0.0
         pgacrp(k)=0.0
         pgacsp(k)=0.0
         pgwet(k)=0.0
         pdry(k)=0.0
         pgsub(k)=0.0
         pgdep(k)=0.0
         pgmlt(k)=0.0
         pgmltevp(k)=0.0
         qschg(k)=0.
         qgchg(k)=0.
      enddo
!
!
! set rs0=episp0*oprez(k)
! episp0=e*saturated pressure at 273.15 k
! e     = 0.622
!
      do k=kts,kte
         rs0(k)=ep2*1000.*svp1/(prez(k)-1000.*svp1)
      end do
!
!***********************************************************************
! calculate precipitation fluxes due to terminal velocities.
!***********************************************************************
!
!- calculate termianl velocity (vt?)  of precipitation q?z
!- find maximum vt? to determine the small delta t
!
!-- rain
!
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.
    do while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
!        if (qrzold(k) .gt. 1.0e-12) then
!!       if (qrzold(k) .gt. 1.0e-8) then
         if (qrz(k) .gt. 1.0e-8) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
!!          tmp1=sqrt(pi*rhowater*xnor/rho(k)/qrzold(k))
            tmp1=sqrt(pi*rhowater*xnor/rho(k)/qrz(k))
            tmp1=sqrt(tmp1)
            vtrold(k)=o6*consta*gambp4*sqrho(k)/tmp1**constb
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-zsfc)/vtrold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtrold(k))
            endif
         else
            vtrold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!- check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv
!
         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the q?z flux
! termi is the qrzold flux pass in the grid box through the upper boundary
! termo is the qrzold flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
!           fluxout=rho(k)*vtrold(k)*qrzold(k)
            fluxout=rho(k)*vtrold(k)*qrz(k)
            flux=(fluxin-fluxout)/rho(k)/dzw(k)
            tmpqrz=qrz(k)
            qrz(k)=qrz(k)+del_tv*flux
            fluxin=fluxout
         enddo
         if (min_q .eq. 1) then
            pptrain=pptrain+fluxin*del_tv
         else
            qrz(min_q-1)=qrz(min_q-1)+del_tv*  &
                          fluxin/rho(min_q-1)/dzw(min_q-1)
         endif
!
      else
         notlast=.false.
      endif
    enddo

!
!-- snow
!
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.

    do while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
!        if (qszold(k) .gt. 1.0e-12) then
!!       if (qszold(k) .gt. 1.0e-8) then
         if (qsz(k) .gt. 1.0e-8) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
!!          tmp1=sqrt(pi*rhosnow*xnos/rho(k)/qszold(k))
            tmp1=sqrt(pi*rhosnow*xnos/rho(k)/qsz(k))
            tmp1=sqrt(tmp1)
            vtsold(k)=o6*constc*gamdp4*sqrho(k)/tmp1**constd
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-zsfc)/vtsold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtsold(k))
            endif
         else
            vtsold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qszold flux
! termi is the qszold flux pass in the grid box through the upper boundary
! termo is the qszold flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
!!          fluxout=rho(k)*vtsold(k)*qszold(k)
            fluxout=rho(k)*vtsold(k)*qsz(k)
            flux=(fluxin-fluxout)/rho(k)/dzw(k)
            qsz(k)=qsz(k)+del_tv*flux
            qsz(k)=amax1(0.,qsz(k))
            fluxin=fluxout
         enddo
         if (min_q .eq. 1) then
            pptsnow=pptsnow+fluxin*del_tv
         else
            qsz(min_q-1)=qsz(min_q-1)+del_tv*  &
                         fluxin/rho(min_q-1)/dzw(min_q-1)
         endif
!
      else
         notlast=.false.
      endif

    enddo
!
!-- grauupel
!
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.
!
    do while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
!        if (qgzold(k) .gt. 1.0e-12) then
!!       if (qgzold(k) .gt. 1.0e-8) then
         if (qgz(k) .gt. 1.0e-8) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
!!          tmp1=sqrt(pi*rhograul*xnog/rho(k)/qgzold(k))
            tmp1=sqrt(pi*rhograul*xnog/rho(k)/qgz(k))
            tmp1=sqrt(tmp1)
            term0=sqrt(4.*grav*rhograul*0.33334/rho(k)/cdrag)
            vtgold(k)=o6*gam4pt5*term0*sqrt(1./tmp1)
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-zsfc)/vtgold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtgold(k))
            endif
         else
            vtgold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qgzold flux
! termi is the qgzold flux pass in the grid box through the upper boundary
! termo is the qgzold flux pass out the grid box through the lower boundary
!
         fluxin=0.
         do k=max_q,min_q,-1
!!          fluxout=rho(k)*vtgold(k)*qgzold(k)
            fluxout=rho(k)*vtgold(k)*qgz(k)
            flux=(fluxin-fluxout)/rho(k)/dzw(k)
            qgz(k)=qgz(k)+del_tv*flux
            qgz(k)=amax1(0.,qgz(k))
            fluxin=fluxout
         enddo
         if (min_q .eq. 1) then
            pptgraul=pptgraul+fluxin*del_tv
         else
            qgz(min_q-1)=qgz(min_q-1)+del_tv*  &
                         fluxin/rho(min_q-1)/dzw(min_q-1)
         endif
!
      else
         notlast=.false.
      endif
!
   enddo

!
!-- cloud ice  (03/21/02) follow vaughan t.j. phillips at gfdl
!
    t_del_tv=0.
    del_tv=dtb
    notlast=.true.
!
    do while (notlast)
!
      min_q=kte
      max_q=kts-1
!
      do k=kts,kte-1
         if (qiz(k) .gt. 1.0e-8) then
            min_q=min0(min_q,k)
            max_q=max0(max_q,k)
            vtiold(k)= 3.29 * (rho(k)* qiz(k))** 0.16  ! heymsfield and donner
            if (k .eq. 1) then
               del_tv=amin1(del_tv,0.9*(zz(k)-zsfc)/vtiold(k))
            else
               del_tv=amin1(del_tv,0.9*(zz(k)-zz(k-1))/vtiold(k))
            endif
         else
            vtiold(k)=0.
         endif
      enddo

      if (max_q .ge. min_q) then
!
!
!- check if the summation of the small delta t >=  big delta t
!             (t_del_tv)          (del_tv)             (dtb)

         t_del_tv=t_del_tv+del_tv

         if ( t_del_tv .ge. dtb ) then
              notlast=.false.
              del_tv=dtb+del_tv-t_del_tv
         endif

! use small delta t to calculate the qgzold flux
! termi is the qgzold flux pass in the grid box through the upper boundary
! termo is the qgzold flux pass out the grid box through the lower boundary
!

         fluxin=0.
         do k=max_q,min_q,-1
            fluxout=rho(k)*vtiold(k)*qiz(k)
            flux=(fluxin-fluxout)/rho(k)/dzw(k)
            qiz(k)=qiz(k)+del_tv*flux
            qiz(k)=amax1(0.,qiz(k))
            fluxin=fluxout
         enddo
         if (min_q .eq. 1) then
            pptice=pptice+fluxin*del_tv
         else
            qiz(min_q-1)=qiz(min_q-1)+del_tv*  &
                         fluxin/rho(min_q-1)/dzw(min_q-1)
         endif
!
      else
         notlast=.false.
      endif
!
   enddo

! microphpysics processes
!
      do 2000 k=kts,kte
!
!***********************************************************************
!*****   diagnose mixing ratios (qrz,qsz), terminal                *****
!*****   velocities (vtr,vts), and slope parameters in size        *****
!*****   distribution(xlambdar,xlambdas) of rain and snow          *****
!*****   follows nagata and ogura, 1991, mwr, 1309-1337. eq (a7)   *****
!***********************************************************************
!
!**** assuming no cloud water can exist in the top two levels due to
!**** radiation consideration
!
!!  if
!!     unsaturated,
!!     no cloud water, rain, ice, snow and graupel
!!  then
!!     skip these processes and jump to line 2000
!
!
        tmp=qiz(k)+qlz(k)+qsz(k)+qrz(k)+qgz(k)*gindex
        if( qvz(k)+qlz(k)+qiz(k) .lt. qsiz(k)  &
            .and. tmp .eq. 0.0 ) go to 2000

!! calculate terminal velocity of rain
!
!       if (qrzold(k) .gt. 1.0e-12) then
!!      if (qrzold(k) .gt. 1.0e-8) then
        if (qrz(k) .gt. 1.0e-8) then
!!          tmp1=sqrt(pi*rhowater*xnor*orho(k)/qrzold(k))
            tmp1=sqrt(pi*rhowater*xnor*orho(k)/qrz(k))
            xlambdar(k)=sqrt(tmp1)
            olambdar(k)=1.0/xlambdar(k)
            vtrold(k)=o6*consta*gambp4*sqrho(k)*olambdar(k)**constb
        else
            vtrold(k)=0.
            olambdar(k)=0.
        endif
!
!       if (qrz(k) .gt. 1.0e-12) then
        if (qrz(k) .gt. 1.0e-8) then
            tmp1=sqrt(pi*rhowater*xnor*orho(k)/qrz(k))
            xlambdar(k)=sqrt(tmp1)
            olambdar(k)=1.0/xlambdar(k)
            vtr(k)=o6*consta*gambp4*sqrho(k)*olambdar(k)**constb
        else
            vtr(k)=0.
            olambdar(k)=0.
        endif
!
!! calculate terminal velocity of snow
!
!       if (qszold(k) .gt. 1.0e-12) then
!!      if (qszold(k) .gt. 1.0e-8) then
        if (qsz(k) .gt. 1.0e-8) then
!!          tmp1=sqrt(pi*rhosnow*xnos*orho(k)/qszold(k))
            tmp1=sqrt(pi*rhosnow*xnos*orho(k)/qsz(k))
            xlambdas(k)=sqrt(tmp1)
            olambdas(k)=1.0/xlambdas(k)
            vtsold(k)=o6*constc*gamdp4*sqrho(k)*olambdas(k)**constd
        else
            vtsold(k)=0.
            olambdas(k)=0.
        endif
!
!       if (qsz(k) .gt. 1.0e-12) then
        if (qsz(k) .gt. 1.0e-8) then
            tmp1=sqrt(pi*rhosnow*xnos*orho(k)/qsz(k))
            xlambdas(k)=sqrt(tmp1)
            olambdas(k)=1.0/xlambdas(k)
            vts(k)=o6*constc*gamdp4*sqrho(k)*olambdas(k)**constd
        else
            vts(k)=0.
            olambdas(k)=0.
        endif
!
!! calculate terminal velocity of graupel
!
!       if (qgzold(k) .gt. 1.0e-12) then
!!      if (qgzold(k) .gt. 1.0e-8) then
        if (qgz(k) .gt. 1.0e-8) then
!!          tmp1=sqrt( pi*rhograul*xnog*orho(k)/qgzold(k))
            tmp1=sqrt( pi*rhograul*xnog*orho(k)/qgz(k))
            xlambdag(k)=sqrt(tmp1)
            olambdag(k)=1.0/xlambdag(k)
            term0=sqrt(4.*grav*rhograul*0.33334*orho(k)*ocdrag)
            vtgold(k)=o6*gam4pt5*term0*sqrt(olambdag(k))
        else
            vtgold(k)=0.
            olambdag(k)=0.
        endif
!
!       if (qgz(k) .gt. 1.0e-12) then
        if (qgz(k) .gt. 1.0e-8) then
            tmp1=sqrt( pi*rhograul*xnog*orho(k)/qgz(k))
            xlambdag(k)=sqrt(tmp1)
            olambdag(k)=1.0/xlambdag(k)
            term0=sqrt(4.*grav*rhograul*0.33334*orho(k)*ocdrag)
            vtg(k)=o6*gam4pt5*term0*sqrt(olambdag(k))
        else
            vtg(k)=0.
            olambdag(k)=0.
        endif
!
!***********************************************************************
!*****  compute viscosity,difusivity,thermal conductivity, and    ******
!*****  schmidt number                                            ******
!***********************************************************************
!c------------------------------------------------------------------
!c      viscmu: dynamic viscosity of air kg/m/s
!c      visc: kinematic viscosity of air = viscmu/rho (m2/s)
!c      avisc=1.49628e-6 kg/m/s=1.49628e-5 g/cm/s
!c      viscmu=1.718e-5 kg/m/s in rh
!c      diffwv: diffusivity of water vapor in air
!c      adiffwv = 8.7602e-5 (8.794e-5 in mm5) kgm/s3
!c              = 8.7602 (8.794 in mm5) gcm/s3
!c      diffwv(k)=2.26e-5 m2/s
!c      schmidt: schmidt number=visc/diffwv
!c      xka: thermal conductivity of air j/m/s/k (kgm/s3/k)
!c      xka(k)=2.43e-2 j/m/s/k in rh
!c      axka=1.4132e3 (1.414e3 in mm5) m2/s2/k = 1.4132e7 cm2/s2/k
!c------------------------------------------------------------------

        viscmu(k)=avisc*tem(k)**1.5/(tem(k)+120.0)
        visc(k)=viscmu(k)*orho(k)
        diffwv(k)=adiffwv*tem(k)**1.81*oprez(k)
        schmidt(k)=visc(k)/diffwv(k)
        xka(k)=axka*viscmu(k)

        if (tem(k) .lt. 273.15) then

!
!***********************************************************************
!*********        snow production processes for t < 0 c       **********
!***********************************************************************
!c
!c (1) ice crystal aggregation to snow (psaut): lin (21)
!c!    psaut=alpha1*(qi-qi0)
!c!    alpha1=1.0e-3*exp(0.025*(t-t0))
!c
!          alpha1=1.0e-3*exp( 0.025*temcc(k) )

           alpha1=1.0e-3*exp( 0.025*temcc(k) )
!
           if(temcc(k) .lt. -20.0) then
             tmp1=-7.6+4.0*exp( -0.2443e-3*(abs(temcc(k))-20)**2.455 )
             qic=1.0e-3*exp(tmp1)*orho(k)
           else
             qic=qi0
           end if
!testing
!          tmp1=amax1( 0.0,alpha1*(qiz(k)-qic) )
!          psaut(k)=amin1( tmp1,qizodt(k) )

           tmp1=odtb*(qiz(k)-qic)*(1.0-exp(-alpha1*dtb))
           psaut(k)=amax1( 0.0,tmp1 )

!c
!c (2) bergeron process transfer of cloud water to snow (psfw)
!c     this process only considered when -31 c < t < 0 c
!c     lin (33) and hsie (17)
!c
!c!
!c!    parama1 and parama2 functions must be user supplied
!c!

! testing
          if( qlz(k) .gt. 1.0e-10 ) then
            temc1=amax1(-30.99,temcc(k))
!           print*,'temc1',temc1,qlz(k)
            a1=parama1( temc1 )
            a2=parama2( temc1 )
            tmp1=1.0-a2
!! change unit from cgs to mks
            a1=a1*0.001**tmp1
!c!   dtberg is the time needed for a crystal to grow from 40 to 50 um
!c !  odtberg=1.0/dtberg
            odtberg=(a1*tmp1)/(xmi50**tmp1-xmi40**tmp1)
!
!c!   compute terminal velocity of a 50 micron ice cystal
!
            vti50=constc*di50**constd*sqrho(k)
!
            eiw=1.0
            save1=a1*xmi50**a2
            save2=0.25*pi*eiw*rho(k)*di50*di50*vti50
!
            tmp2=( save1 + save2*qlz(k) )
!
!! maximum number of 50 micron crystals limited by the amount
!!  of supercool water
!
            xni50mx=qlzodt(k)/tmp2
!
!!   number of 50 micron crystals produced
!
!
            xni50=qiz(k)*( 1.0-exp(-dtb*odtberg) )/xmi50
            xni50=amin1(xni50,xni50mx)
!
            tmp3=odtb*tmp2/save2*( 1.0-exp(-save2*xni50*dtb) )
            psfw(k)=amin1( tmp3,qlzodt(k) )
!testing
!           psfw(k)=0.

!0915     if( temcc(k).gt.-30.99 ) then
!0915        a1=parama1( temcc(k) )
!0915        a2=parama2( temcc(k) )
!0915        tmp1=1.0-a2
!!     change unit from cgs to mks
!0915        a1=a1*0.001**tmp1

!c!    dtberg is the time needed for a crystal to grow from 40 to 50 um
!c!    odtberg=1.0/dtberg
!0915        odtberg=(a1*tmp1)/(xmi50**tmp1-xmi40**tmp1)

!c!    number of 50 micron crystals produced
!0915        xni50=qiz(k)*dtb*odtberg/xmi50

!c!    need to calculate the terminal velocity of a 50 micron
!c!    ice cystal
!0915        vti50=constc*di50**constd*sqrho(k)
!0915        eiw=1.0
!0915        tmp2=xni50*( a1*xmi50**a2 + &
!0915             0.25*qlz(k)*pi*eiw*rho(k)*di50*di50*vti50 )
!0915        psfw(k)=amin1( tmp2,qlzodt(k) )
!0915        psfw(k)=0.
!c
!c (3) reduction of cloud ice by bergeron process (psfi): lin (34)
!c     this process only considered when -31 c < t < 0 c
!c
            tmp1=xni50*xmi50-psfw(k)
            psfi(k)=amin1(tmp1,qizodt(k))
! testing
!           psfi(k)=0.
          end if
!

!0915        tmp1=qiz(k)*odtberg
!0915        psfi(k)=amin1(tmp1,qizodt(k))
! testing
!0915        psfi(k)=0.
!0915     end if
!
          if(qrz(k) .le. 0.0) go to 1000
!
! processes (4) and (5) only need when qrz > 0.0
!
!c
!c (4) cloud ice accretion by rain (praci): lin (25)
!c     may produce snow or graupel
!c
          eri=1.0
!0915     tmp1=qiz(k)*pio4*eri*xnor*consta*sqrho(k)
!0915     tmp2=tmp1*gambp3*olambdar(k)**bp3
!0915     praci(k)=amin1( tmp2,qizodt(k) )

          save1=pio4*eri*xnor*consta*sqrho(k)
          tmp1=save1*gambp3*olambdar(k)**bp3
          praci(k)=qizodt(k)*( 1.0-exp(-tmp1*dtb) )

!c
!c (5) rain accretion by cloud ice (piacr): lin (26)
!c
!0915     tmp2=tmp1*rho(k)*pio6*rhowater*gambp6*oxmi* &
!0915              olambdar(k)**bp6
!0915     piacr(k)=amin1( tmp2,qrzodt(k) )

          tmp2=qiz(k)*save1*rho(k)*pio6*rhowater*gambp6*oxmi* &
                   olambdar(k)**bp6
          piacr(k)=amin1( tmp2,qrzodt(k) )

!
1000      continue
!
          if(qsz(k) .le. 0.0) go to 1200
!
! compute the following processes only when qsz > 0.0
!
!c
!c (6) ice crystal accretion by snow (psaci): lin (22)
!c
          esi=exp( 0.025*temcc(k) )
          save1=pio4*xnos*constc*gamdp3*sqrho(k)* &
               olambdas(k)**dp3
          tmp1=esi*save1
          psaci(k)=qizodt(k)*( 1.0-exp(-tmp1*dtb) )

!0915     tmp1=pio4*xnos*constc*gamdp3*sqrho(k)* &
!0915          olambdas(k)**dp3
!0915     tmp2=qiz(k)*esi*tmp1
!0915     psaci(k)=amin1( tmp2,qizodt(k) )
!c
!c (7) cloud water accretion by snow (psacw): lin (24)
!c
          esw=1.0
          tmp1=esw*save1
          psacw(k)=qlzodt(k)*( 1.0-exp(-tmp1*dtb) )

!0915     tmp2=qlz(k)*esw*tmp1
!0915     psacw(k)=amin1( tmp2,qlzodt(k) )
!c
!c (8) deposition/sublimation of snow (psdep/pssub): lin (31)
!c     includes consideration of ventilation effect
!c
!c     abi=2*pi*(si-1)/rho/(a"+b")
!c
          tmpa=rvapor*xka(k)*tem(k)*tem(k)
          tmpb=xls*xls*rho(k)*qsiz(k)*diffwv(k)
          tmpc=tmpa*qsiz(k)*diffwv(k)
          abi=2.0*pi*(qvoqsiz(k)-1.0)*tmpc/(tmpa+tmpb)
!
!c     vf1s,vf2s=ventilation factors for snow
!c     vf1s=0.78,vf2s=0.31 in lin
!
          tmp1=constc*sqrho(k)*olambdas(k)**dp5/visc(k)
          tmp2=abi*xnos*( vf1s*olambdas(k)*olambdas(k)+ &
                    vf2s*schmidt(k)**0.33334*gamdp5o2*sqrt(tmp1) )
          tmp3=odtb*( qvz(k)-qsiz(k) )
!
          if( tmp2 .le. 0.0) then
            tmp2=amax1( tmp2,tmp3)
            pssub(k)=amax1( tmp2,-qszodt(k) )
            psdep(k)=0.0
          else
            psdep(k)=amin1( tmp2,tmp3 )
            pssub(k)=0.0
          end if

!0915     psdep(k)=amax1(0.0,tmp2)
!0915     pssub(k)=amin1(0.0,tmp2)
!0915     pssub(k)=amax1( pssub(k),-qszodt(k) )
!
          if(qrz(k) .le. 0.0) go to 1200
!
! compute processes (9) and (10) only when qsz > 0.0 and qrz > 0.0
!
!c
!c (9) accretion of snow by rain (pracs): lin (27)
!c
          esr=1.0
          tmpa=olambdar(k)*olambdar(k)
          tmpb=olambdas(k)*olambdas(k)
          tmpc=olambdar(k)*olambdas(k)
          tmp1=pi*pi*esr*xnor*xnos*abs( vtr(k)-vts(k) )*orho(k)
          tmp2=tmpb*tmpb*olambdar(k)*(5.0*tmpb+2.0*tmpc+0.5*tmpa)
          tmp3=tmp1*rhosnow*tmp2
          pracs(k)=amin1( tmp3,qszodt(k) )
!c
!c (10)  accretion of rain by snow (psacr): lin (28)
!c
          tmp3=tmpa*tmpa*olambdas(k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
          tmp4=tmp1*rhowater*tmp3
          psacr(k)=amin1( tmp4,qrzodt(k) )
!
1200      continue
!
        else
!
!***********************************************************************
!*********        snow production processes for t > 0 c       **********
!***********************************************************************
!
         if (qsz(k) .le. 0.0) go to 1400
!c
!c (1) cloud water accretion by snow (psacw): lin (24)
!c
            esw=1.0

            tmp1=esw*pio4*xnos*constc*gamdp3*sqrho(k)* &
                 olambdas(k)**dp3
            psacw(k)=qlzodt(k)*( 1.0-exp(-tmp1*dtb) )

!0915       tmp1=pio4*xnos*constc*gamdp3*sqrho(k)* &
!0915            olambdas(k)**dp3
!0915       tmp2=qlz(k)*esw*tmp1
!0915       psacw(k)=amin1( tmp2,qlzodt(k) )
!c
!c (2) accretion of rain by snow (psacr): lin (28)
!c
            esr=1.0
            tmpa=olambdar(k)*olambdar(k)
            tmpb=olambdas(k)*olambdas(k)
            tmpc=olambdar(k)*olambdas(k)
            tmp1=pi*pi*esr*xnor*xnos*abs( vtr(k)-vts(k) )*orho(k)
            tmp2=tmpa*tmpa*olambdas(k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
            tmp3=tmp1*rhowater*tmp2
            psacr(k)=amin1( tmp3,qrzodt(k) )
!c
!c (3) melting of snow (psmlt): lin (32)
!c     psmlt is negative value
!
            delrs=rs0(k)-qvz(k)
            term1=2.0*pi*orho(k)*( xlv*diffwv(k)*rho(k)*delrs- &
                  xka(k)*temcc(k) )
            tmp1=constc*sqrho(k)*olambdas(k)**dp5/visc(k)
            tmp2=xnos*( vf1s*olambdas(k)*olambdas(k)+  &
                 vf2s*schmidt(k)**0.33334*gamdp5o2*sqrt(tmp1) )
            tmp3=term1*oxlf*tmp2-cwoxlf*temcc(k)*( psacw(k)+psacr(k) )
            tmp4=amin1(0.0,tmp3)
            psmlt(k)=amax1( tmp4,-qszodt(k) )
!c
!c (4) evaporation of melting snow (psmltevp): hr (a27)
!c     but use lin et al. coefficience
!c     psmltevp is a negative value
!c
            tmpa=rvapor*xka(k)*tem(k)*tem(k)
            tmpb=xlv*xlv*rho(k)*qswz(k)*diffwv(k)
            tmpc=tmpa*qswz(k)*diffwv(k)
            tmpd=amin1( 0.0,(qvoqswz(k)-0.90)*qswz(k)*odtb )

!      abr=2.0*pi*(qvoqswz(k)-1.0)*tmpc/(tmpa+tmpb)

            abr=2.0*pi*(qvoqswz(k)-0.90)*tmpc/(tmpa+tmpb)
!
!**** allow evaporation to occur when rh less than 90%
!**** here not using 100% because the evaporation cooling
!**** of temperature is not taking into account yet; hence,
!**** the qsw value is a little bit larger. this will avoid
!**** evaporation can generate cloud.
!
!c    vf1s,vf2s=ventilation factors for snow
!c    vf1s=0.78,vf2s=0.31 in lin
!
            tmp1=constc*sqrho(k)*olambdas(k)**dp5/visc(k)
            tmp2=abr*xnos*( vf1s*olambdas(k)*olambdas(k)+  &
                 vf2s*schmidt(k)**0.33334*gamdp5o2*sqrt(tmp1) )
            tmp3=amin1(0.0,tmp2)
            tmp3=amax1( tmp3,tmpd )
            psmltevp(k)=amax1( tmp3,-qszodt(k) )
1400     continue
!
        end if

!***********************************************************************
!*********           rain production processes                **********
!***********************************************************************
!
!c
!c (1) autoconversion of rain (praut): rh
!c     araut=afa*rho
!c     afa=0.001 rate coefficient for autoconvergence
!c
!c     araut=1.0e-3
!c
         araut=0.001
!testing
!        tmp1=amax1( 0.0,araut*(qlz(k)-ql0) )
!        praut(k)=amin1( tmp1,qlzodt(k) )
         tmp1=odtb*(qlz(k)-ql0)*( 1.0-exp(-araut*dtb) )
         praut(k)=amax1( 0.0,tmp1 )

!c
!c (2) accretion of cloud water by rain (pracw): lin (51)
!c
         erw=1.0
!        tmp1=qlz(k)*pio4*erw*xnor*consta*sqrho(k)
!        tmp2=tmp1*gambp3*olambdar(k)**bp3
!        pracw(k)=amin1( tmp2,qlzodt(k) )

        tmp1=pio4*erw*xnor*consta*sqrho(k)* &
             gambp3*olambdar(k)**bp3
        pracw(k)=qlzodt(k)*( 1.0-exp(-tmp1*dtb) )

!c
!c (3) evaporation of rain (prevp): lin (52)
!c     prevp is negative value
!c
!c     sw=qvoqsw : saturation ratio
!c
         tmpa=rvapor*xka(k)*tem(k)*tem(k)
         tmpb=xlv*xlv*rho(k)*qswz(k)*diffwv(k)
         tmpc=tmpa*qswz(k)*diffwv(k)
         tmpd=amin1(0.0,(qvoqswz(k)-0.90)*qswz(k)*odtb)
!
!      abr=2.0*pi*(qvoqswz(k)-1.0)*tmpc/(tmpa+tmpb)

         abr=2.0*pi*(qvoqswz(k)-0.90)*tmpc/(tmpa+tmpb)
!
!c     vf1r,vf2r=ventilation factors for rain
!c     vf1r=0.78,vf2r=0.31 in rh, lin  and mm5
!
         vf1r=0.78
         vf2r=0.31
         tmp1=consta*sqrho(k)*olambdar(k)**bp5/visc(k)
         tmp2=abr*xnor*( vf1r*olambdar(k)*olambdar(k)+  &
              vf2r*schmidt(k)**0.33334*gambp5o2*sqrt(tmp1) )
         tmp3=amin1( 0.0,tmp2 )
         tmp3=amax1( tmp3,tmpd )
         prevp(k)=amax1( tmp3,-qrzodt(k) )

!
!      if(iout .gt. 0) write(20,*)'tmp1,tmp2,tmp3=',tmp1,tmp2,tmp3
!      if(iout .gt. 0) write(20,*)'qlz,qiz,qrz=',qlz(k),qiz(k),qrz(k)
!      if(iout .gt. 0) write(20,*)'tem,qsz,qvz=',tem(k),qsz(k),qvz(k)



!     if (gindex .eq. 0.) goto 900
!
      if (tem(k) .lt. 273.15) then
!
!
!-- graupel
!***********************************************************************
!*********        graupel production processes for t < 0 c    **********
!***********************************************************************
!c
!c (1) autoconversion of snow to form graupel (pgaut): lin (37)
!c     pgaut=alpha2*(qsz-qs0)
!c     qs0=6.0e-4
!c     alpha2=1.0e-3*exp(0.09*temcc(k))      lin (38)
!
            alpha2=1.0e-3*exp(0.09*temcc(k))
!

! testing
!           tmp1=alpha2*(qsz(k)-qs0)
!           tmp1=amax1(0.0,tmp1)
!           pgaut(k)=amin1( tmp1,qszodt(k) )

            tmp1=odtb*(qsz(k)-qs0)*(1.0-exp(-alpha2*dtb))
            pgaut(k)=amax1( 0.0,tmp1 )

!c
!c (2) freezing of rain to form graupel  (pgfr): lin (45)
!c     positive value
!c     constant in bigg freezing aplume=ap=0.66 /k
!c     constant in raindrop freezing equ. bplume=bp=100./m/m/m/s
!

            if (qrz(k) .gt. 1.e-8 ) then
               bp=100.
               ap=0.66
               tmp1=olambdar(k)*olambdar(k)*olambdar(k)
               tmp2=20.*pi*pi*bp*xnor*rhowater*orho(k)*  &
                    (exp(-ap*temcc(k))-1.0)*tmp1*tmp1*olambdar(k)
               pgfr(k)=amin1( tmp2,qrzodt(k) )
            else
               pgfr(k)=0
            endif

!c
!c       if (qgz(k) = 0.0) skip the other step below about graupel
!c
         if (qgz(k) .eq. 0.0) goto 4000

!c
!c       comparing pgwet(wet process) and pdry(dry process),
!c       we will pick up the small one.
!c

!c       ---------------
!c      | dry processes |
!c       ---------------
!c
!c (3)   accretion of cloud water by graupel  (pgacw): lin (40)
!c       egw=1.0
!c       cdrag=0.6 drag coefficients for hairstone
!c       constg=sqrt(4.*grav*rhograul*0.33334*orho(k)/cdrag)
!c
         egw=1.0
         constg=sqrt(4.*grav*rhograul*0.33334*orho(k)*ocdrag)
         tmp1=pio4*xnog*gam3pt5*constg*olambdag(k)**3.5
         tmp2=qlz(k)*egw*tmp1
         pgacw(k)=amin1( tmp2,qlzodt(k) )
!c
!c (4)   accretion of ice crystal by graupel  (pgaci): lin (41)
!c       egi=1.   for wet growth
!c       egi=0.1  for dry growth
!c
         egi=0.1
         tmp2=qiz(k)*egi*tmp1
         pgaci(k)=amin1( tmp2,qizodt(k) )


!c
!c (5)   accretion of snow by graupel  (pgacs) : lin (29)
!c       compute processes (6) only when qsz > 0.0 and qgz > 0.0
!c
         egs=exp(0.09*temcc(k))
         tmpa=olambdas(k)*olambdas(k)
         tmpb=olambdag(k)*olambdag(k)
         tmpc=olambdas(k)*olambdag(k)
         tmp1=pi*pi*xnos*xnog*abs( vts(k)-vtg(k) )*orho(k)
         tmp2=tmpa*tmpa*olambdag(k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
         tmp3=tmp1*egs*rhosnow*tmp2
         pgacs(k)=amin1( tmp3,qszodt(k) )


!c
!c (6)   accretion of rain by graupel (pgacr): lin (42)
!c       compute processes (6) only when qrz > 0.0 and qgz > 0.0
!c       egr=1.
!c
         egr=1.
         tmpa=olambdar(k)*olambdar(k)
         tmpb=olambdag(k)*olambdag(k)
         tmpc=olambdar(k)*olambdag(k)
         tmp1=pi*pi*xnor*xnog*abs( vtr(k)-vtg(k) )*orho(k)
         tmp2=tmpa*tmpa*olambdag(k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
         tmp3=tmp1*egr*rhowater*tmp2
         pgacr(k)=amin1( tmp3,qrzodt(k) )

!c
!c (7)   calculate total dry process effect pdry(k)
!c
         pdry(k)=pgacw(k)+pgaci(k)+pgacs(k)+pgacr(k)

!c       ---------------
!c      | wet processes |
!c       ---------------
!c
!c (3)   accretion of ice crystal by graupel  (pgacip): lin (41)
!c       egi=1.   for wet growth
!c       egi=0.1  for dry growth
!c
         tmp2=10.*pgaci(k)
         pgacip(k)=amin1( tmp2,qizodt(k) )

!c
!c (4)   accretion of snow by graupel  ((pgacsp) : lin (29)
!c       compute processes (6) only when qsz > 0.0 and qgz > 0.0
!c       egs=exp(0.09*(tem(k)-273.15))  when t < 273.15 k
!c
         tmp3=pgacs(k)*1.0/egs
         pgacsp(k)=amin1( tmp3,qszodt(k) )

!c
!c (5)   wet growth of graupel (pgwet) : lin (43)
!c       may involve pgacs or pgaci and
!c       must include ppgacw or pgacr, or both.
!c       ( the amount of pgacw which is not able
!c       to freeze is shed to rain. )

             term0=constg*olambdag(k)**5.5/visc(k)

!c
!c    vf1s,vf2s=ventilation factors for graupel
!c    vf1s=0.78,vf2s=0.31 in lin
!c    cdrag=0.6  drag coefficient for hairstone
!c    constg2=vf1s*olambdag(k)*olambdag(k)+
!c            vf2s*schmidt(k)**0.33334*gam2pt75*sqrt(term0)

             delrs=rs0(k)-qvz(k)
             tmp0=1./(xlf+cw*temcc(k))
             tmp1=2.*pi*xnog*(rho(k)*xlv*diffwv(k)*delrs-xka(k)*  &
                  temcc(k))*orho(k)*tmp0
             constg2=vf1s*olambdag(k)*olambdag(k)+  &
                     vf2s*schmidt(k)**0.33334*gam2pt75*sqrt(term0)
             tmp3=tmp1*constg2+(pgacip(k)+pgacsp(k))*  &
                  (1-ci*temcc(k)*tmp0)
             tmp3=amax1(0.0,tmp3)
             pgwet(k)=amax1(tmp3,qlzodt(k)+qszodt(k)+qizodt(k) )

!c
!c     comparing pgwet(wet process) and pdry(dry process),
!c     we will apply the small one.
!c     if dry processes then delta4=1.0
!c     if wet processes then delta4=0.0
!
         if ( pdry(k) .lt. pgwet(k) ) then
            delta4=1.0
         else
            delta4=0.0
         endif

!c
!c
!c (6)   pgacrp(k)=pgwet(k)-pgacw(k)-pgacip(k)-pgacsp(k)
!c       if pgacrp(k) > 0. then some of the rain is frozen to hail
!c       if pgacrp(k) < 0. then some of the cloud water collected
!c                            by the hail is unable to freeze and is
!c                            shed as rain.
!c
            pgacrp(k)=pgwet(k)-pgacw(k)-pgacip(k)-pgacsp(k)

!c
!c (8)   deposition/sublimation of graupel  (pgdep/pgsub): lin (46)
!c       includes ventilation effect
!c       constg=sqrt(4.*grav*rhograul*0.33334*orho(k)/cdrag)
!c       constg2=vf1s*olambdag(k)*olambdag(k)+
!c             vf2s*schmidt(k)**0.33334*gam2pt75*constg
!c
!c       abg=2*pi*(si-1)/rho/(a"+b")
!c
            tmpa=rvapor*xka(k)*tem(k)*tem(k)
            tmpb=xls*xls*rho(k)*qsiz(k)*diffwv(k)
            tmpc=tmpa*qsiz(k)*diffwv(k)
            abg=2.0*pi*(qvoqsiz(k)-1.0)*tmpc/(tmpa+tmpb)
!c
!c     vf1s,vf2s=ventilation factors for graupel
!c     vf1s=0.78,vf2s=0.31 in lin
!c     cdrag=0.6  drag coefficient for hairstone
!c
            tmp2=abg*xnog*constg2
            pgdep(k)=amax1(0.0,tmp2)
            pgsub(k)=amin1(0.0,tmp2)
            pgsub(k)=amax1( pgsub(k),-qgzodt(k) )

 4000    continue
        else
!
!***********************************************************************
!*********      graupel production processes for t > 0 c      **********
!***********************************************************************
!
!c
!c (1) accretion of cloud water by graupel (pgacw): lin (40)
!c     egw=1.0
!c     cdrag=0.6 drag coefficients for hairstone
!c     constg=sqrt(4.*grav*rhograul*0.33334*orho(k)/cdrag)

            egw=1.0
            constg=sqrt(4.*grav*rhograul*0.33334*orho(k)*ocdrag)
            tmp1=pio4*xnog*gam3pt5*constg*olambdag(k)**3.5
            tmp2=qlz(k)*egw*tmp1
            pgacw(k)=amin1( tmp2,qlzodt(k) )

!c
!c (2) accretion of rain by graupel (pgacr): lin (42)
!c     compute processes (5) only when qrz > 0.0 and qgz > 0.0
!c     egr=1.
!c
            egr=1.
            tmpa=olambdar(k)*olambdar(k)
            tmpb=olambdag(k)*olambdag(k)
            tmpc=olambdar(k)*olambdag(k)
            tmp1=pi*pi*xnor*xnog*abs( vtr(k)-vtg(k) )*orho(k)
            tmp2=tmpa*tmpa*olambdag(k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
            tmp3=tmp1*egr*rhowater*tmp2
            pgacr(k)=amin1( tmp3,qrzodt(k) )


!c
!c (3) graupel melting to form rain (pgmlt): lin (47)
!c     pgmlt is negative value
!c     constg=sqrt(4.*grav*rhograul*0.33334*orho(k)/cdrag)
!c     constg2=vf1s*olambdag(k)*olambdag(k)+
!c             vf2s*schmidt(k)**0.33334*gam2pt75*constg
!c     cdrag=0.6  drag coefficients for hairstone
!
            delrs=rs0(k)-qvz(k)
            term1=2.0*pi*orho(k)*( xlv*diffwv(k)*rho(k)*delrs- &
                  xka(k)*temcc(k) )
            term0=sqrt(4.*grav*rhograul*0.33334*orho(k)*ocdrag) &
                  *olambdag(k)**5.5/visc(k)

            constg2=vf1s*olambdag(k)*olambdag(k)+ &
                    vf2s*schmidt(k)**0.33334*gam2pt75*sqrt(term0)
            tmp2=xnog*constg2
            tmp3=term1*oxlf*tmp2-cwoxlf*temcc(k)*( pgacw(k)+pgacr(k) )
            tmp4=amin1(0.0,tmp3)
            pgmlt(k)=amax1( tmp4,-qgzodt(k) )


!c
!c (4) evaporation of melting graupel (pgmltevp) : hr (a19)
!c     but use lin et al. coefficience
!c     pgmltevp is a negative value
!c     abg=2.0*pi*(qvoqsiz(k)-1.0)*tmpc/(tmpa+tmpb)
!c
            tmpa=rvapor*xka(k)*tem(k)*tem(k)
            tmpb=xlv*xlv*rho(k)*qswz(k)*diffwv(k)
            tmpc=tmpa*qswz(k)*diffwv(k)
            tmpd=amin1( 0.0,(qvoqswz(k)-0.90)*qswz(k)*odtb )

!c
!c     abg=2*pi*(si-1)/rho/(a"+b")
!c
            abg=2.0*pi*(qvoqswz(k)-0.90)*tmpc/(tmpa+tmpb)
!
!**** allow evaporation to occur when rh less than 90%
!**** here not using 100% because the evaporation cooling
!**** of temperature is not taking into account yet; hence,
!**** the qgw value is a little bit larger. this will avoid
!**** evaporation can generate cloud.
!
!c    vf1s,vf2s=ventilation factors for snow
!c    vf1s=0.78,vf2s=0.31 in lin
!c    constg=sqrt(4.*grav*rhograul*0.33334*orho(k)/cdrag)
!c    constg2=vf1s*olambdag(k)*olambdag(k)+
!c            vf2s*schmidt(k)**0.33334*gam2pt75*constg
!
            tmp2=abg*xnog*constg2
            tmp3=amin1(0.0,tmp2)
            tmp3=amax1( tmp3,tmpd )
            pgmltevp(k)=amax1( tmp3,-qgzodt(k) )

!c
!c (5) accretion of snow by graupel (pgacs) : lin (29)
!c     compute processes (3) only when qsz > 0.0 and qgz > 0.0
!c     egs=1.0
!c
           egs=1.
           tmpa=olambdas(k)*olambdas(k)
           tmpb=olambdag(k)*olambdag(k)
           tmpc=olambdas(k)*olambdag(k)
           tmp1=pi*pi*xnos*xnog*abs( vts(k)-vtg(k) )*orho(k)
           tmp2=tmpa*tmpa*olambdag(k)*(5.0*tmpa+2.0*tmpc+0.5*tmpb)
           tmp3=tmp1*egs*rhosnow*tmp2
           pgacs(k)=amin1( tmp3,qszodt(k) )

        endif


!
  900   continue

!cc
!c
!c**********************************************************************
!c*****     combine all processes together and avoid negative      *****
!c*****     water substances
!***********************************************************************
!c
      if ( temcc(k) .lt. 0.0) then
!,delta4,1.-delta4
!c
!c  gdelta4=gindex*delta4
!c  g1sdelt4=gindex*(1.-delta4)
!c
           gdelta4=gindex*delta4
           g1sdelt4=gindex*(1.-delta4)
!c
!c  combined water vapor depletions
!c
!cc graupel
           tmp=psdep(k)+pgdep(k)*gindex
           if ( tmp .gt. qvzodt(k) ) then
              factor=qvzodt(k)/tmp
              psdep(k)=psdep(k)*factor
              pgdep(k)=pgdep(k)*factor*gindex
           end if
!c
!c  combined cloud water depletions
!c
           tmp=praut(k)+psacw(k)+psfw(k)+pracw(k)+gindex*pgacw(k)
           if ( tmp .gt. qlzodt(k) ) then
              factor=qlzodt(k)/tmp
              praut(k)=praut(k)*factor
              psacw(k)=psacw(k)*factor
              psfw(k)=psfw(k)*factor
              pracw(k)=pracw(k)*factor
!cc graupel
              pgacw(k)=pgacw(k)*factor*gindex
           end if
!c
!c  combined cloud ice depletions
!c
           tmp=psaut(k)+psaci(k)+praci(k)+psfi(k)+pgaci(k)*gdelta4 &
               +pgacip(k)*g1sdelt4
           if (tmp .gt. qizodt(k) ) then
              factor=qizodt(k)/tmp
              psaut(k)=psaut(k)*factor
              psaci(k)=psaci(k)*factor
              praci(k)=praci(k)*factor
              psfi(k)=psfi(k)*factor
!cc graupel
              pgaci(k)=pgaci(k)*factor*gdelta4
              pgacip(k)=pgacip(k)*factor*g1sdelt4
           endif
!c
!c  combined all rain processes
!c
          tmp_r=piacr(k)+psacr(k)-prevp(k)-praut(k)-pracw(k) &
                +pgfr(k)*gindex+pgacr(k)*gdelta4 &
                +pgacrp(k)*g1sdelt4
          if (tmp_r .gt. qrzodt(k) ) then
             factor=qrzodt(k)/tmp_r
             piacr(k)=piacr(k)*factor
             psacr(k)=psacr(k)*factor
             prevp(k)=prevp(k)*factor
!cc graupel
             pgfr(k)=pgfr(k)*factor*gindex
             pgacr(k)=pgacr(k)*factor*gdelta4
             pgacrp(k)=pgacrp(k)*factor*g1sdelt4
          endif

!c
!c     if qrz < 1.0e-4 and qsz < 1.0e-4  then delta2=1.
!c                  (all pracs and psacr become to snow)
!c     if qrz >= 1.0e-4 or qsz >= 1.0e-4 then delta2=0.
!c                 (all pracs and psacr become to graupel)
!c
          if (qrz(k) .lt. 1.0e-4 .and. qsz(k) .lt. 1.0e-4) then
             delta2=1.0
          else
             delta2=0.0
          endif
!
!cc graupel

!c
!c  if qrz(k) < 1.0e-4 then delta3=1. means praci(k) -->  qs
!c                                          piacr(k) -->  qs
!c  if qrz(k) > 1.0e-4 then delta3=0. means praci(k) -->  qg
!c                                          piacr(k) -->  qg : lin (20)

          if (qrz(k) .lt. 1.0e-4) then
             delta3=1.0
          else
             delta3=0.0
          endif
!
!c
!c     if gindex = 0.(no graupel) then delta2=1.0
!c                                     delta3=1.0
!c
          if (gindex .eq. 0.) then
              delta2=1.0
              delta3=1.0
         endif
!
!c
!c   combined all snow processes
!c
          tmp_s=-pssub(k)-(psaut(k)+psaci(k)+psacw(k)+psfw(k)+ &
                 psfi(k)+praci(k)*delta3+piacr(k)*delta3+ &
                 psdep(k))+pgaut(k)*gindex+pgacs(k)*gdelta4+ &
                 pgacsp(k)*g1sdelt4+pracs(k)*(1.-delta2)- &
                 psacr(k)*delta2
          if ( tmp_s .gt. qszodt(k) ) then
             factor=qszodt(k)/tmp_s
             pssub(k)=pssub(k)*factor
             pracs(k)=pracs(k)*factor
!cc graupel
             pgaut(k)=pgaut(k)*factor*gindex
             pgacs(k)=pgacs(k)*factor*gdelta4
             pgacsp(k)=pgacsp(k)*factor*g1sdelt4
          endif

!cc graupel
!

!          if (gindex .eq. 0.) goto 998
!c
!c  combined all graupel processes
!c
           tmp_g=-pgaut(k)-pgfr(k)-pgacw(k)*delta4-pgaci(k)*delta4  &
                 -pgacr(k)*delta4-pgacs(k)*delta4  &
                 -pgwet(k)*(1.-delta4)-pgsub(k)-pgdep(k)  &
                 -psacr(k)*(1-delta2)-pracs(k)*(1-delta2)  &
                 -praci(k)*(1-delta3)-piacr(k)*(1-delta3)
           if (tmp_g .gt. qgzodt(k)) then
               factor=qgzodt(k)/tmp_g
               pgsub(k)=pgsub(k)*factor
           endif

  998      continue
!c
!c  calculate new water substances, thetae, tem, and qvsbar
!c

!cc graupel
         pvapor(k)=-pssub(k)-psdep(k)-prevp(k)-pgsub(k)*gindex &
                   -pgdep(k)*gindex
         qvz(k)=amax1( qvmin,qvz(k)+dtb*pvapor(k) )
         pclw(k)=-praut(k)-pracw(k)-psacw(k)-psfw(k)-pgacw(k)*gindex
         qlz(k)=amax1( 0.0,qlz(k)+dtb*pclw(k) )
         pcli(k)=-psaut(k)-psfi(k)-psaci(k)-praci(k)-pgaci(k)*gdelta4 &
                 -pgacip(k)*g1sdelt4
         qiz(k)=amax1( 0.0,qiz(k)+dtb*pcli(k) )
         tmp_r=piacr(k)+psacr(k)-prevp(k)-praut(k)-pracw(k) &
                +pgfr(k)*gindex+pgacr(k)*gdelta4 &
                +pgacrp(k)*g1sdelt4
 232     format(i2,1x,6(f9.3,1x))
         prain(k)=-tmp_r
         qrz(k)=amax1( 0.0,qrz(k)+dtb*prain(k) )
         tmp_s=-pssub(k)-(psaut(k)+psaci(k)+psacw(k)+psfw(k)+  &
                psfi(k)+praci(k)*delta3+piacr(k)*delta3+  &
                psdep(k))+pgaut(k)*gindex+pgacs(k)*gdelta4+  &
                pgacsp(k)*g1sdelt4+pracs(k)*(1.-delta2)-  &
                psacr(k)*delta2
         psnow(k)=-tmp_s
         qsz(k)=amax1( 0.0,qsz(k)+dtb*psnow(k) )
         qschg(k)=qschg(k)+psnow(k)
         qschg(k)=psnow(k)
!cc graupel
         tmp_g=-pgaut(k)-pgfr(k)-pgacw(k)*delta4-pgaci(k)*delta4 &
               -pgacr(k)*delta4-pgacs(k)*delta4 &
               -pgwet(k)*(1.-delta4)-pgsub(k)-pgdep(k) &
               -psacr(k)*(1-delta2)-pracs(k)*(1-delta2) &
               -praci(k)*(1-delta3)-piacr(k)*(1-delta3)
 252     format(i2,1x,6(f12.9,1x))
 262     format(i2,1x,7(f12.9,1x))
         pgraupel(k)=-tmp_g
         pgraupel(k)=pgraupel(k)*gindex
         qgz(k)=amax1( 0.0,qgz(k)+dtb*pgraupel(k))
!        qgchg(k)=qgchg(k)+pgraupel(k)
         qgchg(k)=pgraupel(k)
         qgz(k)=qgz(k)*gindex

         tmp=ocp/tothz(k)*xlf*(qschg(k)+qgchg(k))
         theiz(k)=theiz(k)+dtb*tmp
         thz(k)=theiz(k)-(xlvocp*qvz(k)-xlfocp*qiz(k))/tothz(k)
         tem(k)=thz(k)*tothz(k)

         temcc(k)=tem(k)-273.15

         if( temcc(k) .lt. -40.0 ) qswz(k)=qsiz(k)
         qlpqi=qlz(k)+qiz(k)
         if ( qlpqi .eq. 0.0 ) then
            qvsbar(k)=qsiz(k)
         else
            qvsbar(k)=( qiz(k)*qsiz(k)+qlz(k)*qswz(k) )/qlpqi
         endif

!
      else
!c
!c  combined cloud water depletions
!c
          tmp=praut(k)+psacw(k)+pracw(k)+pgacw(k)*gindex
          if ( tmp .gt. qlzodt(k) ) then
             factor=qlzodt(k)/tmp
             praut(k)=praut(k)*factor
             psacw(k)=psacw(k)*factor
             pracw(k)=pracw(k)*factor
!cc graupel
             pgacw(k)=pgacw(k)*factor*gindex
          end if
!c
!c  combined all snow processes
!c
          tmp_s=-(psmlt(k)+psmltevp(k))+pgacs(k)*gindex
          if (tmp_s .gt. qszodt(k) ) then
             factor=qszodt(k)/tmp_s
             psmlt(k)=psmlt(k)*factor
             psmltevp(k)=psmltevp(k)*factor
!cc graupel
             pgacs(k)=pgacs(k)*factor*gindex
          endif

!c
!c
!cc graupel
!c
!         if (gindex .eq. 0.) goto 997

!c
!c  combined all graupel processes
!c
            tmp_g=-pgmlt(k)-pgacs(k)-pgmltevp(k)
            if (tmp_g .gt. qgzodt(k)) then
              factor=qgzodt(k)/tmp_g
              pgmltevp(k)=pgmltevp(k)*factor
              pgmlt(k)=pgmlt(k)*factor
            endif
!c
  997     continue

!c
!c  combined all rain processes
!c
          tmp_r=-prevp(k)-(praut(k)+pracw(k)+psacw(k)-psmlt(k)) &
                +pgmlt(k)*gindex-pgacw(k)*gindex
          if (tmp_r .gt. qrzodt(k) ) then
             factor=qrzodt(k)/tmp_r
             prevp(k)=prevp(k)*factor
          endif
!c
!c
!c  calculate new water substances and thetae
!c


          pvapor(k)=-psmltevp(k)-prevp(k)-pgmltevp(k)
          qvz(k)=amax1( qvmin,qvz(k)+dtb*pvapor(k))
          pclw(k)=-praut(k)-pracw(k)-psacw(k)-pgacw(k)*gindex
          qlz(k)=amax1( 0.0,qlz(k)+dtb*pclw(k) )
          pcli(k)=0.0
          qiz(k)=amax1( 0.0,qiz(k)+dtb*pcli(k) )
          tmp_r=-prevp(k)-(praut(k)+pracw(k)+psacw(k)-psmlt(k)) &
                +pgmlt(k)*gindex-pgacw(k)*gindex
 242      format(i2,1x,7(f9.6,1x))
          prain(k)=-tmp_r
          tmpqrz=qrz(k)
          qrz(k)=amax1( 0.0,qrz(k)+dtb*prain(k) )
          tmp_s=-(psmlt(k)+psmltevp(k))+pgacs(k)*gindex
          psnow(k)=-tmp_s
          qsz(k)=amax1( 0.0,qsz(k)+dtb*psnow(k) )
!         qschg(k)=qschg(k)+psnow(k)
          qschg(k)=psnow(k)
!cc graupel

          tmp_g=-pgmlt(k)-pgacs(k)-pgmltevp(k)
!         write(*,272)k,pgmlt(k),pgacs(k),pgmltevp(k),
 272      format(i2,1x,3(f12.9,1x))
          pgraupel(k)=-tmp_g*gindex
          qgz(k)=amax1( 0.0,qgz(k)+dtb*pgraupel(k))
!         qgchg(k)=qgchg(k)+pgraupel(k)
          qgchg(k)=pgraupel(k)
          qgz(k)=qgz(k)*gindex
!
          tmp=ocp/tothz(k)*xlf*(qschg(k)+qgchg(k))
          theiz(k)=theiz(k)+dtb*tmp
          thz(k)=theiz(k)-(xlvocp*qvz(k)-xlfocp*qiz(k))/tothz(k)

          tem(k)=thz(k)*tothz(k)
          temcc(k)=tem(k)-273.15
!         qswz(k)=episp0k*oprez(k)*  &
!                exp( svp2*temcc(k)/(tem(k)-svp3) )
          es=1000.*svp1*exp( svp2*temcc(k)/(tem(k)-svp3) )
          qswz(k)=ep2*es/(prez(k)-es)
          qsiz(k)=qswz(k)
          qvsbar(k)=qswz(k)
!
      end if

!
!***********************************************************************
!**********              saturation adjustment                **********
!***********************************************************************
!
!    allow supersaturation exits linearly from 0% at 500 mb to 50%
!    above 300 mb
!    5.0e-5=1.0/(500mb-300mb)
!
         rsat=1.0+0.5*(50000.0-prez(k))*5.0e-5
         rsat=amax1(1.0,rsat)
         rsat=amin1(1.5,rsat)
         rsat=1.0
         if( qvz(k)+qlz(k)+qiz(k) .lt. rsat*qvsbar(k) ) then

!c
!c   unsaturated
!c
          qvz(k)=qvz(k)+qlz(k)+qiz(k)
          qlz(k)=0.0
          qiz(k)=0.0

          thz(k)=theiz(k)-(xlvocp*qvz(k)-xlfocp*qiz(k))/tothz(k)
          tem(k)=thz(k)*tothz(k)
          temcc(k)=tem(k)-273.15

          go to 1800
!
        else
!c
!c   saturated
!c
!
          pladj(k)=qlz(k)
          piadj(k)=qiz(k)
!

          call satadj(qvz, qlz, qiz, prez, theiz, thz, tothz, kts, kte, &
                      k, xlvocp, xlfocp, episp0k, ep2,svp1,svp2,svp3,svpt0 )

!
          pladj(k)=odtb*(qlz(k)-pladj(k))
          piadj(k)=odtb*(qiz(k)-piadj(k))
!
          pclw(k)=pclw(k)+pladj(k)
          pcli(k)=pcli(k)+piadj(k)
          pvapor(k)=pvapor(k)-( pladj(k)+piadj(k) )
!
          thz(k)=theiz(k)-(xlvocp*qvz(k)-xlfocp*qiz(k))/tothz(k)
          tem(k)=thz(k)*tothz(k)

          temcc(k)=tem(k)-273.15

!         qswz(k)=episp0k*oprez(k)*  &
!                 exp( svp2*temcc(k)/(tem(k)-svp3) )
          es=1000.*svp1*exp( svp2*temcc(k)/(tem(k)-svp3) )
          qswz(k)=ep2*es/(prez(k)-es)
          if (tem(k) .lt. 273.15 ) then
!            qsiz(k)=episp0k*oprez(k)* &
!                   exp( 21.8745584*(tem(k)-273.16)/(tem(k)-7.66) )
             es=1000.*svp1*exp( 21.8745584*(tem(k)-273.16)/(tem(k)-7.66) )
             qsiz(k)=ep2*es/(prez(k)-es)
             if (temcc(k) .lt. -40.0) qswz(k)=qsiz(k)
          else
             qsiz(k)=qswz(k)
          endif
          qlpqi=qlz(k)+qiz(k)
          if ( qlpqi .eq. 0.0 ) then
             qvsbar(k)=qsiz(k)
          else
             qvsbar(k)=( qiz(k)*qsiz(k)+qlz(k)*qswz(k) )/qlpqi
          endif

        end if

!
!***********************************************************************
!*****     melting and freezing of cloud ice and cloud water       *****
!***********************************************************************
        qlpqi=qlz(k)+qiz(k)
        if(qlpqi .le. 0.0) go to 1800
!
!c
!c (1)  homogeneous nucleation when t< -40 c (pihom)
!c
        if(temcc(k) .lt. -40.0) pihom(k)=qlz(k)*odtb
!c
!c (2)  melting of ice crystal when t> 0 c (pimlt)
!c
        if(temcc(k) .gt. 0.0) pimlt(k)=qiz(k)*odtb
!c
!c (3) production of cloud ice by bergeron process (pidw): hsie (p957)
!c     this process only considered when -31 c < t < 0 c
!c
        if(temcc(k) .lt. 0.0 .and. temcc(k) .gt. -31.0) then
!c!
!c!   parama1 and parama2 functions must be user supplied
!c!
          a1=parama1( temcc(k) )
          a2=parama2( temcc(k) )
!! change unit from cgs to mks
          a1=a1*0.001**(1.0-a2)
          xnin=xni0*exp(-bni*temcc(k))
          pidw(k)=xnin*orho(k)*(a1*xmnin**a2)
        end if
!
        pcli(k)=pcli(k)+pihom(k)-pimlt(k)+pidw(k)
        pclw(k)=pclw(k)-pihom(k)+pimlt(k)-pidw(k)
        qlz(k)=amax1( 0.0,qlz(k)+dtb*(-pihom(k)+pimlt(k)-pidw(k)) )
        qiz(k)=amax1( 0.0,qiz(k)+dtb*(pihom(k)-pimlt(k)+pidw(k)) )

!
        call satadj(qvz, qlz, qiz, prez, theiz, thz, tothz, kts, kte, &
                    k, xlvocp, xlfocp, episp0k ,ep2,svp1,svp2,svp3,svpt0)

        thz(k)=theiz(k)-(xlvocp*qvz(k)-xlfocp*qiz(k))/tothz(k)
        tem(k)=thz(k)*tothz(k)

        temcc(k)=tem(k)-273.15

!       qswz(k)=episp0k*oprez(k)* &
!              exp( svp2*temcc(k)/(tem(k)-svp3) )
        es=1000.*svp1*exp( svp2*temcc(k)/(tem(k)-svp3) )
        qswz(k)=ep2*es/(prez(k)-es)

        if (tem(k) .lt. 273.15 ) then
!          qsiz(k)=episp0k*oprez(k)* &
!                 exp( 21.8745584*(tem(k)-273.16)/(tem(k)-7.66) )
           es=1000.*svp1*exp( 21.8745584*(tem(k)-273.16)/(tem(k)-7.66) )
           qsiz(k)=ep2*es/(prez(k)-es)
           if (temcc(k) .lt. -40.0) qswz(k)=qsiz(k)
        else
           qsiz(k)=qswz(k)
        endif
        qlpqi=qlz(k)+qiz(k)
        if ( qlpqi .eq. 0.0 ) then
           qvsbar(k)=qsiz(k)
        else
           qvsbar(k)=( qiz(k)*qsiz(k)+qlz(k)*qswz(k) )/qlpqi
        endif

1800  continue
!
!***********************************************************************
!**********    integrate the productions of rain and snow     **********
!***********************************************************************
!c

2000  continue


!---------------------------------------------------------------------

!
!***********************************************************************
!******      write terms in cloud physics to time series dataset   *****
!***********************************************************************
!
!     open(unit=24,form='formatted',status='new',
!    &     file='cloud.dat')

!9030  format(10e12.6)

!      write(24,*)'tmp'
!      write(24,9030) (tem(k),k=kts+1,kte)
!      write(24,*)'qiz'
!      write(24,9030) (qiz(k),k=kts+1,kte)
!      write(24,*)'qsz'
!      write(24,9030) (qsz(k),k=kts+1,kte)
!      write(24,*)'qrz'
!      write(24,9030) (qrz(k),k=kts+1,kte)
!      write(24,*)'qgz'
!      write(24,9030) (qgz(k),k=kts+1,kte)
!      write(24,*)'qvoqsw'
!      write(24,9030) (qvoqswz(k),k=kts+1,kte)
!      write(24,*)'qvoqsi'
!      write(24,9030) (qvoqsiz(k),k=kts+1,kte)
!      write(24,*)'vtr'
!      write(24,9030) (vtr(k),k=kts+1,kte)
!      write(24,*)'vts'
!      write(24,9030) (vts(k),k=kts+1,kte)
!      write(24,*)'vtg'
!      write(24,9030) (vtg(k),k=kts+1,kte)
!      write(24,*)'pclw'
!      write(24,9030) (pclw(k),k=kts+1,kte)
!      write(24,*)'pvapor'
!      write(24,9030) (pvapor(k),k=kts+1,kte)
!      write(24,*)'pcli'
!      write(24,9030) (pcli(k),k=kts+1,kte)
!      write(24,*)'pimlt'
!      write(24,9030) (pimlt(k),k=kts+1,kte)
!      write(24,*)'pihom'
!      write(24,9030) (pihom(k),k=kts+1,kte)
!      write(24,*)'pidw'
!      write(24,9030) (pidw(k),k=kts+1,kte)
!      write(24,*)'prain'
!      write(24,9030) (prain(k),k=kts+1,kte)
!      write(24,*)'praut'
!      write(24,9030) (praut(k),k=kts+1,kte)
!      write(24,*)'pracw'
!      write(24,9030) (pracw(k),k=kts+1,kte)
!      write(24,*)'prevp'
!      write(24,9030) (prevp(k),k=kts+1,kte)
!      write(24,*)'psnow'
!      write(24,9030) (psnow(k),k=kts+1,kte)
!      write(24,*)'psaut'
!      write(24,9030) (psaut(k),k=kts+1,kte)
!      write(24,*)'psfw'
!      write(24,9030) (psfw(k),k=kts+1,kte)
!      write(24,*)'psfi'
!      write(24,9030) (psfi(k),k=kts+1,kte)
!      write(24,*)'praci'
!      write(24,9030) (praci(k),k=kts+1,kte)
!      write(24,*)'piacr'
!      write(24,9030) (piacr(k),k=kts+1,kte)
!      write(24,*)'psaci'
!      write(24,9030) (psaci(k),k=kts+1,kte)
!      write(24,*)'psacw'
!      write(24,9030) (psacw(k),k=kts+1,kte)
!      write(24,*)'psdep'
!      write(24,9030) (psdep(k),k=kts+1,kte)
!      write(24,*)'pssub'
!      write(24,9030) (pssub(k),k=kts+1,kte)
!      write(24,*)'pracs'
!      write(24,9030) (pracs(k),k=kts+1,kte)
!      write(24,*)'psacr'
!      write(24,9030) (psacr(k),k=kts+1,kte)
!      write(24,*)'psmlt'
!      write(24,9030) (psmlt(k),k=kts+1,kte)
!      write(24,*)'psmltevp'
!      write(24,9030) (psmltevp(k),k=kts+1,kte)
!      write(24,*)'pladj'
!      write(24,9030) (pladj(k),k=kts+1,kte)
!      write(24,*)'piadj'
!      write(24,9030) (piadj(k),k=kts+1,kte)
!      write(24,*)'pgraupel'
!      write(24,9030) (pgraupel(k),k=kts+1,kte)
!      write(24,*)'pgaut'
!      write(24,9030) (pgaut(k),k=kts+1,kte)
!      write(24,*)'pgfr'
!      write(24,9030) (pgfr(k),k=kts+1,kte)
!      write(24,*)'pgacw'
!      write(24,9030) (pgacw(k),k=kts+1,kte)
!      write(24,*)'pgaci'
!      write(24,9030) (pgaci(k),k=kts+1,kte)
!      write(24,*)'pgacr'
!      write(24,9030) (pgacr(k),k=kts+1,kte)
!      write(24,*)'pgacs'
!      write(24,9030) (pgacs(k),k=kts+1,kte)
!      write(24,*)'pgacip'
!      write(24,9030) (pgacip(k),k=kts+1,kte)
!      write(24,*)'pgacrp'
!      write(24,9030) (pgacrp(k),k=kts+1,kte)
!      write(24,*)'pgacsp'
!      write(24,9030) (pgacsp(k),k=kts+1,kte)
!      write(24,*)'pgwet'
!      write(24,9030) (pgwet(k),k=kts+1,kte)
!      write(24,*)'pdry'
!      write(24,9030) (pdry(k),k=kts+1,kte)
!      write(24,*)'pgsub'
!      write(24,9030) (pgsub(k),k=kts+1,kte)
!      write(24,*)'pgdep'
!      write(24,9030) (pgdep(k),k=kts+1,kte)
!      write(24,*)'pgmlt'
!      write(24,9030) (pgmlt(k),k=kts+1,kte)
!      write(24,*)'pgmltevp'
!      write(24,9030) (pgmltevp(k),k=kts+1,kte)



!**** below if qv < qvmin then qv=qvmin, ql=0.0, and qi=0.0
!
      do k=kts+1,kte
         if ( qvz(k) .lt. qvmin ) then
            qlz(k)=0.0
            qiz(k)=0.0
            qvz(k)=amax1( qvmin,qvz(k)+qlz(k)+qiz(k) )
         end if
      enddo
!
  end subroutine clphy1d


!---------------------------------------------------------------------
!                         saturated adjustment
!---------------------------------------------------------------------
      subroutine satadj(qvz, qlz, qiz, prez, theiz, thz, tothz,      &
                        kts, kte, k, xlvocp, xlfocp, episp0k, ep2,svp1,svp2,svp3,svpt0)
!---------------------------------------------------------------------
   implicit none
!---------------------------------------------------------------------
!  this program use newton's method for finding saturated temperature
!  and saturation mixing ratio.
!
! in this saturation adjustment scheme we assume
! (1)  the saturation mixing ratio is the mass weighted average of
!      saturation values over liquid water (qsw), and ice (qsi)
!      following lord et al., 1984 and tao, 1989
!
! (2) the percentage of cloud liquid and cloud ice will
!      be fixed during the saturation calculation
!---------------------------------------------------------------------
!

  integer, intent(in   )             :: kts, kte, k

  real,      dimension( kts:kte ),                                   &
                       intent(inout) :: qvz, qlz, qiz
!
  real,      dimension( kts:kte ),                                   &
                       intent(in   ) :: prez, theiz, tothz

  real,     intent(in   )            :: xlvocp, xlfocp, episp0k
  real,     intent(in   )            :: ep2,svp1,svp2,svp3,svpt0

! local vars

  integer                            :: n

  real, dimension( kts:kte )         :: thz, tem, temcc, qsiz,       &
                                        qswz, qvsbar

  real ::   qsat, qlpqi, ratql, t0, t1, tmp1, ratqi, tsat, absft,    &
            denom1, denom2, dqvsbar, ftsat, dftsat, qpz, qvold,      &
            gindex, es
!
!---------------------------------------------------------------------

      thz(k)=theiz(k)-(xlvocp*qvz(k)-xlfocp*qiz(k))/tothz(k)

      tem(k)=tothz(k)*thz(k)
      if (tem(k) .gt. 273.15) then
!        qsat=episp0k/prez(k)*  &
!            exp( svp2*(tem(k)-273.15)/(tem(k)-svp3) )
         es=1000.*svp1*exp( svp2*(tem(k)-svpt0)/(tem(k)-svp3) )
         qsat=ep2*es/(prez(k)-es)
      else
        qsat=episp0k/prez(k)*  &
             exp( 21.8745584*(tem(k)-273.15)/(tem(k)-7.66) )
      end if
      qpz=qvz(k)+qlz(k)+qiz(k)
      qvold=qvz(k)
      if (qpz .lt. qsat) then
         qvz(k)=qpz
         qiz(k)=0.0
         qlz(k)=0.0
         go to 400
      end if
      qlpqi=qlz(k)+qiz(k)
      if( qlpqi .ge. 1.0e-5) then
        ratql=qlz(k)/qlpqi
        ratqi=qiz(k)/qlpqi
      else
        t0=273.15
!       t1=233.15
        t1=248.15
        tmp1=( t0-tem(k) )/(t0-t1)
        tmp1=amin1(1.0,tmp1)
        tmp1=amax1(0.0,tmp1)
        ratqi=tmp1
        ratql=1.0-tmp1
      end if
!
!
!--  saturation mixing ratios over water and ice
!--  at the outset we will follow bolton 1980 mwr for
!--  the water and murray jas 1967 for the ice
!
!-- dqvsbar=d(qvsbar)/dt
!-- ftsat=f(tsat)
!-- dftsat=d(f(t))/dt
!
!  first guess of tsat

      tsat=tem(k)
      absft=1.0
!
      do 200 n=1,20
         denom1=1.0/(tsat-svp3)
         denom2=1.0/(tsat-7.66)
!        qswz(k)=episp0k/prez(k)*  &
!                exp( svp2*denom1*(tsat-273.15) )
         es=1000.*svp1*exp( svp2*denom1*(tsat-svpt0) )
         qswz(k)=ep2*es/(prez(k)-es)
         if (tem(k) .lt. 273.15) then
!           qsiz(k)=episp0k/prez(k)*  &
!                   exp( 21.8745584*denom2*(tsat-273.15) )
            es=1000.*svp1*exp( 21.8745584*denom2*(tsat-273.15) )
            qsiz(k)=ep2*es/(prez(k)-es)
            if (tem(k) .lt. 233.15) qswz(k)=qsiz(k)
         else
            qsiz(k)=qswz(k)
         endif
         qvsbar(k)=ratql*qswz(k)+ratqi*qsiz(k)
!
!        if( absft .lt. 0.01 .and. n .gt. 3 ) go to 300
         if( absft .lt. 0.01 ) go to 300
!
         dqvsbar=ratql*qswz(k)*svp2*243.5*denom1*denom1+  &
                 ratqi*qsiz(k)*21.8745584*265.5*denom2*denom2
         ftsat=tsat+(xlvocp+ratqi*xlfocp)*qvsbar(k)-  &
               tothz(k)*theiz(k)-xlfocp*ratqi*(qvz(k)+qlz(k)+qiz(k))
         dftsat=1.0+(xlvocp+ratqi*xlfocp)*dqvsbar
         tsat=tsat-ftsat/dftsat
         absft=abs(ftsat)

200   continue
9020  format(1x,'point can not converge, absft,n=',e12.5,i5)
!
300   continue
      if( qpz .gt. qvsbar(k) ) then
        qvz(k)=qvsbar(k)
        qiz(k)=ratqi*( qpz-qvz(k) )
        qlz(k)=ratql*( qpz-qvz(k) )
      else
        qvz(k)=qpz
        qiz(k)=0.0
        qlz(k)=0.0
      end if
 400  continue

   end subroutine satadj


!----------------------------------------------------------------
   real function parama1(temp)
!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
!  this program calculate the parameter for crystal growth rate
!  in bergeron process
!----------------------------------------------------------------

      real, intent (in   )   :: temp
      real, dimension(32)    :: a1
      integer                :: i1, i1p1
      real                   :: ratio

      data a1/0.100e-10,0.7939e-7,0.7841e-6,0.3369e-5,0.4336e-5, &
              0.5285e-5,0.3728e-5,0.1852e-5,0.2991e-6,0.4248e-6, &
              0.7434e-6,0.1812e-5,0.4394e-5,0.9145e-5,0.1725e-4, &
              0.3348e-4,0.1725e-4,0.9175e-5,0.4412e-5,0.2252e-5, &
              0.9115e-6,0.4876e-6,0.3473e-6,0.4758e-6,0.6306e-6, &
              0.8573e-6,0.7868e-6,0.7192e-6,0.6513e-6,0.5956e-6, &
              0.5333e-6,0.4834e-6/

      i1=int(-temp)+1
      i1p1=i1+1
      ratio=-(temp)-float(i1-1)
      parama1=a1(i1)+ratio*( a1(i1p1)-a1(i1) )

   end function parama1

!----------------------------------------------------------------
   real function parama2(temp)
!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
!  this program calculate the parameter for crystal growth rate
!  in bergeron process
!----------------------------------------------------------------

      real, intent (in   )   :: temp
      real, dimension(32)    :: a2
      integer                :: i1, i1p1
      real                   :: ratio

      data a2/0.0100,0.4006,0.4831,0.5320,0.5307,0.5319,0.5249, &
              0.4888,0.3849,0.4047,0.4318,0.4771,0.5183,0.5463, &
              0.5651,0.5813,0.5655,0.5478,0.5203,0.4906,0.4447, &
              0.4126,0.3960,0.4149,0.4320,0.4506,0.4483,0.4460, &
              0.4433,0.4413,0.4382,0.4361/
      i1=int(-temp)+1
      i1p1=i1+1
      ratio=-(temp)-float(i1-1)
      parama2=a2(i1)+ratio*( a2(i1p1)-a2(i1) )

   end function parama2

!----------------------------------------------------------------
   real function ggamma(x)
!----------------------------------------------------------------
   implicit none
!----------------------------------------------------------------
      real, intent(in   ) :: x
      real, dimension(8)  :: b
      integer             ::j, k1
      real                ::pf, g1to2 ,temp

      data b/-.577191652,.988205891,-.897056937,.918206857,  &
             -.756704078,.482199394,-.193527818,.035868343/

      pf=1.
      temp=x
      do 10 j=1,200
      if (temp .le. 2) go to 20
      temp=temp-1.
   10 pf=pf*temp
  100 format(//,5x,'module_mp_lin: input to gamma function too large, x=',e12.5)
      write(wrf_err_message,100)x
      !call endrun(wrf_err_message)
      call endrun('ggamma in mp-lin')
   20 g1to2=1.
      temp=temp - 1.
      do 30 k1=1,8
   30 g1to2=g1to2 + b(k1)*temp**k1
      ggamma=pf*g1to2

      end function ggamma

!----------------------------------------------------------------

end module module_mp_lin_v2

