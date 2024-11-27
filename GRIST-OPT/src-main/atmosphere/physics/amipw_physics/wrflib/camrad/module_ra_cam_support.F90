! adopted from wrf-v3.4.1, add AMIPW_PHYSICS for GRIST
module module_ra_cam_support
#ifdef AMIPW_PHYSICS
  use grist_handle_error, only: endrun
  use grist_mpi
#else
  use module_cam_support, only: endrun
#endif

  implicit none
      integer, parameter :: r8 = 8
      real(r8), parameter:: inf = 1.e20 ! cam sets this differently in infnan.f90
      integer, parameter:: bigint = o'17777777777'           ! largest possible 32-bit integer 

      integer :: ixcldliq 
      integer :: ixcldice
!     integer :: levsiz    ! size of level dimension on dataset
      integer, parameter :: nbands = 2          ! number of spectral bands
      integer, parameter :: naer_all = 12 + 1
      integer, parameter :: naer = 10 + 1
      integer, parameter :: bnd_nbr_lw=7 
      integer, parameter :: ndstsz = 4    ! number of dust size bins
      integer :: idxsul
      integer :: idxsslt
      integer :: idxdustfirst
      integer :: idxcarbonfirst
      integer :: idxocpho
      integer :: idxbcpho
      integer :: idxocphi
      integer :: idxbcphi
      integer :: idxbg  
      integer :: idxvolc

  integer :: mxaerl                            ! maximum level of background aerosol

! indices to sections of array that represent
! groups of aerosols

  integer, parameter :: &
      numdust         = 4, &
      numcarbon      = 4

! portion of each species group to use in computation
! of relative radiative forcing.

  real(r8) :: sulscl_rf  = 0._r8 !
  real(r8) :: carscl_rf  = 0._r8
  real(r8) :: ssltscl_rf = 0._r8
  real(r8) :: dustscl_rf = 0._r8
  real(r8) :: bgscl_rf   = 0._r8
  real(r8) :: volcscl_rf = 0._r8

! "background" aerosol species mmr.
  real(r8) :: tauback = 0._r8

! portion of each species group to use in computation
! of aerosol forcing in driving the climate
  real(r8) :: sulscl  = 1._r8
  real(r8) :: carscl  = 1._r8
  real(r8) :: ssltscl = 1._r8
  real(r8) :: dustscl = 1._r8
  real(r8) :: volcscl = 1._r8

!from volcrad.f90 module
     integer, parameter :: idx_lw_0500_0650=3
     integer, parameter :: idx_lw_0650_0800=4
     integer, parameter :: idx_lw_0800_1000=5
     integer, parameter :: idx_lw_1000_1200=6
     integer, parameter :: idx_lw_1200_2000=7

! first two values represent the overlap of volcanics with the non-window
! (0-800, 1200-2200 cm^-1) and window (800-1200 cm^-1) regions.|  coefficients
! were derived using crm_volc_minimize.pro with spectral flux optimization
! on first iteration, total heating rate on subsequent iterations (2-9).
! five profiles for hls, hlw, mls, mlw, and tro conditions were given equal
! weight.  rms heating rate errors for a visible stratospheric optical
! depth of 1.0 are 0.02948 k/day.
!
      real(r8) :: abs_cff_mss_aer(bnd_nbr_lw) = &
         (/ 70.257384, 285.282943, &
         1.0273851e+02, 6.3073303e+01, 1.2039569e+02, &
         3.6343643e+02, 2.7138528e+02 /)

!from radae.f90 module
      real(r8), parameter:: min_tp_h2o = 160.0        ! min t_p for pre-calculated abs/emis
      real(r8), parameter:: max_tp_h2o = 349.999999   ! max t_p for pre-calculated abs/emis
      real(r8), parameter:: dtp_h2o = 21.111111111111 ! difference in adjacent elements of tp_h2o
      real(r8), parameter:: min_te_h2o = -120.0       ! min t_e-t_p for pre-calculated abs/emis
      real(r8), parameter:: max_te_h2o = 79.999999    ! max t_e-t_p for pre-calculated abs/emis
      real(r8), parameter:: dte_h2o  = 10.0           ! difference in adjacent elements of te_h2o
      real(r8), parameter:: min_rh_h2o = 0.0          ! min rh for pre-calculated abs/emis
      real(r8), parameter:: max_rh_h2o = 1.19999999   ! max rh for pre-calculated abs/emis
      real(r8), parameter:: drh_h2o = 0.2             ! difference in adjacent elements of rh
      real(r8), parameter:: min_lu_h2o = -8.0         ! min log_10(u) for pre-calculated abs/emis
      real(r8), parameter:: min_u_h2o  = 1.0e-8       ! min pressure-weighted path-length
      real(r8), parameter:: max_lu_h2o =  3.9999999   ! max log_10(u) for pre-calculated abs/emis
      real(r8), parameter:: dlu_h2o  = 0.5            ! difference in adjacent elements of lu_h2o
      real(r8), parameter:: min_lp_h2o = -3.0         ! min log_10(p) for pre-calculated abs/emis
      real(r8), parameter:: min_p_h2o = 1.0e-3        ! min log_10(p) for pre-calculated abs/emis
      real(r8), parameter:: max_lp_h2o = -0.0000001   ! max log_10(p) for pre-calculated abs/emis
      real(r8), parameter:: dlp_h2o = 0.3333333333333 ! difference in adjacent elements of lp_h2o
      integer, parameter :: n_u = 25   ! number of u in abs/emis tables
      integer, parameter :: n_p = 10   ! number of p in abs/emis tables
      integer, parameter :: n_tp = 10  ! number of t_p in abs/emis tables
      integer, parameter :: n_te = 21  ! number of t_e in abs/emis tables
      integer, parameter :: n_rh = 7   ! number of rh in abs/emis tables
      real(r8):: c16,c17,c26,c27,c28,c29,c30,c31
      real(r8):: fwcoef      ! farwing correction constant
      real(r8):: fwc1,fwc2   ! farwing correction constants
      real(r8):: fc1         ! farwing correction constant
      real(r8):: amco2 ! molecular weight of co2   (g/mol)
      real(r8):: amd   ! molecular weight of dry air (g/mol)
      real(r8):: p0    ! standard pressure (dynes/cm**2)

! these are now allocatable. jm 20090612
  real(r8), allocatable, dimension(:,:,:,:,:)  :: ah2onw   ! (n_p, n_tp, n_u, n_te, n_rh)   ! absorptivity (non-window)
  real(r8), allocatable, dimension(:,:,:,:,:)  :: eh2onw   ! (n_p, n_tp, n_u, n_te, n_rh)   ! emissivity   (non-window)
  real(r8), allocatable, dimension(:,:,:,:,:)  :: ah2ow    ! (n_p, n_tp, n_u, n_te, n_rh)    ! absorptivity (window, for adjacent layers)
  real(r8), allocatable, dimension(:,:,:,:,:)  :: cn_ah2ow ! (n_p, n_tp, n_u, n_te, n_rh)    ! continuum transmission for absorptivity (window)
  real(r8), allocatable, dimension(:,:,:,:,:)  :: cn_eh2ow ! (n_p, n_tp, n_u, n_te, n_rh)    ! continuum transmission for emissivity   (window)
  real(r8), allocatable, dimension(:,:,:,:,:)  :: ln_ah2ow ! (n_p, n_tp, n_u, n_te, n_rh)    ! line-only transmission for absorptivity (window)
  real(r8), allocatable, dimension(:,:,:,:,:)  :: ln_eh2ow ! (n_p, n_tp, n_u, n_te, n_rh)    ! line-only transmission for emissivity   (window)

!
! constant coefficients for water vapor overlap with trace gases.
! reference: ramanathan, v. and  p.downey, 1986: a nonisothermal
!            emissivity and absorptivity formulation for water vapor
!            journal of geophysical research, vol. 91., d8, pp 8649-8666
!
  real(r8):: coefh(2,4) = reshape(  &
         (/ (/5.46557e+01,-7.30387e-02/), &
            (/1.09311e+02,-1.46077e-01/), &
            (/5.11479e+01,-6.82615e-02/), &
            (/1.02296e+02,-1.36523e-01/) /), (/2,4/) )
!
  real(r8):: coefj(3,2) = reshape( &
            (/ (/2.82096e-02,2.47836e-04,1.16904e-06/), &
               (/9.27379e-02,8.04454e-04,6.88844e-06/) /), (/3,2/) )
!
  real(r8):: coefk(3,2) = reshape( &
            (/ (/2.48852e-01,2.09667e-03,2.60377e-06/) , &
               (/1.03594e+00,6.58620e-03,4.04456e-06/) /), (/3,2/) )

  integer, parameter :: ntemp = 192 ! number of temperatures in h2o sat. table for tp
  real(r8) :: estblh2o(0:ntemp)       ! saturation vapor pressure for h2o for tp rang
  integer, parameter :: o_fa = 6   ! degree+1 of poly of t_e for absorptivity as u->inf.
  integer, parameter :: o_fe = 6   ! degree+1 of poly of t_e for emissivity as u->inf.

!-----------------------------------------------------------------------------
! data for f in c/h/e fit -- value of a and e as u->infinity
! new c/lt/e fit (hitran 2k, ckd 2.4) -- no change
!     these values are determined by integrals of planck functions or
!     derivatives of planck functions only.
!-----------------------------------------------------------------------------
!
! fa/fe coefficients for 2 bands (0-800 & 1200-2200, 800-1200 cm^-1)
!
! coefficients of polynomial for f_a in t_e
!
  real(r8), parameter:: fat(o_fa,nbands) = reshape( (/ &
       (/-1.06665373e-01,  2.90617375e-02, -2.70642049e-04,   &   ! 0-800&1200-2200 cm^-1
          1.07595511e-06, -1.97419681e-09,  1.37763374e-12/), &   !   0-800&1200-2200 cm^-1
       (/ 1.10666537e+00, -2.90617375e-02,  2.70642049e-04,   &   ! 800-1200 cm^-1
         -1.07595511e-06,  1.97419681e-09, -1.37763374e-12/) /) & !   800-1200 cm^-1
       , (/o_fa,nbands/) )
!
! coefficients of polynomial for f_e in t_e
!
  real(r8), parameter:: fet(o_fe,nbands) = reshape( (/ &
      (/3.46148163e-01,  1.51240299e-02, -1.21846479e-04,   &   ! 0-800&1200-2200 cm^-1
        4.04970123e-07, -6.15368936e-10,  3.52415071e-13/), &   !   0-800&1200-2200 cm^-1
      (/6.53851837e-01, -1.51240299e-02,  1.21846479e-04,   &   ! 800-1200 cm^-1
       -4.04970123e-07,  6.15368936e-10, -3.52415071e-13/) /) & !   800-1200 cm^-1
      , (/o_fa,nbands/) )


      real(r8) ::  gravit     ! acceleration of gravity (cgs)
      real(r8) ::  rga        ! 1./gravit
      real(r8) ::  gravmks    ! acceleration of gravity (mks)
      real(r8) ::  cpair      ! specific heat of dry air
      real(r8) ::  epsilo     ! ratio of mol. wght of h2o to dry air
      real(r8) ::  epsqs      ! ratio of mol. wght of h2o to dry air
      real(r8) ::  sslp       ! standard sea-level pressure
      real(r8) ::  stebol     ! stefan-boltzmann's constant
      real(r8) ::  rgsslp     ! 0.5/(gravit*sslp)
      real(r8) ::  dpfo3      ! voigt correction factor for o3
      real(r8) ::  dpfco2     ! voigt correction factor for co2
      real(r8) ::  dayspy     ! number of days per 1 year
      real(r8) ::  pie        ! 3.14.....
      real(r8) ::  mwdry      ! molecular weight dry air ~ kg/kmole (shr_const_mwdair)
      real(r8) ::  scon       ! solar constant (not used in wrf)
      real(r8) ::  co2mmr
real(r8) ::   mwco2              ! molecular weight of carbon dioxide
real(r8) ::   mwh2o              ! molecular weight water vapor (shr_const_mwwv)
real(r8) ::   mwch4              ! molecular weight ch4
real(r8) ::   mwn2o              ! molecular weight n2o
real(r8) ::   mwf11              ! molecular weight cfc11
real(r8) ::   mwf12              ! molecular weight cfc12
real(r8) ::   cappa              ! r/cp
real(r8) ::   rair               ! gas constant for dry air (j/k/kg)
real(r8) ::   tmelt              ! freezing t of fresh water ~ k
real(r8) ::   r_universal        ! universal gas constant ~ j/k/kmole
real(r8) ::   latvap             ! latent heat of evaporation ~ j/kg
real(r8) ::   latice             ! latent heat of fusion ~ j/kg
real(r8) ::   zvir               ! r_v/r_d - 1.
  integer plenest  ! length of saturation vapor pressure table
  parameter (plenest=250)
! 
! table of saturation vapor pressure values es from tmin degrees
! to tmax+1 degrees k in one degree increments.  ttrice defines the
! transition region where es is a combination of ice & water values
!
real(r8) estbl(plenest)      ! table values of saturation vapor pressure
real(r8) tmin       ! min temperature (k) for table
real(r8) tmax       ! max temperature (k) for table
real(r8) pcf(6)     ! polynomial coeffs -> es transition water to ice
!real(r8), allocatable :: pin(:)           ! ozone pressure level (levsiz)
!real(r8), allocatable :: ozmix(:,:,:)     ! mixing ratio
!real(r8), allocatable, target :: abstot_3d(:,:,:,:) ! non-adjacent layer absorptivites
!real(r8), allocatable, target :: absnxt_3d(:,:,:,:) ! nearest layer absorptivities
!real(r8), allocatable, target :: emstot_3d(:,:,:)   ! total emissivity

!from aer_optics.f90 module
integer, parameter :: idxvis = 8     ! index to visible band
integer, parameter :: nrh = 1000   ! number of relative humidity values for look-up-table
integer, parameter :: nspint = 19   ! number of spectral intervals

! these are now allocatable,  jm 20090612
real(r8), allocatable, dimension(:,:) :: ksul    ! (nrh, nspint)    ! sulfate specific extinction  ( m^2 g-1 )
real(r8), allocatable, dimension(:,:) :: wsul    ! (nrh, nspint)    ! sulfate single scattering albedo
real(r8), allocatable, dimension(:,:) :: gsul    ! (nrh, nspint)    ! sulfate asymmetry parameter
real(r8), allocatable, dimension(:,:) :: ksslt   ! (nrh, nspint)   ! sea-salt specific extinction  ( m^2 g-1 )
real(r8), allocatable, dimension(:,:) :: wsslt   ! (nrh, nspint)   ! sea-salt single scattering albedo
real(r8), allocatable, dimension(:,:) :: gsslt   ! (nrh, nspint)   ! sea-salt asymmetry parameter
real(r8), allocatable, dimension(:,:) :: kcphil  ! (nrh, nspint)  ! hydrophilic carbon specific extinction  ( m^2 g-1 )
real(r8), allocatable, dimension(:,:) :: wcphil  ! (nrh, nspint)  ! hydrophilic carbon single scattering albedo
real(r8), allocatable, dimension(:,:) :: gcphil  ! (nrh, nspint)  ! hydrophilic carbon asymmetry parameter

real(r8) :: kbg(nspint)          ! background specific extinction  ( m^2 g-1 )
real(r8) :: wbg(nspint)          ! background single scattering albedo
real(r8) :: gbg(nspint)          ! background asymmetry parameter
real(r8) :: kcphob(nspint)       ! hydrophobic carbon specific extinction  ( m^2 g-1 )
real(r8) :: wcphob(nspint)       ! hydrophobic carbon single scattering albedo
real(r8) :: gcphob(nspint)       ! hydrophobic carbon asymmetry parameter
real(r8) :: kcb(nspint)          ! black carbon specific extinction  ( m^2 g-1 )
real(r8) :: wcb(nspint)          ! black carbon single scattering albedo
real(r8) :: gcb(nspint)          ! black carbon asymmetry parameter
real(r8) :: kvolc(nspint)        ! volcanic specific extinction  ( m^2 g-1)
real(r8) :: wvolc(nspint)        ! volcanic single scattering albedo
real(r8) :: gvolc(nspint)        ! volcanic asymmetry parameter

real(r8) :: kdst(ndstsz, nspint) ! dust specific extinction  ( m^2 g-1 )
real(r8) :: wdst(ndstsz, nspint) ! dust single scattering albedo
real(r8) :: gdst(ndstsz, nspint) ! dust asymmetry parameter
!
!from comozp.f90 module
      real(r8) cplos    ! constant for ozone path length integral
      real(r8) cplol    ! constant for ozone path length integral

!from ghg_surfvals.f90 module
   real(r8) :: co2vmr = 3.550e-4         ! co2   volume mixing ratio
   real(r8) :: n2ovmr = 0.311e-6         ! n2o   volume mixing ratio
   real(r8) :: ch4vmr = 1.714e-6         ! ch4   volume mixing ratio
   real(r8) :: f11vmr = 0.280e-9         ! cfc11 volume mixing ratio
   real(r8) :: f12vmr = 0.503e-9         ! cfc12 volume mixing ratio

integer, parameter :: cyr = 233  ! number of years of co2 data

   integer  :: yrdata(cyr) = &
 (/ 1869, 1870, 1871, 1872, 1873, 1874, 1875, &
    1876, 1877, 1878, 1879, 1880, 1881, 1882, &
    1883, 1884, 1885, 1886, 1887, 1888, 1889, &
    1890, 1891, 1892, 1893, 1894, 1895, 1896, &
    1897, 1898, 1899, 1900, 1901, 1902, 1903, &
    1904, 1905, 1906, 1907, 1908, 1909, 1910, &
    1911, 1912, 1913, 1914, 1915, 1916, 1917, &
    1918, 1919, 1920, 1921, 1922, 1923, 1924, &
    1925, 1926, 1927, 1928, 1929, 1930, 1931, &
    1932, 1933, 1934, 1935, 1936, 1937, 1938, &
    1939, 1940, 1941, 1942, 1943, 1944, 1945, &
    1946, 1947, 1948, 1949, 1950, 1951, 1952, &
    1953, 1954, 1955, 1956, 1957, 1958, 1959, &
    1960, 1961, 1962, 1963, 1964, 1965, 1966, &
    1967, 1968, 1969, 1970, 1971, 1972, 1973, &
    1974, 1975, 1976, 1977, 1978, 1979, 1980, &
    1981, 1982, 1983, 1984, 1985, 1986, 1987, &
    1988, 1989, 1990, 1991, 1992, 1993, 1994, &
    1995, 1996, 1997, 1998, 1999, 2000, 2001, &
    2002, 2003, 2004, 2005, 2006, 2007, 2008, &
    2009, 2010, 2011, 2012, 2013, 2014, 2015, &
    2016, 2017, 2018, 2019, 2020, 2021, 2022, &
    2023, 2024, 2025, 2026, 2027, 2028, 2029, &
    2030, 2031, 2032, 2033, 2034, 2035, 2036, &
    2037, 2038, 2039, 2040, 2041, 2042, 2043, &
    2044, 2045, 2046, 2047, 2048, 2049, 2050, &
    2051, 2052, 2053, 2054, 2055, 2056, 2057, &
    2058, 2059, 2060, 2061, 2062, 2063, 2064, &
    2065, 2066, 2067, 2068, 2069, 2070, 2071, &
    2072, 2073, 2074, 2075, 2076, 2077, 2078, &
    2079, 2080, 2081, 2082, 2083, 2084, 2085, &
    2086, 2087, 2088, 2089, 2090, 2091, 2092, &
    2093, 2094, 2095, 2096, 2097, 2098, 2099, &
    2100, 2101                               /)

! a2 future scenario
    real(r8)  :: co2(cyr) = &
 (/ 289.263, 289.263, 289.416, 289.577, 289.745, 289.919, 290.102, &
    290.293, 290.491, 290.696, 290.909, 291.129, 291.355, 291.587, 291.824, &
    292.066, 292.313, 292.563, 292.815, 293.071, 293.328, 293.586, 293.843, &
    294.098, 294.35, 294.598, 294.842, 295.082, 295.32, 295.558, 295.797,   &
    296.038, 296.284, 296.535, 296.794, 297.062, 297.338, 297.62, 297.91,   &
    298.204, 298.504, 298.806, 299.111, 299.419, 299.729, 300.04, 300.352,  &
    300.666, 300.98, 301.294, 301.608, 301.923, 302.237, 302.551, 302.863,  &
    303.172, 303.478, 303.779, 304.075, 304.366, 304.651, 304.93, 305.206,  &
    305.478, 305.746, 306.013, 306.28, 306.546, 306.815, 307.087, 307.365,  &
    307.65, 307.943, 308.246, 308.56, 308.887, 309.228, 309.584, 309.956,   &
    310.344, 310.749, 311.172, 311.614, 312.077, 312.561, 313.068, 313.599, &
    314.154, 314.737, 315.347, 315.984, 316.646, 317.328, 318.026, 318.742, &
    319.489, 320.282, 321.133, 322.045, 323.021, 324.06, 325.155, 326.299,  &
    327.484, 328.698, 329.933, 331.194, 332.499, 333.854, 335.254, 336.69,  &
    338.15, 339.628, 341.125, 342.65, 344.206, 345.797, 347.397, 348.98,    &
    350.551, 352.1, 354.3637, 355.7772, 357.1601, 358.5306, 359.9046,       &
    361.4157, 363.0445, 364.7761, 366.6064, 368.5322, 370.534, 372.5798,    &
    374.6564, 376.7656, 378.9087, 381.0864, 383.2994, 385.548, 387.8326,    &
    390.1536, 392.523, 394.9625, 397.4806, 400.075, 402.7444, 405.4875,     &
    408.3035, 411.1918, 414.1518, 417.1831, 420.2806, 423.4355, 426.6442,   &
    429.9076, 433.2261, 436.6002, 440.0303, 443.5168, 447.06, 450.6603,     &
    454.3059, 457.9756, 461.6612, 465.3649, 469.0886, 472.8335, 476.6008,   &
    480.3916, 484.2069, 488.0473, 491.9184, 495.8295, 499.7849, 503.7843,   &
    507.8278, 511.9155, 516.0476, 520.2243, 524.4459, 528.7127, 533.0213,   &
    537.3655, 541.7429, 546.1544, 550.6005, 555.0819, 559.5991, 564.1525,   &
    568.7429, 573.3701, 578.0399, 582.7611, 587.5379, 592.3701, 597.2572,   &
    602.1997, 607.1975, 612.2507, 617.3596, 622.524, 627.7528, 633.0616,    &
    638.457, 643.9384, 649.505, 655.1568, 660.8936, 666.7153, 672.6219,     &
    678.6133, 684.6945, 690.8745, 697.1569, 703.5416, 710.0284, 716.6172,   &
    723.308, 730.1008, 736.9958, 743.993, 751.0975, 758.3183, 765.6594,     &
    773.1207, 780.702, 788.4033, 796.2249, 804.1667, 812.2289, 820.4118,    &
    828.6444, 828.6444 /)

      integer  :: ntoplw      ! top level to solve for longwave cooling (wrf sets this to 1 for model top below 10 mb)

      logical :: masterproc = .true.
      logical :: ozncyc            ! true => cycle ozone dataset
!     logical :: dosw              ! true => shortwave calculation this timestep
!     logical :: dolw              ! true => longwave calculation this timestep
      logical :: indirect          ! true => include indirect radiative effects of sulfate aerosols
!     logical :: doabsems          ! true => abs/emiss calculation this timestep
      logical :: radforce   = .false.          ! true => calculate aerosol shortwave forcing
      logical :: trace_gas=.false.             ! set true for chemistry
      logical :: strat_volcanic   = .false.    ! true => volcanic aerosol mass available

    real(r8) retab(95)
    !
    !       tabulated values of re(t) in the temperature interval
    !       180 k -- 274 k; hexagonal columns assumed:
    !
    data retab / 						&
         5.92779, 6.26422, 6.61973, 6.99539, 7.39234,	&
         7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,	&
         10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319, 	&
         15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,	&
         20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,	&
         27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943, 	&
         31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601, 	&
         34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,	&
         38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,	&
         42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,	&
         50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,	&
         65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,	&
         93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424, 	&
         124.954, 130.630, 136.457, 142.446, 148.608, 154.956,	&
         161.503, 168.262, 175.248, 182.473, 189.952, 197.699,	&
         205.728, 214.055, 222.694, 231.661, 240.971, 250.639/	
    !
    save retab
contains



subroutine sortarray(n, ain, indxa) 
!-----------------------------------------------
!
! purpose:
!       sort an array
! alogrithm:
!       based on shell's sorting method.
!
! author: t. craig
!-----------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none
!
!  arguments
!
   integer , intent(in) :: n             ! total number of elements
   integer , intent(inout) :: indxa(n)   ! array of integers
   real(r8), intent(inout) :: ain(n)     ! array to sort
!
!  local variables
!
   integer :: i, j                ! loop indices
   integer :: ni                  ! starting increment
   integer :: itmp                ! temporary index
   real(r8):: atmp                ! temporary value to swap
 
   ni = 1 
   do while(.true.) 
      ni = 3*ni + 1 
      if (ni <= n) cycle  
      exit  
   end do 
 
   do while(.true.) 
      ni = ni/3 
      do i = ni + 1, n 
         atmp = ain(i) 
         itmp = indxa(i) 
         j = i 
         do while(.true.) 
            if (ain(j-ni) <= atmp) exit  
            ain(j) = ain(j-ni) 
            indxa(j) = indxa(j-ni) 
            j = j - ni 
            if (j > ni) cycle  
            exit  
         end do 
         ain(j) = atmp 
         indxa(j) = itmp 
      end do 
      if (ni > 1) cycle  
      exit  
   end do 
   return  
 
end subroutine sortarray
subroutine trcab(lchnk   ,ncol    ,pcols, pverp,               &
                 k1      ,k2      ,ucfc11  ,ucfc12  ,un2o0   , &
                 un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
                 uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
                 bch4    ,to3co2  ,pnm     ,dw      ,pnew    , &
                 s2c     ,uptype  ,dplh2o  ,abplnk1 ,tco2    , &
                 th2o    ,to3     ,abstrc  , &
                 aer_trn_ttl)
!----------------------------------------------------------------------- 
! 
! purpose: 
! calculate absorptivity for non nearest layers for ch4, n2o, cfc11 and
! cfc12.
! 
! method: 
! see ccm3 description for equations.
! 
! author: j. kiehl
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use volcrad

   implicit none

!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                    ! chunk identifier
   integer, intent(in) :: ncol                     ! number of atmospheric columns
   integer, intent(in) :: pcols, pverp
   integer, intent(in) :: k1,k2                    ! level indices
!
   real(r8), intent(in) :: to3co2(pcols)           ! pressure weighted temperature
   real(r8), intent(in) :: pnm(pcols,pverp)        ! interface pressures
   real(r8), intent(in) :: ucfc11(pcols,pverp)     ! cfc11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)     ! cfc12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)      ! n2o path length
!
   real(r8), intent(in) :: un2o1(pcols,pverp)      ! n2o path length (hot band)
   real(r8), intent(in) :: uch4(pcols,pverp)       ! ch4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)     ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)     ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)     ! co2 9.4 micron band path length
!
   real(r8), intent(in) :: uco221(pcols,pverp)     ! co2 10.4 micron band path length
   real(r8), intent(in) :: uco222(pcols,pverp)     ! co2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)     ! co2 10.4 micron band path length
   real(r8), intent(in) :: bn2o0(pcols,pverp)      ! pressure factor for n2o
   real(r8), intent(in) :: bn2o1(pcols,pverp)      ! pressure factor for n2o
!
   real(r8), intent(in) :: bch4(pcols,pverp)       ! pressure factor for ch4
   real(r8), intent(in) :: dw(pcols)               ! h2o path length
   real(r8), intent(in) :: pnew(pcols)             ! pressure
   real(r8), intent(in) :: s2c(pcols,pverp)        ! continuum path length
   real(r8), intent(in) :: uptype(pcols,pverp)     ! p-type h2o path length
!
   real(r8), intent(in) :: dplh2o(pcols)           ! p squared h2o path length
   real(r8), intent(in) :: abplnk1(14,pcols,pverp) ! planck factor
   real(r8), intent(in) :: tco2(pcols)             ! co2 transmission factor
   real(r8), intent(in) :: th2o(pcols)             ! h2o transmission factor
   real(r8), intent(in) :: to3(pcols)              ! o3 transmission factor

   real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,bnd_nbr_lw) ! aer trn.

!
!  output arguments
!
   real(r8), intent(out) :: abstrc(pcols)           ! total trace gas absorptivity
!
!--------------------------local variables------------------------------
!
   integer  i,l                     ! loop counters

   real(r8) sqti(pcols)             ! square root of mean temp
   real(r8) du1                     ! cfc11 path length
   real(r8) du2                     ! cfc12 path length
   real(r8) acfc1                   ! cfc11 absorptivity 798 cm-1
   real(r8) acfc2                   ! cfc11 absorptivity 846 cm-1
!
   real(r8) acfc3                   ! cfc11 absorptivity 933 cm-1
   real(r8) acfc4                   ! cfc11 absorptivity 1085 cm-1
   real(r8) acfc5                   ! cfc12 absorptivity 889 cm-1
   real(r8) acfc6                   ! cfc12 absorptivity 923 cm-1
   real(r8) acfc7                   ! cfc12 absorptivity 1102 cm-1
!
   real(r8) acfc8                   ! cfc12 absorptivity 1161 cm-1
   real(r8) du01                    ! n2o path length
   real(r8) dbeta01                 ! n2o pressure factor
   real(r8) dbeta11                 !         "
   real(r8) an2o1                   ! absorptivity of 1285 cm-1 n2o band
!
   real(r8) du02                    ! n2o path length
   real(r8) dbeta02                 ! n2o pressure factor
   real(r8) an2o2                   ! absorptivity of 589 cm-1 n2o band
   real(r8) du03                    ! n2o path length
   real(r8) dbeta03                 ! n2o pressure factor
!
   real(r8) an2o3                   ! absorptivity of 1168 cm-1 n2o band
   real(r8) duch4                   ! ch4 path length
   real(r8) dbetac                  ! ch4 pressure factor
   real(r8) ach4                    ! absorptivity of 1306 cm-1 ch4 band
   real(r8) du11                    ! co2 path length
!
   real(r8) du12                    !       "
   real(r8) du13                    !       "
   real(r8) dbetc1                  ! co2 pressure factor
   real(r8) dbetc2                  ! co2 pressure factor
   real(r8) aco21                   ! absorptivity of 1064 cm-1 band
!
   real(r8) du21                    ! co2 path length
   real(r8) du22                    !       "
   real(r8) du23                    !       "
   real(r8) aco22                   ! absorptivity of 961 cm-1 band
   real(r8) tt(pcols)               ! temp. factor for h2o overlap factor
!
   real(r8) psi1                    !                 "
   real(r8) phi1                    !                 "
   real(r8) p1                      ! h2o overlap factor
   real(r8) w1                      !        "
   real(r8) ds2c(pcols)             ! continuum path length
!
   real(r8) duptyp(pcols)           ! p-type path length
   real(r8) tw(pcols,6)             ! h2o transmission factor
   real(r8) g1(6)                   !         "
   real(r8) g2(6)                   !         "
   real(r8) g3(6)                   !         "
!
   real(r8) g4(6)                   !         "
   real(r8) ab(6)                   ! h2o temp. factor
   real(r8) bb(6)                   !         "
   real(r8) abp(6)                  !         "
   real(r8) bbp(6)                  !         "
!
   real(r8) tcfc3                   ! transmission for cfc11 band
   real(r8) tcfc4                   ! transmission for cfc11 band
   real(r8) tcfc6                   ! transmission for cfc12 band
   real(r8) tcfc7                   ! transmission for cfc12 band
   real(r8) tcfc8                   ! transmission for cfc12 band
!
   real(r8) tlw                     ! h2o transmission
   real(r8) tch4                    ! ch4 transmission
!
!--------------------------data statements------------------------------
!
   data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,0.0321962/
   data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
   data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
   data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,0.0161418/
   data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,2.2915e-2/
   data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,-9.5743e-5,-1.0304e-4/
   data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,1.9892e-2/
   data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,-1.0115e-4,-8.8061e-5/
!
!--------------------------statement functions--------------------------
!
   real(r8) func, u, b
   func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
!
!------------------------------------------------------------------------
!
   do i = 1,ncol
      sqti(i) = sqrt(to3co2(i))
!
! h2o transmission
!
      tt(i) = abs(to3co2(i) - 250.0)
      ds2c(i) = abs(s2c(i,k1) - s2c(i,k2))
      duptyp(i) = abs(uptype(i,k1) - uptype(i,k2))
   end do
!
   do l = 1,6
      do i = 1,ncol
         psi1 = exp(abp(l)*tt(i) + bbp(l)*tt(i)*tt(i))
         phi1 = exp(ab(l)*tt(i) + bb(l)*tt(i)*tt(i))
         p1 = pnew(i)*(psi1/phi1)/sslp
         w1 = dw(i)*phi1
         tw(i,l) = exp(-g1(l)*p1*(sqrt(1.0 + g2(l)*(w1/p1)) - 1.0) - &
                   g3(l)*ds2c(i)-g4(l)*duptyp(i))
      end do
   end do
!
   do i=1,ncol
      tw(i,1)=tw(i,1)*(0.7*aer_trn_ttl(i,k1,k2,idx_lw_0650_0800)+&! l=1: 0750--0820 cm-1
                       0.3*aer_trn_ttl(i,k1,k2,idx_lw_0800_1000)) 
      tw(i,2)=tw(i,2)*aer_trn_ttl(i,k1,k2,idx_lw_0800_1000) ! l=2: 0820--0880 cm-1
      tw(i,3)=tw(i,3)*aer_trn_ttl(i,k1,k2,idx_lw_0800_1000) ! l=3: 0880--0900 cm-1
      tw(i,4)=tw(i,4)*aer_trn_ttl(i,k1,k2,idx_lw_0800_1000) ! l=4: 0900--1000 cm-1
      tw(i,5)=tw(i,5)*aer_trn_ttl(i,k1,k2,idx_lw_1000_1200) ! l=5: 1000--1120 cm-1
      tw(i,6)=tw(i,6)*aer_trn_ttl(i,k1,k2,idx_lw_1000_1200) ! l=6: 1120--1170 cm-1
   end do                    ! end loop over lon
   do i = 1,ncol
      du1 = abs(ucfc11(i,k1) - ucfc11(i,k2))
      du2 = abs(ucfc12(i,k1) - ucfc12(i,k2))
!
! cfc transmissions
!
      tcfc3 = exp(-175.005*du1)
      tcfc4 = exp(-1202.18*du1)
      tcfc6 = exp(-5786.73*du2)
      tcfc7 = exp(-2873.51*du2)
      tcfc8 = exp(-2085.59*du2)
!
! absorptivity for cfc11 bands
!
      acfc1 =  50.0*(1.0 - exp(-54.09*du1))*tw(i,1)*abplnk1(7,i,k2)
      acfc2 =  60.0*(1.0 - exp(-5130.03*du1))*tw(i,2)*abplnk1(8,i,k2)
      acfc3 =  60.0*(1.0 - tcfc3)*tw(i,4)*tcfc6*abplnk1(9,i,k2)
      acfc4 = 100.0*(1.0 - tcfc4)*tw(i,5)*abplnk1(10,i,k2)
!
! absorptivity for cfc12 bands
!
      acfc5 = 45.0*(1.0 - exp(-1272.35*du2))*tw(i,3)*abplnk1(11,i,k2)
      acfc6 = 50.0*(1.0 - tcfc6)* tw(i,4) * abplnk1(12,i,k2)
      acfc7 = 80.0*(1.0 - tcfc7)* tw(i,5) * tcfc4*abplnk1(13,i,k2)
      acfc8 = 70.0*(1.0 - tcfc8)* tw(i,6) * abplnk1(14,i,k2)
!
! emissivity for ch4 band 1306 cm-1
!
      tlw = exp(-1.0*sqrt(dplh2o(i)))
      tlw=tlw*aer_trn_ttl(i,k1,k2,idx_lw_1200_2000)
      duch4 = abs(uch4(i,k1) - uch4(i,k2))
      dbetac = abs(bch4(i,k1) - bch4(i,k2))/duch4
      ach4 = 6.00444*sqti(i)*log(1.0 + func(duch4,dbetac))*tlw*abplnk1(3,i,k2)
      tch4 = 1.0/(1.0 + 0.02*func(duch4,dbetac))
!
! absorptivity for n2o bands
!
      du01 = abs(un2o0(i,k1) - un2o0(i,k2))
      du11 = abs(un2o1(i,k1) - un2o1(i,k2))
      dbeta01 = abs(bn2o0(i,k1) - bn2o0(i,k2))/du01
      dbeta11 = abs(bn2o1(i,k1) - bn2o1(i,k2))/du11
!
! 1285 cm-1 band
!
      an2o1 = 2.35558*sqti(i)*log(1.0 + func(du01,dbeta01) &
              + func(du11,dbeta11))*tlw*tch4*abplnk1(4,i,k2)
      du02 = 0.100090*du01
      du12 = 0.0992746*du11
      dbeta02 = 0.964282*dbeta01
!
! 589 cm-1 band
!
      an2o2 = 2.65581*sqti(i)*log(1.0 + func(du02,dbeta02) + &
              func(du12,dbeta02))*th2o(i)*tco2(i)*abplnk1(5,i,k2)
      du03 = 0.0333767*du01
      dbeta03 = 0.982143*dbeta01
!
! 1168 cm-1 band
!
      an2o3 = 2.54034*sqti(i)*log(1.0 + func(du03,dbeta03))* &
              tw(i,6)*tcfc8*abplnk1(6,i,k2)
!
! emissivity for 1064 cm-1 band of co2
!
      du11 = abs(uco211(i,k1) - uco211(i,k2))
      du12 = abs(uco212(i,k1) - uco212(i,k2))
      du13 = abs(uco213(i,k1) - uco213(i,k2))
      dbetc1 = 2.97558*abs(pnm(i,k1) + pnm(i,k2))/(2.0*sslp*sqti(i))
      dbetc2 = 2.0*dbetc1
      aco21 = 3.7571*sqti(i)*log(1.0 + func(du11,dbetc1) &
              + func(du12,dbetc2) + func(du13,dbetc2)) &
              *to3(i)*tw(i,5)*tcfc4*tcfc7*abplnk1(2,i,k2)
!
! emissivity for 961 cm-1 band
!
      du21 = abs(uco221(i,k1) - uco221(i,k2))
      du22 = abs(uco222(i,k1) - uco222(i,k2))
      du23 = abs(uco223(i,k1) - uco223(i,k2))
      aco22 = 3.8443*sqti(i)*log(1.0 + func(du21,dbetc1) &
              + func(du22,dbetc1) + func(du23,dbetc2)) &
              *tw(i,4)*tcfc3*tcfc6*abplnk1(1,i,k2)
!
! total trace gas absorptivity
!
      abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
                  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 + &
                  aco21 + aco22
   end do
!
   return
!
end subroutine trcab



subroutine trcabn(lchnk   ,ncol    ,pcols, pverp,               &
                  k2      ,kn      ,ucfc11  ,ucfc12  ,un2o0   , &
                  un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
                  uco221  ,uco222  ,uco223  ,tbar    ,bplnk   , &
                  winpl   ,pinpl   ,tco2    ,th2o    ,to3     , &
                  uptype  ,dw      ,s2c     ,up2     ,pnew    , &
                  abstrc  ,uinpl   , &
                  aer_trn_ngh)
!----------------------------------------------------------------------- 
! 
! purpose: 
! calculate nearest layer absorptivity due to ch4, n2o, cfc11 and cfc12
! 
! method: 
! equations in ccm3 description
! 
! author: j. kiehl
! 
!-----------------------------------------------------------------------
!
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use volcrad

   implicit none
 
!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: pcols, pverp
   integer, intent(in) :: k2                    ! level index
   integer, intent(in) :: kn                    ! level index
!
   real(r8), intent(in) :: tbar(pcols,4)        ! pressure weighted temperature
   real(r8), intent(in) :: ucfc11(pcols,pverp)  ! cfc11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)  ! cfc12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)   ! n2o path length
   real(r8), intent(in) :: un2o1(pcols,pverp)   ! n2o path length (hot band)
!
   real(r8), intent(in) :: uch4(pcols,pverp)    ! ch4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco221(pcols,pverp)  ! co2 10.4 micron band path length
!
   real(r8), intent(in) :: uco222(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8), intent(in) :: bplnk(14,pcols,4)    ! weighted planck fnc. for absorptivity
   real(r8), intent(in) :: winpl(pcols,4)       ! fractional path length
   real(r8), intent(in) :: pinpl(pcols,4)       ! pressure factor for subdivided layer
!
   real(r8), intent(in) :: tco2(pcols)          ! co2 transmission
   real(r8), intent(in) :: th2o(pcols)          ! h2o transmission
   real(r8), intent(in) :: to3(pcols)           ! o3 transmission
   real(r8), intent(in) :: dw(pcols)            ! h2o path length
   real(r8), intent(in) :: pnew(pcols)          ! pressure factor
!
   real(r8), intent(in) :: s2c(pcols,pverp)     ! h2o continuum factor
   real(r8), intent(in) :: uptype(pcols,pverp)  ! p-type path length
   real(r8), intent(in) :: up2(pcols)           ! p squared path length
   real(r8), intent(in) :: uinpl(pcols,4)       ! nearest layer subdivision factor
   real(r8), intent(in) :: aer_trn_ngh(pcols,bnd_nbr_lw) 
                             ! [fraction] total transmission between 
                             !            nearest neighbor sub-levels
!
!  output arguments
!
   real(r8), intent(out) :: abstrc(pcols)        ! total trace gas absorptivity

!
!--------------------------local variables------------------------------
!
   integer i,l                   ! loop counters
!
   real(r8) sqti(pcols)          ! square root of mean temp
   real(r8) rsqti(pcols)         ! reciprocal of sqti
   real(r8) du1                  ! cfc11 path length
   real(r8) du2                  ! cfc12 path length
   real(r8) acfc1                ! absorptivity of cfc11 798 cm-1 band
!
   real(r8) acfc2                ! absorptivity of cfc11 846 cm-1 band
   real(r8) acfc3                ! absorptivity of cfc11 933 cm-1 band
   real(r8) acfc4                ! absorptivity of cfc11 1085 cm-1 band
   real(r8) acfc5                ! absorptivity of cfc11 889 cm-1 band
   real(r8) acfc6                ! absorptivity of cfc11 923 cm-1 band
!
   real(r8) acfc7                ! absorptivity of cfc11 1102 cm-1 band
   real(r8) acfc8                ! absorptivity of cfc11 1161 cm-1 band
   real(r8) du01                 ! n2o path length
   real(r8) dbeta01              ! n2o pressure factors
   real(r8) dbeta11              !        "
!
   real(r8)  an2o1               ! absorptivity of the 1285 cm-1 n2o band
   real(r8) du02                 ! n2o path length
   real(r8) dbeta02              ! n2o pressure factor
   real(r8) an2o2                ! absorptivity of the 589 cm-1 n2o band
   real(r8) du03                 ! n2o path length
!
   real(r8) dbeta03              ! n2o pressure factor
   real(r8) an2o3                ! absorptivity of the 1168 cm-1 n2o band
   real(r8) duch4                ! ch4 path length
   real(r8) dbetac               ! ch4 pressure factor
   real(r8) ach4                 ! absorptivity of the 1306 cm-1 ch4 band
!
   real(r8) du11                 ! co2 path length
   real(r8) du12                 !       "
   real(r8) du13                 !       "
   real(r8) dbetc1               ! co2 pressure factor
   real(r8) dbetc2               ! co2 pressure factor
!
   real(r8) aco21                ! absorptivity of the 1064 cm-1 co2 band
   real(r8) du21                 ! co2 path length
   real(r8) du22                 !       "
   real(r8) du23                 !       "
   real(r8) aco22                ! absorptivity of the 961 cm-1 co2 band
!
   real(r8) tt(pcols)            ! temp. factor for h2o overlap
   real(r8) psi1                 !          "
   real(r8) phi1                 !          "
   real(r8) p1                   ! factor for h2o overlap
   real(r8) w1                   !          "
!
   real(r8) ds2c(pcols)          ! continuum path length
   real(r8) duptyp(pcols)        ! p-type path length
   real(r8) tw(pcols,6)          ! h2o transmission overlap
   real(r8) g1(6)                ! h2o overlap factor
   real(r8) g2(6)                !         "
!
   real(r8) g3(6)                !         "
   real(r8) g4(6)                !         "
   real(r8) ab(6)                ! h2o temp. factor
   real(r8) bb(6)                !         "
   real(r8) abp(6)               !         "
!
   real(r8) bbp(6)               !         "
   real(r8) tcfc3                ! transmission of cfc11 band
   real(r8) tcfc4                ! transmission of cfc11 band
   real(r8) tcfc6                ! transmission of cfc12 band
   real(r8) tcfc7                !         "
!
   real(r8) tcfc8                !         "
   real(r8) tlw                  ! h2o transmission
   real(r8) tch4                 ! ch4 transmission
!
!--------------------------data statements------------------------------
!
   data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,0.0321962/
   data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
   data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
   data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,0.0161418/
   data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,2.2915e-2/
   data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,-9.5743e-5,-1.0304e-4/
   data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,1.9892e-2/
   data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,-1.0115e-4,-8.8061e-5/
!
!--------------------------statement functions--------------------------
!
   real(r8) func, u, b
   func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
!
!------------------------------------------------------------------
!
   do i = 1,ncol
      sqti(i) = sqrt(tbar(i,kn))
      rsqti(i) = 1. / sqti(i)
!
! h2o transmission
!
      tt(i) = abs(tbar(i,kn) - 250.0)
      ds2c(i) = abs(s2c(i,k2+1) - s2c(i,k2))*uinpl(i,kn)
      duptyp(i) = abs(uptype(i,k2+1) - uptype(i,k2))*uinpl(i,kn)
   end do
!
   do l = 1,6
      do i = 1,ncol
         psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
         phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
         p1 = pnew(i) * (psi1/phi1) / sslp
         w1 = dw(i) * winpl(i,kn) * phi1
         tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0+g2(l)*(w1/p1))-1.0) &
                   - g3(l)*ds2c(i)-g4(l)*duptyp(i))
      end do
   end do
!
   do i=1,ncol
      tw(i,1)=tw(i,1)*(0.7*aer_trn_ngh(i,idx_lw_0650_0800)+&! l=1: 0750--0820 cm-1
                       0.3*aer_trn_ngh(i,idx_lw_0800_1000))
      tw(i,2)=tw(i,2)*aer_trn_ngh(i,idx_lw_0800_1000) ! l=2: 0820--0880 cm-1
      tw(i,3)=tw(i,3)*aer_trn_ngh(i,idx_lw_0800_1000) ! l=3: 0880--0900 cm-1
      tw(i,4)=tw(i,4)*aer_trn_ngh(i,idx_lw_0800_1000) ! l=4: 0900--1000 cm-1
      tw(i,5)=tw(i,5)*aer_trn_ngh(i,idx_lw_1000_1200) ! l=5: 1000--1120 cm-1
      tw(i,6)=tw(i,6)*aer_trn_ngh(i,idx_lw_1000_1200) ! l=6: 1120--1170 cm-1
   end do                    ! end loop over lon

   do i = 1,ncol
!
      du1 = abs(ucfc11(i,k2+1) - ucfc11(i,k2)) * winpl(i,kn)
      du2 = abs(ucfc12(i,k2+1) - ucfc12(i,k2)) * winpl(i,kn)
!
! cfc transmissions
!
      tcfc3 = exp(-175.005*du1)
      tcfc4 = exp(-1202.18*du1)
      tcfc6 = exp(-5786.73*du2)
      tcfc7 = exp(-2873.51*du2)
      tcfc8 = exp(-2085.59*du2)
!
! absorptivity for cfc11 bands
!
      acfc1 = 50.0*(1.0 - exp(-54.09*du1)) * tw(i,1)*bplnk(7,i,kn)
      acfc2 = 60.0*(1.0 - exp(-5130.03*du1))*tw(i,2)*bplnk(8,i,kn)
      acfc3 = 60.0*(1.0 - tcfc3)*tw(i,4)*tcfc6 * bplnk(9,i,kn)
      acfc4 = 100.0*(1.0 - tcfc4)* tw(i,5) * bplnk(10,i,kn)
!
! absorptivity for cfc12 bands
!
      acfc5 = 45.0*(1.0 - exp(-1272.35*du2))*tw(i,3)*bplnk(11,i,kn)
      acfc6 = 50.0*(1.0 - tcfc6)*tw(i,4)*bplnk(12,i,kn)
      acfc7 = 80.0*(1.0 - tcfc7)* tw(i,5)*tcfc4 *bplnk(13,i,kn)
      acfc8 = 70.0*(1.0 - tcfc8)*tw(i,6)*bplnk(14,i,kn)
!
! absorptivity for ch4 band 1306 cm-1
!
      tlw = exp(-1.0*sqrt(up2(i)))
      tlw=tlw*aer_trn_ngh(i,idx_lw_1200_2000)
      duch4 = abs(uch4(i,k2+1) - uch4(i,k2)) * winpl(i,kn)
      dbetac = 2.94449 * pinpl(i,kn) * rsqti(i) / sslp
      ach4 = 6.00444*sqti(i)*log(1.0 + func(duch4,dbetac)) * tlw * bplnk(3,i,kn)
      tch4 = 1.0/(1.0 + 0.02*func(duch4,dbetac))
!
! absorptivity for n2o bands
!
      du01 = abs(un2o0(i,k2+1) - un2o0(i,k2)) * winpl(i,kn)
      du11 = abs(un2o1(i,k2+1) - un2o1(i,k2)) * winpl(i,kn)
      dbeta01 = 19.399 *  pinpl(i,kn) * rsqti(i) / sslp
      dbeta11 = dbeta01
!
! 1285 cm-1 band
!
      an2o1 = 2.35558*sqti(i)*log(1.0 + func(du01,dbeta01) &
              + func(du11,dbeta11)) * tlw * tch4 * bplnk(4,i,kn)
      du02 = 0.100090*du01
      du12 = 0.0992746*du11
      dbeta02 = 0.964282*dbeta01
!
! 589 cm-1 band
!
      an2o2 = 2.65581*sqti(i)*log(1.0 + func(du02,dbeta02) &
              +  func(du12,dbeta02)) * tco2(i) * th2o(i) * bplnk(5,i,kn)
      du03 = 0.0333767*du01
      dbeta03 = 0.982143*dbeta01
!
! 1168 cm-1 band
!
      an2o3 = 2.54034*sqti(i)*log(1.0 + func(du03,dbeta03)) * &
              tw(i,6) * tcfc8 * bplnk(6,i,kn)
!
! absorptivity for 1064 cm-1 band of co2
!
      du11 = abs(uco211(i,k2+1) - uco211(i,k2)) * winpl(i,kn)
      du12 = abs(uco212(i,k2+1) - uco212(i,k2)) * winpl(i,kn)
      du13 = abs(uco213(i,k2+1) - uco213(i,k2)) * winpl(i,kn)
      dbetc1 = 2.97558 * pinpl(i,kn) * rsqti(i) / sslp
      dbetc2 = 2.0 * dbetc1
      aco21 = 3.7571*sqti(i)*log(1.0 + func(du11,dbetc1) &
              + func(du12,dbetc2) + func(du13,dbetc2)) &
              * to3(i) * tw(i,5) * tcfc4 * tcfc7 * bplnk(2,i,kn)
!
! absorptivity for 961 cm-1 band of co2
!
      du21 = abs(uco221(i,k2+1) - uco221(i,k2)) * winpl(i,kn)
      du22 = abs(uco222(i,k2+1) - uco222(i,k2)) * winpl(i,kn)
      du23 = abs(uco223(i,k2+1) - uco223(i,k2)) * winpl(i,kn)
      aco22 = 3.8443*sqti(i)*log(1.0 + func(du21,dbetc1) &
              + func(du22,dbetc1) + func(du23,dbetc2)) &
              * tw(i,4) * tcfc3 * tcfc6 * bplnk(1,i,kn)
!
! total trace gas absorptivity
!
      abstrc(i) = acfc1 + acfc2 + acfc3 + acfc4 + acfc5 + acfc6 + &
                  acfc7 + acfc8 + an2o1 + an2o2 + an2o3 + ach4 + &
                  aco21 + aco22
   end do
!
   return
!
end subroutine trcabn



subroutine trcems(lchnk   ,ncol    ,pcols, pverp,               &
                  k       ,co2t    ,pnm     ,ucfc11  ,ucfc12  , &
                  un2o0   ,un2o1   ,bn2o0   ,bn2o1   ,uch4    , &
                  bch4    ,uco211  ,uco212  ,uco213  ,uco221  , &
                  uco222  ,uco223  ,uptype  ,w       ,s2c     , &
                  up2     ,emplnk  ,th2o    ,tco2    ,to3     , &
                  emstrc  , &
                 aer_trn_ttl)
!----------------------------------------------------------------------- 
! 
! purpose: 
!  calculate emissivity for ch4, n2o, cfc11 and cfc12 bands.
! 
! method: 
!  see ccm3 description for equations.
! 
! author: j. kiehl
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use volcrad

   implicit none

!
!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: pcols, pverp

   real(r8), intent(in) :: co2t(pcols,pverp)    ! pressure weighted temperature
   real(r8), intent(in) :: pnm(pcols,pverp)     ! interface pressure
   real(r8), intent(in) :: ucfc11(pcols,pverp)  ! cfc11 path length
   real(r8), intent(in) :: ucfc12(pcols,pverp)  ! cfc12 path length
   real(r8), intent(in) :: un2o0(pcols,pverp)   ! n2o path length
!
   real(r8), intent(in) :: un2o1(pcols,pverp)   ! n2o path length (hot band)
   real(r8), intent(in) :: uch4(pcols,pverp)    ! ch4 path length
   real(r8), intent(in) :: uco211(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco212(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8), intent(in) :: uco213(pcols,pverp)  ! co2 9.4 micron band path length
!
   real(r8), intent(in) :: uco221(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8), intent(in) :: uco222(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8), intent(in) :: uco223(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8), intent(in) :: uptype(pcols,pverp)  ! continuum path length
   real(r8), intent(in) :: bn2o0(pcols,pverp)   ! pressure factor for n2o
!
   real(r8), intent(in) :: bn2o1(pcols,pverp)   ! pressure factor for n2o
   real(r8), intent(in) :: bch4(pcols,pverp)    ! pressure factor for ch4
   real(r8), intent(in) :: emplnk(14,pcols)     ! emissivity planck factor
   real(r8), intent(in) :: th2o(pcols)          ! water vapor overlap factor
   real(r8), intent(in) :: tco2(pcols)          ! co2 overlap factor
!
   real(r8), intent(in) :: to3(pcols)           ! o3 overlap factor
   real(r8), intent(in) :: s2c(pcols,pverp)     ! h2o continuum path length
   real(r8), intent(in) :: w(pcols,pverp)       ! h2o path length
   real(r8), intent(in) :: up2(pcols)           ! pressure squared h2o path length
!
   integer, intent(in) :: k                 ! level index

   real(r8), intent(in) :: aer_trn_ttl(pcols,pverp,pverp,bnd_nbr_lw) ! aer trn.

!
!  output arguments
!
   real(r8), intent(out) :: emstrc(pcols,pverp)  ! total trace gas emissivity

!
!--------------------------local variables------------------------------
!
   integer i,l               ! loop counters
!
   real(r8) sqti(pcols)          ! square root of mean temp
   real(r8) ecfc1                ! emissivity of cfc11 798 cm-1 band
   real(r8) ecfc2                !     "      "    "   846 cm-1 band
   real(r8) ecfc3                !     "      "    "   933 cm-1 band
   real(r8) ecfc4                !     "      "    "   1085 cm-1 band
!
   real(r8) ecfc5                !     "      "  cfc12 889 cm-1 band
   real(r8) ecfc6                !     "      "    "   923 cm-1 band
   real(r8) ecfc7                !     "      "    "   1102 cm-1 band
   real(r8) ecfc8                !     "      "    "   1161 cm-1 band
   real(r8) u01                  ! n2o path length
!
   real(r8) u11                  ! n2o path length
   real(r8) beta01               ! n2o pressure factor
   real(r8) beta11               ! n2o pressure factor
   real(r8) en2o1                ! emissivity of the 1285 cm-1 n2o band
   real(r8) u02                  ! n2o path length
!
   real(r8) u12                  ! n2o path length
   real(r8) beta02               ! n2o pressure factor
   real(r8) en2o2                ! emissivity of the 589 cm-1 n2o band
   real(r8) u03                  ! n2o path length
   real(r8) beta03               ! n2o pressure factor
!
   real(r8) en2o3                ! emissivity of the 1168 cm-1 n2o band
   real(r8) betac                ! ch4 pressure factor
   real(r8) ech4                 ! emissivity of 1306 cm-1 ch4 band
   real(r8) betac1               ! co2 pressure factor
   real(r8) betac2               ! co2 pressure factor
!
   real(r8) eco21                ! emissivity of 1064 cm-1 co2 band
   real(r8) eco22                ! emissivity of 961 cm-1 co2 band
   real(r8) tt(pcols)            ! temp. factor for h2o overlap factor
   real(r8) psi1                 ! narrow band h2o temp. factor
   real(r8) phi1                 !             "
!
   real(r8) p1                   ! h2o line overlap factor
   real(r8) w1                   !          "
   real(r8) tw(pcols,6)          ! h2o transmission overlap
   real(r8) g1(6)                ! h2o overlap factor
   real(r8) g2(6)                !          "
!
   real(r8) g3(6)                !          "
   real(r8) g4(6)                !          "
   real(r8) ab(6)                !          "
   real(r8) bb(6)                !          "
   real(r8) abp(6)               !          "
!
   real(r8) bbp(6)               !          "
   real(r8) tcfc3                ! transmission for cfc11 band
   real(r8) tcfc4                !          "
   real(r8) tcfc6                ! transmission for cfc12 band
   real(r8) tcfc7                !          "
!
   real(r8) tcfc8                !          "
   real(r8) tlw                  ! h2o overlap factor
   real(r8) tch4                 ! ch4 overlap factor
!
!--------------------------data statements------------------------------
!
   data g1 /0.0468556,0.0397454,0.0407664,0.0304380,0.0540398,0.0321962/
   data g2 /14.4832,4.30242,5.23523,3.25342,0.698935,16.5599/
   data g3 /26.1898,18.4476,15.3633,12.1927,9.14992,8.07092/
   data g4 /0.0261782,0.0369516,0.0307266,0.0243854,0.0182932,0.0161418/
   data ab /3.0857e-2,2.3524e-2,1.7310e-2,2.6661e-2,2.8074e-2,2.2915e-2/
   data bb /-1.3512e-4,-6.8320e-5,-3.2609e-5,-1.0228e-5,-9.5743e-5,-1.0304e-4/
   data abp/2.9129e-2,2.4101e-2,1.9821e-2,2.6904e-2,2.9458e-2,1.9892e-2/
   data bbp/-1.3139e-4,-5.5688e-5,-4.6380e-5,-8.0362e-5,-1.0115e-4,-8.8061e-5/
!
!--------------------------statement functions--------------------------
!
   real(r8) func, u, b
   func(u,b) = u/sqrt(4.0 + u*(1.0 + 1.0 / b))
!
!-----------------------------------------------------------------------
!
   do i = 1,ncol
      sqti(i) = sqrt(co2t(i,k))
!
! transmission for h2o
!
      tt(i) = abs(co2t(i,k) - 250.0)
   end do
!
   do l = 1,6
      do i = 1,ncol
         psi1 = exp(abp(l)*tt(i)+bbp(l)*tt(i)*tt(i))
         phi1 = exp(ab(l)*tt(i)+bb(l)*tt(i)*tt(i))
         p1 = pnm(i,k) * (psi1/phi1) / sslp
         w1 = w(i,k) * phi1
         tw(i,l) = exp(- g1(l)*p1*(sqrt(1.0+g2(l)*(w1/p1))-1.0) &
                   - g3(l)*s2c(i,k)-g4(l)*uptype(i,k))
      end do
   end do

!     overlap h2o tranmission with straer continuum in 6 trace gas 
!                 subbands

      do i=1,ncol
         tw(i,1)=tw(i,1)*(0.7*aer_trn_ttl(i,k,1,idx_lw_0650_0800)+&! l=1: 0750--0820 cm-1
                          0.3*aer_trn_ttl(i,k,1,idx_lw_0800_1000))
         tw(i,2)=tw(i,2)*aer_trn_ttl(i,k,1,idx_lw_0800_1000) ! l=2: 0820--0880 cm-1
         tw(i,3)=tw(i,3)*aer_trn_ttl(i,k,1,idx_lw_0800_1000) ! l=3: 0880--0900 cm-1
         tw(i,4)=tw(i,4)*aer_trn_ttl(i,k,1,idx_lw_0800_1000) ! l=4: 0900--1000 cm-1
         tw(i,5)=tw(i,5)*aer_trn_ttl(i,k,1,idx_lw_1000_1200) ! l=5: 1000--1120 cm-1
         tw(i,6)=tw(i,6)*aer_trn_ttl(i,k,1,idx_lw_1000_1200) ! l=6: 1120--1170 cm-1
      end do                    ! end loop over lon
!
   do i = 1,ncol
!
! transmission due to cfc bands
!
      tcfc3 = exp(-175.005*ucfc11(i,k))
      tcfc4 = exp(-1202.18*ucfc11(i,k))
      tcfc6 = exp(-5786.73*ucfc12(i,k))
      tcfc7 = exp(-2873.51*ucfc12(i,k))
      tcfc8 = exp(-2085.59*ucfc12(i,k))
!
! emissivity for cfc11 bands
!
      ecfc1 = 50.0*(1.0 - exp(-54.09*ucfc11(i,k))) * tw(i,1) * emplnk(7,i)
      ecfc2 = 60.0*(1.0 - exp(-5130.03*ucfc11(i,k)))* tw(i,2) * emplnk(8,i)
      ecfc3 = 60.0*(1.0 - tcfc3)*tw(i,4)*tcfc6*emplnk(9,i)
      ecfc4 = 100.0*(1.0 - tcfc4)*tw(i,5)*emplnk(10,i)
!
! emissivity for cfc12 bands
!
      ecfc5 = 45.0*(1.0 - exp(-1272.35*ucfc12(i,k)))*tw(i,3)*emplnk(11,i)
      ecfc6 = 50.0*(1.0 - tcfc6)*tw(i,4)*emplnk(12,i)
      ecfc7 = 80.0*(1.0 - tcfc7)*tw(i,5)* tcfc4 * emplnk(13,i)
      ecfc8 = 70.0*(1.0 - tcfc8)*tw(i,6) * emplnk(14,i)
!
! emissivity for ch4 band 1306 cm-1
!
      tlw = exp(-1.0*sqrt(up2(i)))

!     overlap h2o vibration rotation band with straer continuum 
!             for ch4 1306 cm-1 and n2o 1285 cm-1 bands

            tlw=tlw*aer_trn_ttl(i,k,1,idx_lw_1200_2000)
      betac = bch4(i,k)/uch4(i,k)
      ech4 = 6.00444*sqti(i)*log(1.0 + func(uch4(i,k),betac)) *tlw * emplnk(3,i)
      tch4 = 1.0/(1.0 + 0.02*func(uch4(i,k),betac))
!
! emissivity for n2o bands
!
      u01 = un2o0(i,k)
      u11 = un2o1(i,k)
      beta01 = bn2o0(i,k)/un2o0(i,k)
      beta11 = bn2o1(i,k)/un2o1(i,k)
!
! 1285 cm-1 band
!
      en2o1 = 2.35558*sqti(i)*log(1.0 + func(u01,beta01) + &
              func(u11,beta11))*tlw*tch4*emplnk(4,i)
      u02 = 0.100090*u01
      u12 = 0.0992746*u11
      beta02 = 0.964282*beta01
!
! 589 cm-1 band
!
      en2o2 = 2.65581*sqti(i)*log(1.0 + func(u02,beta02) + &
              func(u12,beta02)) * tco2(i) * th2o(i) * emplnk(5,i)
      u03 = 0.0333767*u01
      beta03 = 0.982143*beta01
!
! 1168 cm-1 band
!
      en2o3 = 2.54034*sqti(i)*log(1.0 + func(u03,beta03)) * &
              tw(i,6) * tcfc8 * emplnk(6,i)
!
! emissivity for 1064 cm-1 band of co2
!
      betac1 = 2.97558*pnm(i,k) / (sslp*sqti(i))
      betac2 = 2.0 * betac1
      eco21 = 3.7571*sqti(i)*log(1.0 + func(uco211(i,k),betac1) &
              + func(uco212(i,k),betac2) + func(uco213(i,k),betac2)) &
              * to3(i) * tw(i,5) * tcfc4 * tcfc7 * emplnk(2,i)
!
! emissivity for 961 cm-1 band
!
      eco22 = 3.8443*sqti(i)*log(1.0 + func(uco221(i,k),betac1) &
              + func(uco222(i,k),betac1) + func(uco223(i,k),betac2)) &
              * tw(i,4) * tcfc3 * tcfc6 * emplnk(1,i)
!
! total trace gas emissivity
!
      emstrc(i,k) = ecfc1 + ecfc2 + ecfc3 + ecfc4 + ecfc5 +ecfc6 + &
                    ecfc7 + ecfc8 + en2o1 + en2o2 + en2o3 + ech4 + &
                    eco21 + eco22
   end do
!
   return
!
end subroutine trcems

subroutine trcmix(lchnk   ,ncol     ,pcols, pver, &
                  pmid    ,clat, n2o      ,ch4     ,          &
                  cfc11   , cfc12   )
!----------------------------------------------------------------------- 
! 
! purpose: 
! specify zonal mean mass mixing ratios of ch4, n2o, cfc11 and
! cfc12
! 
! method: 
! distributions assume constant mixing ratio in the troposphere
! and a decrease of mixing ratio in the stratosphere. tropopause
! defined by ptrop. the scale height of the particular trace gas
! depends on latitude. this assumption produces a more realistic
! stratospheric distribution of the various trace gases.
! 
! author: j. kiehl
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use phys_grid,    only: get_rlat_all_p
!  use physconst,    only: mwdry, mwch4, mwn2o, mwf11, mwf12
!  use ghg_surfvals, only: ch4vmr, n2ovmr, f11vmr, f12vmr

   implicit none

!-----------------------------arguments---------------------------------
!
! input
!
   integer, intent(in) :: lchnk                    ! chunk identifier
   integer, intent(in) :: ncol                     ! number of atmospheric columns
   integer, intent(in) :: pcols, pver

   real(r8), intent(in) :: pmid(pcols,pver)        ! model pressures
   real(r8), intent(in) :: clat(pcols)             ! latitude in radians for columns
!
! output
!
   real(r8), intent(out) :: n2o(pcols,pver)         ! nitrous oxide mass mixing ratio
   real(r8), intent(out) :: ch4(pcols,pver)         ! methane mass mixing ratio
   real(r8), intent(out) :: cfc11(pcols,pver)       ! cfc11 mass mixing ratio
   real(r8), intent(out) :: cfc12(pcols,pver)       ! cfc12 mass mixing ratio

!
!--------------------------local variables------------------------------

   real(r8) :: rmwn2o       ! ratio of molecular weight n2o   to dry air
   real(r8) :: rmwch4       ! ratio of molecular weight ch4   to dry air
   real(r8) :: rmwf11       ! ratio of molecular weight cfc11 to dry air
   real(r8) :: rmwf12       ! ratio of molecular weight cfc12 to dry air
!
   integer i                ! longitude loop index
   integer k                ! level index
!
!  real(r8) clat(pcols)         ! latitude in radians for columns
   real(r8) coslat(pcols)       ! cosine of latitude
   real(r8) dlat                ! latitude in degrees
   real(r8) ptrop               ! pressure level of tropopause
   real(r8) pratio              ! pressure divided by ptrop
!
   real(r8) xn2o                ! pressure scale height for n2o
   real(r8) xch4                ! pressure scale height for ch4
   real(r8) xcfc11              ! pressure scale height for cfc11
   real(r8) xcfc12              ! pressure scale height for cfc12
!
   real(r8) ch40                ! tropospheric mass mixing ratio for ch4
   real(r8) n2o0                ! tropospheric mass mixing ratio for n2o
   real(r8) cfc110              ! tropospheric mass mixing ratio for cfc11
   real(r8) cfc120              ! tropospheric mass mixing ratio for cfc12
!
!-----------------------------------------------------------------------
   rmwn2o = mwn2o/mwdry      ! ratio of molecular weight n2o   to dry air
   rmwch4 = mwch4/mwdry      ! ratio of molecular weight ch4   to dry air
   rmwf11 = mwf11/mwdry      ! ratio of molecular weight cfc11 to dry air
   rmwf12 = mwf12/mwdry      ! ratio of molecular weight cfc12 to dry air
!
! get latitudes
!
!  call get_rlat_all_p(lchnk, ncol, clat)
   do i = 1, ncol
      coslat(i) = cos(clat(i))
   end do
!
! set tropospheric mass mixing ratios
!
   ch40   = rmwch4 * ch4vmr
   n2o0   = rmwn2o * n2ovmr
   cfc110 = rmwf11 * f11vmr
   cfc120 = rmwf12 * f12vmr

   do i = 1, ncol
      coslat(i) = cos(clat(i))
   end do
!
   do k = 1,pver
      do i = 1,ncol
!
!        set stratospheric scale height factor for gases
         dlat = abs(57.2958 * clat(i))
         if(dlat.le.45.0) then
            xn2o = 0.3478 + 0.00116 * dlat
            xch4 = 0.2353
            xcfc11 = 0.7273 + 0.00606 * dlat
            xcfc12 = 0.4000 + 0.00222 * dlat
         else
            xn2o = 0.4000 + 0.013333 * (dlat - 45)
            xch4 = 0.2353 + 0.0225489 * (dlat - 45)
            xcfc11 = 1.00 + 0.013333 * (dlat - 45)
            xcfc12 = 0.50 + 0.024444 * (dlat - 45)
         end if
!
!        pressure of tropopause
         ptrop = 250.0e2 - 150.0e2*coslat(i)**2.0
!
!        determine output mass mixing ratios
         if (pmid(i,k) >= ptrop) then
            ch4(i,k) = ch40
            n2o(i,k) = n2o0
            cfc11(i,k) = cfc110
            cfc12(i,k) = cfc120
         else
            pratio = pmid(i,k)/ptrop
            ch4(i,k) = ch40 * (pratio)**xch4
            n2o(i,k) = n2o0 * (pratio)**xn2o
            cfc11(i,k) = cfc110 * (pratio)**xcfc11
            cfc12(i,k) = cfc120 * (pratio)**xcfc12
         end if
      end do
   end do
!
   return
!
end subroutine trcmix

subroutine trcplk(lchnk   ,ncol    ,pcols, pver, pverp,         &
                  tint    ,tlayr   ,tplnke  ,emplnk  ,abplnk1 , &
                  abplnk2 )
!----------------------------------------------------------------------- 
! 
! purpose: 
!   calculate planck factors for absorptivity and emissivity of
!   ch4, n2o, cfc11 and cfc12
! 
! method: 
!   planck function and derivative evaluated at the band center.
! 
! author: j. kiehl
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid

   implicit none
!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                ! chunk identifier
   integer, intent(in) :: ncol                 ! number of atmospheric columns
   integer, intent(in) :: pcols, pver, pverp

   real(r8), intent(in) :: tint(pcols,pverp)   ! interface temperatures
   real(r8), intent(in) :: tlayr(pcols,pverp)  ! k-1 level temperatures
   real(r8), intent(in) :: tplnke(pcols)       ! top layer temperature
!
! output arguments
!
   real(r8), intent(out) :: emplnk(14,pcols)         ! emissivity planck factor
   real(r8), intent(out) :: abplnk1(14,pcols,pverp)  ! non-nearest layer plack factor
   real(r8), intent(out) :: abplnk2(14,pcols,pverp)  ! nearest layer factor

!
!--------------------------local variables------------------------------
!
   integer wvl                   ! wavelength index
   integer i,k                   ! loop counters
!
   real(r8) f1(14)                   ! planck function factor
   real(r8) f2(14)                   !        "
   real(r8) f3(14)                   !        "
!
!--------------------------data statements------------------------------
!
   data f1 /5.85713e8,7.94950e8,1.47009e9,1.40031e9,1.34853e8, &
            1.05158e9,3.35370e8,3.99601e8,5.35994e8,8.42955e8, &
            4.63682e8,5.18944e8,8.83202e8,1.03279e9/
   data f2 /2.02493e11,3.04286e11,6.90698e11,6.47333e11, &
            2.85744e10,4.41862e11,9.62780e10,1.21618e11, &
            1.79905e11,3.29029e11,1.48294e11,1.72315e11, &
            3.50140e11,4.31364e11/
   data f3 /1383.0,1531.0,1879.0,1849.0,848.0,1681.0, &
            1148.0,1217.0,1343.0,1561.0,1279.0,1328.0, &
            1586.0,1671.0/
!
!-----------------------------------------------------------------------
!
! calculate emissivity planck factor
!
   do wvl = 1,14
      do i = 1,ncol
         emplnk(wvl,i) = f1(wvl)/(tplnke(i)**4.0*(exp(f3(wvl)/tplnke(i))-1.0))
      end do
   end do
!
! calculate absorptivity planck factor for tint and tlayr temperatures
!
   do wvl = 1,14
      do k = ntoplw, pverp
         do i = 1, ncol
!
! non-nearlest layer function
!
            abplnk1(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tint(i,k)))  &
                               /(tint(i,k)**5.0*(exp(f3(wvl)/tint(i,k))-1.0)**2.0)
!
! nearest layer function
!
            abplnk2(wvl,i,k) = (f2(wvl)*exp(f3(wvl)/tlayr(i,k))) &
                               /(tlayr(i,k)**5.0*(exp(f3(wvl)/tlayr(i,k))-1.0)**2.0)
         end do
      end do
   end do
!
   return
end subroutine trcplk

subroutine trcpth(lchnk   ,ncol    ,pcols, pver, pverp,         &
                  tnm     ,pnm     ,cfc11   ,cfc12   ,n2o     , &
                  ch4     ,qnm     ,ucfc11  ,ucfc12  ,un2o0   , &
                  un2o1   ,uch4    ,uco211  ,uco212  ,uco213  , &
                  uco221  ,uco222  ,uco223  ,bn2o0   ,bn2o1   , &
                  bch4    ,uptype  )
!----------------------------------------------------------------------- 
! 
! purpose: 
! calculate path lengths and pressure factors for ch4, n2o, cfc11
! and cfc12.
! 
! method: 
! see ccm3 description for details
! 
! author: j. kiehl
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use ghg_surfvals, only: co2mmr

   implicit none

!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: pcols, pver, pverp

   real(r8), intent(in) :: tnm(pcols,pver)      ! model level temperatures
   real(r8), intent(in) :: pnm(pcols,pverp)     ! pres. at model interfaces (dynes/cm2)
   real(r8), intent(in) :: qnm(pcols,pver)      ! h2o specific humidity
   real(r8), intent(in) :: cfc11(pcols,pver)    ! cfc11 mass mixing ratio
!
   real(r8), intent(in) :: cfc12(pcols,pver)    ! cfc12 mass mixing ratio
   real(r8), intent(in) :: n2o(pcols,pver)      ! n2o mass mixing ratio
   real(r8), intent(in) :: ch4(pcols,pver)      ! ch4 mass mixing ratio

!
! output arguments
!
   real(r8), intent(out) :: ucfc11(pcols,pverp)  ! cfc11 path length
   real(r8), intent(out) :: ucfc12(pcols,pverp)  ! cfc12 path length
   real(r8), intent(out) :: un2o0(pcols,pverp)   ! n2o path length
   real(r8), intent(out) :: un2o1(pcols,pverp)   ! n2o path length (hot band)
   real(r8), intent(out) :: uch4(pcols,pverp)    ! ch4 path length
!
   real(r8), intent(out) :: uco211(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8), intent(out) :: uco212(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8), intent(out) :: uco213(pcols,pverp)  ! co2 9.4 micron band path length
   real(r8), intent(out) :: uco221(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8), intent(out) :: uco222(pcols,pverp)  ! co2 10.4 micron band path length
!
   real(r8), intent(out) :: uco223(pcols,pverp)  ! co2 10.4 micron band path length
   real(r8), intent(out) :: bn2o0(pcols,pverp)   ! pressure factor for n2o
   real(r8), intent(out) :: bn2o1(pcols,pverp)   ! pressure factor for n2o
   real(r8), intent(out) :: bch4(pcols,pverp)    ! pressure factor for ch4
   real(r8), intent(out) :: uptype(pcols,pverp)  ! p-type continuum path length

!
!---------------------------local variables-----------------------------
!
   integer   i               ! longitude index
   integer   k               ! level index
!
   real(r8) co2fac(pcols,1)      ! co2 factor
   real(r8) alpha1(pcols)        ! stimulated emission term
   real(r8) alpha2(pcols)        ! stimulated emission term
   real(r8) rt(pcols)            ! reciprocal of local temperature
   real(r8) rsqrt(pcols)         ! reciprocal of sqrt of temp
!
   real(r8) pbar(pcols)          ! mean pressure
   real(r8) dpnm(pcols)          ! difference in pressure
   real(r8) diff                 ! diffusivity factor
!
!--------------------------data statements------------------------------
!
   data diff /1.66/
!
!-----------------------------------------------------------------------
!
!  calculate path lengths for the trace gases at model top
!
   do i = 1,ncol
      ucfc11(i,ntoplw) = 1.8 * cfc11(i,ntoplw) * pnm(i,ntoplw) * rga
      ucfc12(i,ntoplw) = 1.8 * cfc12(i,ntoplw) * pnm(i,ntoplw) * rga
      un2o0(i,ntoplw) = diff * 1.02346e5 * n2o(i,ntoplw) * pnm(i,ntoplw) * rga / sqrt(tnm(i,ntoplw))
      un2o1(i,ntoplw) = diff * 2.01909 * un2o0(i,ntoplw) * exp(-847.36/tnm(i,ntoplw))
      uch4(i,ntoplw)  = diff * 8.60957e4 * ch4(i,ntoplw) * pnm(i,ntoplw) * rga / sqrt(tnm(i,ntoplw))
      co2fac(i,1)     = diff * co2mmr * pnm(i,ntoplw) * rga
      alpha1(i) = (1.0 - exp(-1540.0/tnm(i,ntoplw)))**3.0/sqrt(tnm(i,ntoplw))
      alpha2(i) = (1.0 - exp(-1360.0/tnm(i,ntoplw)))**3.0/sqrt(tnm(i,ntoplw))
      uco211(i,ntoplw) = 3.42217e3 * co2fac(i,1) * alpha1(i) * exp(-1849.7/tnm(i,ntoplw))
      uco212(i,ntoplw) = 6.02454e3 * co2fac(i,1) * alpha1(i) * exp(-2782.1/tnm(i,ntoplw))
      uco213(i,ntoplw) = 5.53143e3 * co2fac(i,1) * alpha1(i) * exp(-3723.2/tnm(i,ntoplw))
      uco221(i,ntoplw) = 3.88984e3 * co2fac(i,1) * alpha2(i) * exp(-1997.6/tnm(i,ntoplw))
      uco222(i,ntoplw) = 3.67108e3 * co2fac(i,1) * alpha2(i) * exp(-3843.8/tnm(i,ntoplw))
      uco223(i,ntoplw) = 6.50642e3 * co2fac(i,1) * alpha2(i) * exp(-2989.7/tnm(i,ntoplw))
      bn2o0(i,ntoplw) = diff * 19.399 * pnm(i,ntoplw)**2.0 * n2o(i,ntoplw) * &
                   1.02346e5 * rga / (sslp*tnm(i,ntoplw))
      bn2o1(i,ntoplw) = bn2o0(i,ntoplw) * exp(-847.36/tnm(i,ntoplw)) * 2.06646e5
      bch4(i,ntoplw) = diff * 2.94449 * ch4(i,ntoplw) * pnm(i,ntoplw)**2.0 * rga * &
                  8.60957e4 / (sslp*tnm(i,ntoplw))
      uptype(i,ntoplw) = diff * qnm(i,ntoplw) * pnm(i,ntoplw)**2.0 *  &
                    exp(1800.0*(1.0/tnm(i,ntoplw) - 1.0/296.0)) * rga / sslp
   end do
!
! calculate trace gas path lengths through model atmosphere
!
   do k = ntoplw,pver
      do i = 1,ncol
         rt(i) = 1./tnm(i,k)
         rsqrt(i) = sqrt(rt(i))
         pbar(i) = 0.5 * (pnm(i,k+1) + pnm(i,k)) / sslp
         dpnm(i) = (pnm(i,k+1) - pnm(i,k)) * rga
         alpha1(i) = diff * rsqrt(i) * (1.0 - exp(-1540.0/tnm(i,k)))**3.0
         alpha2(i) = diff * rsqrt(i) * (1.0 - exp(-1360.0/tnm(i,k)))**3.0
         ucfc11(i,k+1) = ucfc11(i,k) +  1.8 * cfc11(i,k) * dpnm(i)
         ucfc12(i,k+1) = ucfc12(i,k) +  1.8 * cfc12(i,k) * dpnm(i)
         un2o0(i,k+1) = un2o0(i,k) + diff * 1.02346e5 * n2o(i,k) * rsqrt(i) * dpnm(i)
         un2o1(i,k+1) = un2o1(i,k) + diff * 2.06646e5 * n2o(i,k) * &
                        rsqrt(i) * exp(-847.36/tnm(i,k)) * dpnm(i)
         uch4(i,k+1) = uch4(i,k) + diff * 8.60957e4 * ch4(i,k) * rsqrt(i) * dpnm(i)
         uco211(i,k+1) = uco211(i,k) + 1.15*3.42217e3 * alpha1(i) * &
                         co2mmr * exp(-1849.7/tnm(i,k)) * dpnm(i)
         uco212(i,k+1) = uco212(i,k) + 1.15*6.02454e3 * alpha1(i) * &
                         co2mmr * exp(-2782.1/tnm(i,k)) * dpnm(i)
         uco213(i,k+1) = uco213(i,k) + 1.15*5.53143e3 * alpha1(i) * &
                         co2mmr * exp(-3723.2/tnm(i,k)) * dpnm(i)
         uco221(i,k+1) = uco221(i,k) + 1.15*3.88984e3 * alpha2(i) * &
                         co2mmr * exp(-1997.6/tnm(i,k)) * dpnm(i)
         uco222(i,k+1) = uco222(i,k) + 1.15*3.67108e3 * alpha2(i) * &
                         co2mmr * exp(-3843.8/tnm(i,k)) * dpnm(i)
         uco223(i,k+1) = uco223(i,k) + 1.15*6.50642e3 * alpha2(i) * &
                         co2mmr * exp(-2989.7/tnm(i,k)) * dpnm(i)
         bn2o0(i,k+1) = bn2o0(i,k) + diff * 19.399 * pbar(i) * rt(i) &
                        * 1.02346e5 * n2o(i,k) * dpnm(i)
         bn2o1(i,k+1) = bn2o1(i,k) + diff * 19.399 * pbar(i) * rt(i) &
                        * 2.06646e5 * exp(-847.36/tnm(i,k)) * n2o(i,k)*dpnm(i)
         bch4(i,k+1) = bch4(i,k) + diff * 2.94449 * rt(i) * pbar(i) &
                       * 8.60957e4 * ch4(i,k) * dpnm(i)
         uptype(i,k+1) = uptype(i,k) + diff *qnm(i,k) * &
                         exp(1800.0*(1.0/tnm(i,k) - 1.0/296.0)) * pbar(i) * dpnm(i)
      end do
   end do
!
   return
end subroutine trcpth



subroutine aqsat(t       ,p       ,es      ,qs        ,ii      , &
                 ilen    ,kk      ,kstart  ,kend      )
!----------------------------------------------------------------------- 
! 
! purpose: 
! utility procedure to look up and return saturation vapor pressure from
! precomputed table, calculate and return saturation specific humidity
! (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
! this routine is useful for evaluating only a selected region in the
! vertical.
! 
! method: 
! <describe the algorithm(s) used in the routine.> 
! <also include any applicable external references.> 
! 
! author: j. hack
! 
!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: ii             ! i dimension of arrays t, p, es, qs
   integer, intent(in) :: kk             ! k dimension of arrays t, p, es, qs
   integer, intent(in) :: ilen           ! length of vectors in i direction which
   integer, intent(in) :: kstart         ! starting location in k direction
   integer, intent(in) :: kend           ! ending location in k direction
   real(r8), intent(in) :: t(ii,kk)          ! temperature
   real(r8), intent(in) :: p(ii,kk)          ! pressure
!
! output arguments
!
   real(r8), intent(out) :: es(ii,kk)         ! saturation vapor pressure
   real(r8), intent(out) :: qs(ii,kk)         ! saturation specific humidity
!
!---------------------------local workspace-----------------------------
!
   real(r8) omeps             ! 1 - 0.622
   integer i, k           ! indices
!
!-----------------------------------------------------------------------
!
   omeps = 1.0 - epsqs
   do k=kstart,kend
      do i=1,ilen
         es(i,k) = estblf(t(i,k))
!
! saturation specific humidity
!
         qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
!
! the following check is to avoid the generation of negative values
! that can occur in the upper stratosphere and mesosphere
!
         qs(i,k) = min(1.0_r8,qs(i,k))
!
         if (qs(i,k) < 0.0) then
            qs(i,k) = 1.0
            es(i,k) = p(i,k)
         end if
      end do
   end do
!
   return
end subroutine aqsat
!===============================================================================
  subroutine cldefr(lchnk   ,ncol    ,pcols, pver, pverp, &
       landfrac,t       ,rel     ,rei     ,ps      ,pmid    , landm, icefrac, snowh)
!----------------------------------------------------------------------- 
! 
! purpose: 
! compute cloud water and ice particle size 
! 
! method: 
! use empirical formulas to construct effective radii
! 
! author: j.t. kiehl, b. a. boville, p. rasch
! 
!-----------------------------------------------------------------------

    implicit none
!------------------------------arguments--------------------------------
!
! input arguments
!
    integer, intent(in) :: lchnk                 ! chunk identifier
    integer, intent(in) :: ncol                  ! number of atmospheric columns
    integer, intent(in) :: pcols, pver, pverp

    real(r8), intent(in) :: landfrac(pcols)      ! land fraction
    real(r8), intent(in) :: icefrac(pcols)       ! ice fraction
    real(r8), intent(in) :: t(pcols,pver)        ! temperature
    real(r8), intent(in) :: ps(pcols)            ! surface pressure
    real(r8), intent(in) :: pmid(pcols,pver)     ! midpoint pressures
    real(r8), intent(in) :: landm(pcols)
    real(r8), intent(in) :: snowh(pcols)         ! snow depth over land, water equivalent (m)
!
! output arguments
!
    real(r8), intent(out) :: rel(pcols,pver)      ! liquid effective drop size (microns)
    real(r8), intent(out) :: rei(pcols,pver)      ! ice effective drop size (microns)
!

!++pjr
! following kiehl
         call reltab(ncol, pcols, pver, t, landfrac, landm, icefrac, rel, snowh)

! following kristjansson and mitchell
         call reitab(ncol, pcols, pver, t, rei)
!--pjr
!
!
    return
  end subroutine cldefr


subroutine background(lchnk, ncol, pint, pcols, pverr, pverrp, mmr)
!-----------------------------------------------------------------------
!
! purpose:
! set global mean tropospheric aerosol background (or tuning) field
!
! method:
! specify aerosol mixing ratio.
! aerosol mass mixing ratio
! is specified so that the column visible aerosol optical depth is a
! specified global number (tauback). this means that the actual mixing
! ratio depends on pressure thickness of the lowest three atmospheric
! layers near the surface.
!
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use aer_optics, only: kbg,idxvis
!  use physconst, only: gravit
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!#include <ptrrgrid.h>
!------------------------------arguments--------------------------------
!
! input arguments
!
   integer, intent(in) :: lchnk                 ! chunk identifier
   integer, intent(in) :: ncol                  ! number of atmospheric columns
   integer, intent(in) :: pcols,pverr,pverrp

   real(r8), intent(in) :: pint(pcols,pverrp)   ! interface pressure (mks)
!
! output arguments
!
   real(r8), intent(out) :: mmr(pcols,pverr)    ! "background" aerosol mass mixing ratio
!
!---------------------------local variables-----------------------------
!
   integer i          ! longitude index
   integer k          ! level index
!
   real(r8) mass2mmr  ! factor to convert mass to mass mixing ratio
   real(r8) mass      ! mass of "background" aerosol as specified by tauback
!
!-----------------------------------------------------------------------
!
   do i=1,ncol
      mass2mmr =  gravmks / (pint(i,pverrp)-pint(i,pverrp-mxaerl))
      do k=1,pverr
!
! compute aerosol mass mixing ratio for specified levels (1.e3 factor is
! for units conversion of the extinction coefficiant from m2/g to m2/kg)
!
        if ( k >= pverrp-mxaerl ) then
! kaervs is not consistent with the values in aer_optics
! this ?should? be changed.
! rhfac is also implemented differently
            mass = tauback / (1.e3 * kbg(idxvis))
            mmr(i,k) = mass2mmr*mass
         else
            mmr(i,k) = 0._r8
         endif
!
      enddo
   enddo
!
   return
end subroutine background

subroutine scale_aerosols(aerosolt, pcols, pver, ncol, lchnk, scale)
!-----------------------------------------------------------------
! scale each species as determined by scale factors
!-----------------------------------------------------------------
  integer, intent(in) :: ncol, lchnk ! number of columns and chunk index
  integer, intent(in) :: pcols, pver
  real(r8), intent(in) :: scale(naer_all) ! scale each aerosol by this amount
  real(r8), intent(inout) :: aerosolt(pcols, pver, naer_all) ! aerosols
  integer m

  do m = 1, naer_all
     aerosolt(:ncol, :, m) = scale(m)*aerosolt(:ncol, :, m)
  end do

  return
end subroutine scale_aerosols

subroutine get_int_scales(scales)
  real(r8), intent(out)::scales(naer_all)  ! scale each aerosol by this amount
  integer i                                  ! index through species

!initialize
  scales = 1.

  scales(idxbg) = 1._r8
  scales(idxsul) = sulscl 
  scales(idxsslt) = ssltscl  
  
  do i = idxcarbonfirst, idxcarbonfirst+numcarbon-1
    scales(i) = carscl
  enddo
  
  do i = idxdustfirst, idxdustfirst+numdust-1
    scales(i) = dustscl
  enddo

  scales(idxvolc) = volcscl

  return
end subroutine get_int_scales

subroutine vert_interpolate (match_ps, aerosolc, m_hybi, paerlev, naer_c, pint, n, aerosol_mmr, pcols, pver, pverp, ncol, c)
!--------------------------------------------------------------------
! input: match surface pressure, cam interface pressure,
!        month index, number of columns, chunk index
!
! output: aerosol mass mixing ratio (aerosol_mmr)
!
! method:
!         interpolate column mass (cumulative) from match onto
!           cam's vertical grid (pressure coordinate)
!         convert back to mass mixing ratio
!
!--------------------------------------------------------------------

!  use physconst,     only: gravit

   integer, intent(in)  :: paerlev,naer_c,pcols,pver,pverp
   real(r8), intent(out) :: aerosol_mmr(pcols,pver,naer)  ! aerosol mmr from match
   real(r8), intent(in) :: match_ps(pcols)                ! surface pressure at a particular month
   real(r8), intent(in) :: pint(pcols,pverp)              ! interface pressure from cam
   real(r8), intent(in) :: aerosolc(pcols,paerlev,naer_c)
   real(r8), intent(in) :: m_hybi(paerlev)

   integer, intent(in) :: ncol,c                          ! chunk index and number of columns
   integer, intent(in) :: n                               ! prv or nxt month index
!
! local workspace
!
   integer m                           ! index to aerosol species
   integer kupper(pcols)               ! last upper bound for interpolation
   integer i, k, kk, kkstart, kount    ! loop vars for interpolation
   integer isv, ksv, msv               ! loop indices to save

   logical bad                         ! indicates a bad point found
   logical lev_interp_comp             ! interpolation completed for a level

   real(r8) aerosol(pcols,pverp,naer)  ! cumulative mass of aerosol in column beneath upper
                                       ! interface of level in column at particular month
   real(r8) dpl, dpu                   ! lower and upper intepolation factors
   real(r8) v_coord                    ! vertical coordinate
   real(r8) m_to_mmr                   ! mass to mass mixing ratio conversion factor
   real(r8) aer_diff                   ! temp var for difference between aerosol masses

!  call t_startf ('vert_interpolate')
!
! initialize index array
!
   do i=1,ncol
      kupper(i) = 1
   end do
!
! assign total mass to topmost level
!
   
   do i=1,ncol
   do m=1,naer
   aerosol(i,1,m) = aerosolc(i,1,m)
   enddo
   enddo
!
! at every pressure level, interpolate onto that pressure level
!
   do k=2,pver
!
! top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = paerlev
      do i=1,ncol
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! store level indices for interpolation
!
! for the pressure interpolation should be comparing
! pint(column,lev) with m_hybi(lev)*m_ps_cam_col(month,column,chunk)
!
      lev_interp_comp = .false.
      do kk=kkstart,paerlev-1
         if(.not.lev_interp_comp) then
         do i=1,ncol
            v_coord = pint(i,k)
            if (m_hybi(kk)*match_ps(i) .lt. v_coord .and. v_coord .le. m_hybi(kk+1)*match_ps(i)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! if all indices for this level have been found, do the interpolation and
! go to the next level
!
! interpolate in pressure.
!
         if (kount.eq.ncol) then
            do i=1,ncol
             do m=1,naer
               dpu = pint(i,k) - m_hybi(kupper(i))*match_ps(i)
               dpl = m_hybi(kupper(i)+1)*match_ps(i) - pint(i,k)
               aerosol(i,k,m) = &
                    (aerosolc(i,kupper(i)  ,m)*dpl + &
                     aerosolc(i,kupper(i)+1,m)*dpu)/(dpl + dpu)
             enddo
            enddo !i
            lev_interp_comp = .true.
         end if
         end if
      end do
!
! if we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and

! must extrapolate from the bottom or top pressure level for at least some
! of the longitude points.
!

      if(.not.lev_interp_comp) then
         do i=1,ncol
          do m=1,naer 
            if (pint(i,k) .lt. m_hybi(1)*match_ps(i)) then
               aerosol(i,k,m) =  aerosolc(i,1,m)
            else if (pint(i,k) .gt. m_hybi(paerlev)*match_ps(i)) then
               aerosol(i,k,m) = 0.0
            else
               dpu = pint(i,k) - m_hybi(kupper(i))*match_ps(i)
               dpl = m_hybi(kupper(i)+1)*match_ps(i) - pint(i,k)
               aerosol(i,k,m) = &
                    (aerosolc(i,kupper(i)  ,m)*dpl + &
                     aerosolc(i,kupper(i)+1,m)*dpu)/(dpl + dpu)
            end if
          enddo
         end do

         if (kount.gt.ncol) then
            call endrun ('vert_interpolate: bad data: non-monotonicity suspected in dependent variable')
         end if
      end if
   end do

!  call t_startf ('vi_checks')
!
! aerosol mass beneath lowest interface (pverp) must be 0
!
   aerosol(1:ncol,pverp,:) = 0.
!
! set mass in layer to zero whenever it is less than
!   1.e-40 kg/m^2 in the layer
!
   do m = 1, naer
      do k = 1, pver
         do i = 1, ncol
            if (aerosol(i,k,m) < 1.e-40_r8) aerosol(i,k,m) = 0.
         end do
      end do
   end do
!
! set mass in layer to zero whenever it is less than
!   10^-15 relative to column total mass
! convert back to mass mixing ratios.
! exit if mmr is negative
!
   do m = 1, naer
      do k = 1, pver
         do i = 1, ncol
            aer_diff = aerosol(i,k,m) - aerosol(i,k+1,m)
            if( abs(aer_diff) < 1e-15*aerosol(i,1,m)) then
               aer_diff = 0.
            end if
            m_to_mmr = gravmks / (pint(i,k+1)-pint(i,k))
            aerosol_mmr(i,k,m)= aer_diff * m_to_mmr
            if (aerosol_mmr(i,k,m) < 0) then
               write(6,*)'vert_interpolate: mmr < 0, m, col, lev, mmr',m, i, k, aerosol_mmr(i,k,m)
               write(6,*)'vert_interpolate: aerosol(k),(k+1)',aerosol(i,k,m),aerosol(i,k+1,m)
               write(6,*)'vert_interpolate: pint(k+1),(k)',pint(i,k+1),pint(i,k)
               write(6,*)'n,c',n,c
               call endrun()
            end if
         end do
      end do
   end do

!  call t_stopf ('vi_checks')
!  call t_stopf ('vert_interpolate')

   return
end subroutine vert_interpolate


!===============================================================================
  subroutine cldems(lchnk   ,ncol    ,pcols, pver, pverp, clwp    ,fice    ,rei     ,emis    )
!----------------------------------------------------------------------- 
! 
! purpose: 
! compute cloud emissivity using cloud liquid water path (g/m**2)
! 
! method: 
! <describe the algorithm(s) used in the routine.> 
! <also include any applicable external references.> 
! 
! author: j.t. kiehl
! 
!-----------------------------------------------------------------------

    implicit none
!------------------------------parameters-------------------------------
!
    real(r8) kabsl                  ! longwave liquid absorption coeff (m**2/g)
    parameter (kabsl = 0.090361)
!
!------------------------------arguments--------------------------------
!
! input arguments
!
    integer, intent(in) :: lchnk                   ! chunk identifier
    integer, intent(in) :: ncol                    ! number of atmospheric columns
    integer, intent(in) :: pcols, pver, pverp

    real(r8), intent(in) :: clwp(pcols,pver)       ! cloud liquid water path (g/m**2)
    real(r8), intent(in) :: rei(pcols,pver)        ! ice effective drop size (microns)
    real(r8), intent(in) :: fice(pcols,pver)       ! fractional ice content within cloud
!
! output arguments
!
    real(r8), intent(out) :: emis(pcols,pver)       ! cloud emissivity (fraction)
!
!---------------------------local workspace-----------------------------
!
    integer i,k                 ! longitude, level indices
    real(r8) kabs                   ! longwave absorption coeff (m**2/g)
    real(r8) kabsi                  ! ice absorption coefficient
!
!-----------------------------------------------------------------------
!
    do k=1,pver
       do i=1,ncol
          kabsi = 0.005 + 1./rei(i,k)
          kabs = kabsl*(1.-fice(i,k)) + kabsi*fice(i,k)
          emis(i,k) = 1. - exp(-1.66*kabs*clwp(i,k))
       end do
    end do
!
    return
  end subroutine cldems

!===============================================================================
  subroutine cldovrlap(lchnk   ,ncol    ,pcols, pver, pverp, pint    ,cld     ,nmxrgn  ,pmxrgn  )
!----------------------------------------------------------------------- 
! 
! purpose: 
! partitions each column into regions with clouds in neighboring layers.
! this information is used to implement maximum overlap in these regions
! with random overlap between them.
! on output,
!    nmxrgn contains the number of regions in each column
!    pmxrgn contains the interface pressures for the lower boundaries of
!           each region! 
! method: 

! 
! author: w. collins
! 
!-----------------------------------------------------------------------

    implicit none
!
! input arguments
!
    integer, intent(in) :: lchnk                ! chunk identifier
    integer, intent(in) :: ncol                 ! number of atmospheric columns
    integer, intent(in) :: pcols, pver, pverp

    real(r8), intent(in) :: pint(pcols,pverp)   ! interface pressure
    real(r8), intent(in) :: cld(pcols,pver)     ! fractional cloud cover
!
! output arguments
!
    real(r8), intent(out) :: pmxrgn(pcols,pverp)! maximum values of pressure for each
!    maximally overlapped region.
!    0->pmxrgn(i,1) is range of pressure for
!    1st region,pmxrgn(i,1)->pmxrgn(i,2) for
!    2nd region, etc
    integer nmxrgn(pcols)                    ! number of maximally overlapped regions
!
!---------------------------local variables-----------------------------
!
    integer i                    ! longitude index
    integer k                    ! level index
    integer n                    ! max-overlap region counter

    real(r8) pnm(pcols,pverp)    ! interface pressure

    logical cld_found            ! flag for detection of cloud
    logical cld_layer(pver)      ! flag for cloud in layer
!
!------------------------------------------------------------------------
!

    do i = 1, ncol
       cld_found = .false.
       cld_layer(:) = cld(i,:) > 0.0_r8
       pmxrgn(i,:) = 0.0
       pnm(i,:)=pint(i,:)*10.
       n = 1
       do k = 1, pver
          if (cld_layer(k) .and.  .not. cld_found) then
             cld_found = .true.
          else if ( .not. cld_layer(k) .and. cld_found) then
             cld_found = .false.
             if (count(cld_layer(k:pver)) == 0) then
                exit
             endif
             pmxrgn(i,n) = pnm(i,k)
             n = n + 1
          endif
       end do
       pmxrgn(i,n) = pnm(i,pverp)
       nmxrgn(i) = n
    end do

    return
  end subroutine cldovrlap

!===============================================================================
  subroutine cldclw(lchnk   ,ncol    ,pcols, pver, pverp, zi      ,clwp    ,tpw     ,hl      )
!----------------------------------------------------------------------- 
! 
! purpose: 
! evaluate cloud liquid water path clwp (g/m**2)
! 
! method: 
! <describe the algorithm(s) used in the routine.> 
! <also include any applicable external references.> 
! 
! author: j.t. kiehl
! 
!-----------------------------------------------------------------------

    implicit none

!
! input arguments
!
    integer, intent(in) :: lchnk                 ! chunk identifier
    integer, intent(in) :: ncol                  ! number of atmospheric columns
    integer, intent(in) :: pcols, pver, pverp

    real(r8), intent(in) :: zi(pcols,pverp)      ! height at layer interfaces(m)
    real(r8), intent(in) :: tpw(pcols)           ! total precipitable water (mm)
!
! output arguments
!
    real(r8) clwp(pcols,pver)     ! cloud liquid water path (g/m**2)
    real(r8) hl(pcols)            ! liquid water scale height
    real(r8) rhl(pcols)           ! 1/hl

!
!---------------------------local workspace-----------------------------
!
    integer i,k               ! longitude, level indices
    real(r8) clwc0                ! reference liquid water concentration (g/m**3)
    real(r8) emziohl(pcols,pverp) ! exp(-zi/hl)
!
!-----------------------------------------------------------------------
!
! set reference liquid water concentration
!
    clwc0 = 0.21
!
! diagnose liquid water scale height from precipitable water
!
    do i=1,ncol
       hl(i)  = 700.0*log(max(tpw(i)+1.0_r8,1.0_r8))
       rhl(i) = 1.0/hl(i)
    end do
!
! evaluate cloud liquid water path (vertical integral of exponential fn)
!
    do k=1,pverp
       do i=1,ncol
          emziohl(i,k) = exp(-zi(i,k)*rhl(i))
       end do
    end do
    do k=1,pver
       do i=1,ncol
          clwp(i,k) = clwc0*hl(i)*(emziohl(i,k+1) - emziohl(i,k))
       end do
    end do
!
    return
  end subroutine cldclw


!===============================================================================
  subroutine reltab(ncol, pcols, pver, t, landfrac, landm, icefrac, rel, snowh)
!----------------------------------------------------------------------- 
! 
! purpose: 
! compute cloud water size
! 
! method: 
! analytic formula following the formulation originally developed by j. t. kiehl
! 
! author: phil rasch
! 
!-----------------------------------------------------------------------
!   use physconst,          only: tmelt
    implicit none
!------------------------------arguments--------------------------------
!
! input arguments
!
    integer, intent(in) :: ncol
    integer, intent(in) :: pcols, pver
    real(r8), intent(in) :: landfrac(pcols)      ! land fraction
    real(r8), intent(in) :: icefrac(pcols)       ! ice fraction
    real(r8), intent(in) :: snowh(pcols)         ! snow depth over land, water equivalent (m)
    real(r8), intent(in) :: landm(pcols)         ! land fraction ramping to zero over ocean
    real(r8), intent(in) :: t(pcols,pver)        ! temperature

!
! output arguments
!
    real(r8), intent(out) :: rel(pcols,pver)      ! liquid effective drop size (microns)
!
!---------------------------local workspace-----------------------------
!
    integer i,k               ! lon, lev indices
    real(r8) rliqland         ! liquid drop size if over land
    real(r8) rliqocean        ! liquid drop size if over ocean
    real(r8) rliqice          ! liquid drop size if over sea ice
!
!-----------------------------------------------------------------------
!
    rliqocean = 14.0_r8
    rliqice   = 14.0_r8
    rliqland  = 8.0_r8
    do k=1,pver
       do i=1,ncol
! jrm reworked effective radius algorithm
          ! start with temperature-dependent value appropriate for continental air
          ! note: findmcnew has a pressure dependence here
          rel(i,k) = rliqland + (rliqocean-rliqland) * min(1.0_r8,max(0.0_r8,(tmelt-t(i,k))*0.05))
          ! modify for snow depth over land
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * min(1.0_r8,max(0.0_r8,snowh(i)*10.))
          ! ramp between polluted value over land to clean value over ocean.
          rel(i,k) = rel(i,k) + (rliqocean-rel(i,k)) * min(1.0_r8,max(0.0_r8,1.0-landm(i)))
          ! ramp between the resultant value and a sea ice value in the presence of ice.
          rel(i,k) = rel(i,k) + (rliqice-rel(i,k)) * min(1.0_r8,max(0.0_r8,icefrac(i)))
! end jrm
       end do
    end do
  end subroutine reltab
!===============================================================================
  subroutine reitab(ncol, pcols, pver, t, re)
    !

    integer, intent(in) :: ncol, pcols, pver
    real(r8), intent(out) :: re(pcols,pver)
    real(r8), intent(in) :: t(pcols,pver)
    real(r8) corr
    integer i
    integer k
    integer index
    !
    do k=1,pver
       do i=1,ncol
          index = int(t(i,k)-179.)
          index = min(max(index,1),94)
          corr = t(i,k) - int(t(i,k))
          re(i,k) = retab(index)*(1.-corr)		&
               +retab(index+1)*corr
          !           re(i,k) = amax1(amin1(re(i,k),30.),10.)
       end do
    end do
    !
    return
  end subroutine reitab
  
  function exp_interpol(x, f, y) result(g)

    ! purpose:
    !   interpolates f(x) to point y
    !   assuming f(x) = f(x0) exp a(x - x0)
    !   where a = ( ln f(x1) - ln f(x0) ) / (x1 - x0)
    !   x0 <= x <= x1
    !   assumes x is monotonically increasing

    ! author: d. fillmore

!   use shr_kind_mod, only: r8 => shr_kind_r8

    implicit none

    real(r8), intent(in), dimension(:) :: x  ! grid points
    real(r8), intent(in), dimension(:) :: f  ! grid function values
    real(r8), intent(in) :: y                ! interpolation point
    real(r8) :: g                            ! interpolated function value

    integer :: k  ! interpolation point index
    integer :: n  ! length of x
    real(r8) :: a

    n = size(x)

    ! find k such that x(k) < y =< x(k+1)
    ! set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)

    if (y <= x(1)) then
      k = 1
    else if (y >= x(n)) then
      k = n - 1
    else
      k = 1
      do while (y > x(k+1) .and. k < n)
        k = k + 1
      end do
    end if

    ! interpolate
    a = (  log( f(k+1) / f(k) )  ) / ( x(k+1) - x(k) )
    g = f(k) * exp( a * (y - x(k)) )

  end function exp_interpol

  function lin_interpol(x, f, y) result(g)
    
    ! purpose:
    !   interpolates f(x) to point y
    !   assuming f(x) = f(x0) + a * (x - x0)
    !   where a = ( f(x1) - f(x0) ) / (x1 - x0)
    !   x0 <= x <= x1
    !   assumes x is monotonically increasing

    ! author: d. fillmore

!   use shr_kind_mod, only: r8 => shr_kind_r8

    implicit none
    
    real(r8), intent(in), dimension(:) :: x  ! grid points
    real(r8), intent(in), dimension(:) :: f  ! grid function values
    real(r8), intent(in) :: y                ! interpolation point
    real(r8) :: g                            ! interpolated function value
    
    integer :: k  ! interpolation point index
    integer :: n  ! length of x
    real(r8) :: a

    n = size(x)

    ! find k such that x(k) < y =< x(k+1)
    ! set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)

    if (y <= x(1)) then 
      k = 1 
    else if (y >= x(n)) then
      k = n - 1
    else 
      k = 1 
      do while (y > x(k+1) .and. k < n)
        k = k + 1
      end do
    end if

    ! interpolate
    a = (  f(k+1) - f(k) ) / ( x(k+1) - x(k) )
    g = f(k) + a * (y - x(k))

  end function lin_interpol

  function lin_interpol2(x, f, y) result(g)

    ! purpose:
    !   interpolates f(x) to point y
    !   assuming f(x) = f(x0) + a * (x - x0)
    !   where a = ( f(x1) - f(x0) ) / (x1 - x0)
    !   x0 <= x <= x1
    !   assumes x is monotonically increasing

    ! author: d. fillmore ::  j. done changed from r8 to r4

    implicit none

    real, intent(in), dimension(:) :: x  ! grid points
    real, intent(in), dimension(:) :: f  ! grid function values
    real, intent(in) :: y                ! interpolation point
    real :: g                            ! interpolated function value

    integer :: k  ! interpolation point index
    integer :: n  ! length of x
    real    :: a

    n = size(x)

    ! find k such that x(k) < y =< x(k+1)
    ! set k = 1 if y <= x(1)  and  k = n-1 if y > x(n)

    if (y <= x(1)) then
      k = 1
    else if (y >= x(n)) then
      k = n - 1
    else
      k = 1
      do while (y > x(k+1) .and. k < n)
        k = k + 1
      end do
    end if

    ! interpolate
    a = (  f(k+1) - f(k) ) / ( x(k+1) - x(k) )
    g = f(k) + a * (y - x(k))

  end function lin_interpol2    


subroutine getfactors (cycflag, np1, cdayminus, cdayplus, cday, &
                       fact1, fact2)
!---------------------------------------------------------------------------
!
! purpose: determine time interpolation factors (normally for a boundary dataset)
!          for linear interpolation.
!
! method:  assume 365 days per year.  output variable fact1 will be the weight to
!          apply to data at calendar time "cdayminus", and fact2 the weight to apply
!          to data at time "cdayplus".  combining these values will produce a result
!          valid at time "cday".  output arguments fact1 and fact2 will be between
!          0 and 1, and fact1 + fact2 = 1 to roundoff.
!
! author:  jim rosinski
!
!---------------------------------------------------------------------------
   implicit none
!
! arguments
!
   logical, intent(in) :: cycflag             ! flag indicates whether dataset is being cycled yearly

   integer, intent(in) :: np1                 ! index points to forward time slice matching cdayplus

   real(r8), intent(in) :: cdayminus          ! calendar day of rearward time slice
   real(r8), intent(in) :: cdayplus           ! calendar day of forward time slice
   real(r8), intent(in) :: cday               ! calenar day to be interpolated to
   real(r8), intent(out) :: fact1             ! time interpolation factor to apply to rearward time slice
   real(r8), intent(out) :: fact2             ! time interpolation factor to apply to forward time slice

!  character(len=*), intent(in) :: str        ! string to be added to print in case of error (normally the callers name)
!
! local workspace
!
   real(r8) :: deltat                         ! time difference (days) between cdayminus and cdayplus
   real(r8), parameter :: daysperyear = 365.  ! number of days in a year
!
! initial sanity checks
!
!  if (np1 == 1 .and. .not. cycflag) then
!     call endrun ('getfactors:'//str//' cycflag false and forward month index = jan. not allowed')
!  end if

!  if (np1 < 1) then
!     call endrun ('getfactors:'//str//' input arg np1 must be > 0')
!  end if

   if (cycflag) then
      if ((cday < 1.) .or. (cday > (daysperyear+1.))) then
         write(6,*) 'getfactors:', ' bad cday=',cday
         call endrun ()
      end if
   else
      if (cday < 1.) then
         write(6,*) 'getfactors:',  ' bad cday=',cday
         call endrun ()
      end if
   end if
!
! determine time interpolation factors.  account for december-january
! interpolation if dataset is being cycled yearly.
!
   if (cycflag .and. np1 == 1) then                     ! dec-jan interpolation
      deltat = cdayplus + daysperyear - cdayminus
      if (cday > cdayplus) then                         ! we are in december
         fact1 = (cdayplus + daysperyear - cday)/deltat
         fact2 = (cday - cdayminus)/deltat
      else                                              ! we are in january
         fact1 = (cdayplus - cday)/deltat
         fact2 = (cday + daysperyear - cdayminus)/deltat
      end if
   else
      deltat = cdayplus - cdayminus
      fact1 = (cdayplus - cday)/deltat
      fact2 = (cday - cdayminus)/deltat
   end if

   if (.not. validfactors (fact1, fact2)) then
      write(6,*) 'getfactors: ', ' bad fact1 and/or fact2=', fact1, fact2
      call endrun ()
   end if

   return
end subroutine getfactors

logical function validfactors (fact1, fact2)
!---------------------------------------------------------------------------
!
! purpose: check sanity of time interpolation factors to within 32-bit roundoff
!
!---------------------------------------------------------------------------
   implicit none

   real(r8), intent(in) :: fact1, fact2           ! time interpolation factors

   validfactors = .true.
   if (abs(fact1+fact2-1.) > 1.e-6 .or. &
       fact1 > 1.000001 .or. fact1 < -1.e-6 .or. &
       fact2 > 1.000001 .or. fact2 < -1.e-6) then

      validfactors = .false.
   end if

   return
end function validfactors

subroutine get_rf_scales(scales)

  real(r8), intent(out)::scales(naer_all)  ! scale aerosols by this amount

  integer i                                  ! loop index

  scales(idxbg) = bgscl_rf
  scales(idxsul) = sulscl_rf
  scales(idxsslt) = ssltscl_rf

  do i = idxcarbonfirst, idxcarbonfirst+numcarbon-1
    scales(i) = carscl_rf
  enddo

  do i = idxdustfirst, idxdustfirst+numdust-1
    scales(i) = dustscl_rf
  enddo

  scales(idxvolc) = volcscl_rf

end subroutine get_rf_scales

function psi(tpx,iband)
!    
! history: first version for hitran 1996 (c/h/e)
!          current version for hitran 2000 (c/lt/e)
! short function for hulst-curtis-godson temperature factors for
!   computing effective h2o path
! line data for h2o: hitran 2000, plus h2o patches v11.0 for 1341 missing
!                    lines between 500 and 2820 cm^-1.
!                    see cfa-www.harvard.edu/hitran
! isotopes of h2o: all
! line widths: air-broadened only (self set to 0)
! code for line strengths and widths: genln3
! reference: edwards, d.p., 1992: genln2, a general line-by-line atmospheric
!                     transmittance and radiance model, version 3.0 description
!                     and users guide, ncar/tn-367+str, 147 pp.
!     
! note: functions have been normalized by dividing by their values at
!       a path temperature of 160k
!
! spectral intervals:     
!   1 = 0-800 cm^-1 and 1200-2200 cm^-1
!   2 = 800-1200 cm^-1      
!
! formulae: goody and yung, atmospheric radiation: theoretical basis,
!           2nd edition, oxford university press, 1989.
! psi: function for pressure along path
!      eq. 6.30, p. 228
!
   real(r8),intent(in):: tpx      ! path temperature
   integer, intent(in):: iband    ! band to process
   real(r8) psi                   ! psi for given band
   real(r8),parameter ::  psi_r0(nbands) = (/ 5.65308452e-01, -7.30087891e+01/)
   real(r8),parameter ::  psi_r1(nbands) = (/ 4.07519005e-03,  1.22199547e+00/)
   real(r8),parameter ::  psi_r2(nbands) = (/-1.04347237e-05, -7.12256227e-03/)
   real(r8),parameter ::  psi_r3(nbands) = (/ 1.23765354e-08,  1.47852825e-05/)

   psi = (((psi_r3(iband) * tpx) + psi_r2(iband)) * tpx + psi_r1(iband)) * tpx + psi_r0(iband)
end function psi

function phi(tpx,iband)
!
! history: first version for hitran 1996 (c/h/e)
!          current version for hitran 2000 (c/lt/e)
! short function for hulst-curtis-godson temperature factors for
!   computing effective h2o path
! line data for h2o: hitran 2000, plus h2o patches v11.0 for 1341 missing
!                    lines between 500 and 2820 cm^-1.
!                    see cfa-www.harvard.edu/hitran
! isotopes of h2o: all
! line widths: air-broadened only (self set to 0)
! code for line strengths and widths: genln3
! reference: edwards, d.p., 1992: genln2, a general line-by-line atmospheric
!                     transmittance and radiance model, version 3.0 description
!                     and users guide, ncar/tn-367+str, 147 pp.
!
! note: functions have been normalized by dividing by their values at
!       a path temperature of 160k
!
! spectral intervals:
!   1 = 0-800 cm^-1 and 1200-2200 cm^-1
!   2 = 800-1200 cm^-1
!
! formulae: goody and yung, atmospheric radiation: theoretical basis,
!           2nd edition, oxford university press, 1989.
! phi: function for h2o path
!      eq. 6.25, p. 228
!
   real(r8),intent(in):: tpx      ! path temperature
   integer, intent(in):: iband    ! band to process
   real(r8) phi                   ! phi for given band
   real(r8),parameter ::  phi_r0(nbands) = (/ 9.60917711e-01, -2.21031342e+01/)
   real(r8),parameter ::  phi_r1(nbands) = (/ 4.86076751e-04,  4.24062610e-01/)
   real(r8),parameter ::  phi_r2(nbands) = (/-1.84806265e-06, -2.95543415e-03/)
   real(r8),parameter ::  phi_r3(nbands) = (/ 2.11239959e-09,  7.52470896e-06/)

   phi = (((phi_r3(iband) * tpx) + phi_r2(iband)) * tpx + phi_r1(iband)) &
          * tpx + phi_r0(iband)
end function phi

function fh2oself( temp )
!
! short function for h2o self-continuum temperature factor in
!   calculation of effective h2o self-continuum path length
!
! h2o continuum: ckd 2.4
! code for continuum: genln3
! reference: edwards, d.p., 1992: genln2, a general line-by-line atmospheric
!                     transmittance and radiance model, version 3.0 description
!                     and users guide, ncar/tn-367+str, 147 pp.
!
! in genln, the temperature scaling of the self-continuum is handled
!    by exponential interpolation/extrapolation from observations at
!    260k and 296k by:
!
!         tfac =  (t(ipath) - 296.0)/(260.0 - 296.0)
!         csfft = csff296*(csff260/csff296)**tfac
!
! for 800-1200 cm^-1, (csff260/csff296) ranges from ~2.1 to ~1.9
!     with increasing wavenumber.  the ratio <csff260>/<csff296>,
!     where <> indicates average over wavenumber, is ~2.07
!
! fh2oself is (<csff260>/<csff296>)**tfac
!
   real(r8),intent(in) :: temp     ! path temperature
   real(r8) fh2oself               ! mean ratio of self-continuum at temp and 296k

   fh2oself = 2.0727484**((296.0 - temp) / 36.0)
end function fh2oself

! from wv_saturation.f90

subroutine esinti(epslon  ,latvap  ,latice  ,rh2o    ,cpair   ,tmelt   )
!----------------------------------------------------------------------- 
! 
! purpose: 
! initialize es lookup tables
! 
! method: 
! <describe the algorithm(s) used in the routine.> 
! <also include any applicable external references.> 
! 
! author: j. hack
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use wv_saturation, only: gestbl
   implicit none
!------------------------------arguments--------------------------------
!
! input arguments
!
   real(r8), intent(in) :: epslon          ! ratio of h2o to dry air molecular weights
   real(r8), intent(in) :: latvap          ! latent heat of vaporization
   real(r8), intent(in) :: latice          ! latent heat of fusion
   real(r8), intent(in) :: rh2o            ! gas constant for water vapor
   real(r8), intent(in) :: cpair           ! specific heat of dry air
   real(r8), intent(in) :: tmelt           ! melting point of water (k)
!
!---------------------------local workspace-----------------------------
!
   real(r8) tmn             ! minimum temperature entry in table
   real(r8) tmx             ! maximum temperature entry in table
   real(r8) trice           ! trans range from es over h2o to es over ice
   logical ip           ! ice phase (true or false)
!
!-----------------------------------------------------------------------
!
! specify control parameters first
!
   tmn   = 173.16
   tmx   = 375.16
   trice =  20.00
   ip    = .true.
!
! call gestbl to build saturation vapor pressure table.
!
   call gestbl(tmn     ,tmx     ,trice   ,ip      ,epslon  , &
               latvap  ,latice  ,rh2o    ,cpair   ,tmelt )
!
   return
end subroutine esinti

subroutine gestbl(tmn     ,tmx     ,trice   ,ip      ,epsil   , &
                  latvap  ,latice  ,rh2o    ,cpair   ,tmeltx   )
!-----------------------------------------------------------------------
!
! purpose:
! builds saturation vapor pressure table for later lookup procedure.
!
! method:
! uses goff & gratch (1946) relationships to generate the table
! according to a set of free parameters defined below.  auxiliary
! routines are also included for making rapid estimates (well with 1%)
! of both es and d(es)/dt for the particular table configuration.
!
! author: j. hack
!
!-----------------------------------------------------------------------
!  use pmgrid, only: masterproc
   implicit none
!------------------------------arguments--------------------------------
!
! input arguments
!
   real(r8), intent(in) :: tmn           ! minimum temperature entry in es lookup table
   real(r8), intent(in) :: tmx           ! maximum temperature entry in es lookup table
   real(r8), intent(in) :: epsil         ! ratio of h2o to dry air molecular weights
   real(r8), intent(in) :: trice         ! transition range from es over range to es over ice
   real(r8), intent(in) :: latvap        ! latent heat of vaporization
   real(r8), intent(in) :: latice        ! latent heat of fusion
   real(r8), intent(in) :: rh2o          ! gas constant for water vapor
   real(r8), intent(in) :: cpair         ! specific heat of dry air
   real(r8), intent(in) :: tmeltx        ! melting point of water (k)
!
!---------------------------local variables-----------------------------
!
   real(r8) t             ! temperature
   real(r8) rgasv 
   real(r8) cp
   real(r8) hlatf
   real(r8) ttrice
   real(r8) hlatv
   integer n          ! increment counter
   integer lentbl     ! calculated length of lookup table
   integer itype      ! ice phase: 0 -> no ice phase
!            1 -> ice phase, no transition
!           -x -> ice phase, x degree transition
   logical ip         ! ice phase logical flag
   logical icephs
!
!-----------------------------------------------------------------------
!
! set es table parameters
!
   tmin   = tmn       ! minimum temperature entry in table
   tmax   = tmx       ! maximum temperature entry in table
   ttrice = trice     ! trans. range from es over h2o to es over ice
   icephs = ip        ! ice phase (true or false)
!
! set physical constants required for es calculation
!
   epsqs  = epsil
   hlatv  = latvap
   hlatf  = latice
   rgasv  = rh2o
   cp     = cpair
   tmelt  = tmeltx
!
   lentbl = int(tmax-tmin+2.000001)
   if (lentbl .gt. plenest) then
      write(6,9000) tmax, tmin, plenest
      call endrun ('gestbl')    ! abnormal termination
   end if
!
! begin building es table.
! check whether ice phase requested.
! if so, set appropriate transition range for temperature
!
   if (icephs) then
      if (ttrice /= 0.0) then
         itype = -ttrice
      else
         itype = 1
      end if
   else
      itype = 0
   end if
!
   t = tmin - 1.0
   do n=1,lentbl
      t = t + 1.0
      call gffgch(t,estbl(n),itype)
   end do
!
   do n=lentbl+1,plenest
      estbl(n) = -99999.0
   end do
!
! table complete -- set coefficients for polynomial approximation of
! difference between saturation vapor press over water and saturation
! pressure over ice for -ttrice < t < 0 (degrees c). note: polynomial
! is valid in the range -40 < t < 0 (degrees c).
!
!                  --- degree 5 approximation ---
!
   pcf(1) =  5.04469588506e-01
   pcf(2) = -5.47288442819e+00
   pcf(3) = -3.67471858735e-01
   pcf(4) = -8.95963532403e-03
   pcf(5) = -7.78053686625e-05
!
!                  --- degree 6 approximation ---
!
!-----pcf(1) =  7.63285250063e-02
!-----pcf(2) = -5.86048427932e+00
!-----pcf(3) = -4.38660831780e-01
!-----pcf(4) = -1.37898276415e-02
!-----pcf(5) = -2.14444472424e-04
!-----pcf(6) = -1.36639103771e-06
!
   if (masterproc) then
      write(6,*)' *** saturation vapor pressure table completed ***'
   end if

   return
!
9000 format('gestbl: fatal error *********************************',/, &
            ' tmax and tmin require a larger dimension on the length', &
            ' of the saturation vapor pressure table estbl(plenest)',/, &
            ' tmax, tmin, and plenest => ', 2f7.2, i3)
!
end subroutine gestbl

subroutine gffgch(t       ,es      ,itype   )
!----------------------------------------------------------------------- 
! 
! purpose: 
! computes saturation vapor pressure over water and/or over ice using
! goff & gratch (1946) relationships. 
! <say what the routine does> 
! 
! method: 
! t (temperature), and itype are input parameters, while es (saturation
! vapor pressure) is an output parameter.  the input parameter itype
! serves two purposes: a value of zero indicates that saturation vapor
! pressures over water are to be returned (regardless of temperature),
! while a value of one indicates that saturation vapor pressures over
! ice should be returned when t is less than freezing degrees.  if itype
! is negative, its absolute value is interpreted to define a temperature
! transition region below freezing in which the returned
! saturation vapor pressure is a weighted average of the respective ice
! and water value.  that is, in the temperature range 0 => -itype
! degrees c, the saturation vapor pressures are assumed to be a weighted
! average of the vapor pressure over supercooled water and ice (all
! water at 0 c; all ice at -itype c).  maximum transition range => 40 c
! 
! author: j. hack
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use physconst, only: tmelt
!  use abortutils, only: endrun
    
   implicit none
!------------------------------arguments--------------------------------
!
! input arguments
!
   real(r8), intent(in) :: t          ! temperature
!
! output arguments
!
   integer, intent(inout) :: itype   ! flag for ice phase and associated transition

   real(r8), intent(out) :: es         ! saturation vapor pressure
!
!---------------------------local variables-----------------------------
!
   real(r8) e1         ! intermediate scratch variable for es over water
   real(r8) e2         ! intermediate scratch variable for es over water
   real(r8) eswtr      ! saturation vapor pressure over water
   real(r8) f          ! intermediate scratch variable for es over water
   real(r8) f1         ! intermediate scratch variable for es over water
   real(r8) f2         ! intermediate scratch variable for es over water
   real(r8) f3         ! intermediate scratch variable for es over water
   real(r8) f4         ! intermediate scratch variable for es over water
   real(r8) f5         ! intermediate scratch variable for es over water
   real(r8) ps         ! reference pressure (mb)
   real(r8) t0         ! reference temperature (freezing point of water)
   real(r8) term1      ! intermediate scratch variable for es over ice
   real(r8) term2      ! intermediate scratch variable for es over ice
   real(r8) term3      ! intermediate scratch variable for es over ice
   real(r8) tr         ! transition range for es over water to es over ice
   real(r8) ts         ! reference temperature (boiling point of water)
   real(r8) weight     ! intermediate scratch variable for es transition
   integer itypo   ! intermediate scratch variable for holding itype
!
!-----------------------------------------------------------------------
!
! check on whether there is to be a transition region for es
!
   if (itype < 0) then
      tr    = abs(float(itype))
      itypo = itype
      itype = 1
   else
      tr    = 0.0
      itypo = itype
   end if
   if (tr > 40.0) then
      write(6,900) tr
      call endrun ('gffgch')                ! abnormal termination
   end if
!
   if(t < (tmelt - tr) .and. itype == 1) go to 10
!
! water
!
   ps = 1013.246
   ts = 373.16
   e1 = 11.344*(1.0 - t/ts)
   e2 = -3.49149*(ts/t - 1.0)
   f1 = -7.90298*(ts/t - 1.0)
   f2 = 5.02808*log10(ts/t)
   f3 = -1.3816*(10.0**e1 - 1.0)/10000000.0
   f4 = 8.1328*(10.0**e2 - 1.0)/1000.0
   f5 = log10(ps)
   f  = f1 + f2 + f3 + f4 + f5
   es = (10.0**f)*100.0
   eswtr = es
!
   if(t >= tmelt .or. itype == 0) go to 20
!
! ice
!
10 continue
   t0    = tmelt
   term1 = 2.01889049/(t0/t)
   term2 = 3.56654*log(t0/t)
   term3 = 20.947031*(t0/t)
   es    = 575.185606e10*exp(-(term1 + term2 + term3))
!
   if (t < (tmelt - tr)) go to 20
!
! weighted transition between water and ice
!
   weight = min((tmelt - t)/tr,1.0_r8)
   es = weight*es + (1.0 - weight)*eswtr
!
20 continue
   itype = itypo
   return
!
900 format('gffgch: fatal error ******************************',/, &
           'transition range for water to ice saturation vapor', &
           ' pressure, tr, exceeds maximum allowable value of', &
           ' 40.0 degrees c',/, ' tr = ',f7.2)
!
end subroutine gffgch

   real(r8) function estblf( td )
!
! saturation vapor pressure table lookup
!
   real(r8), intent(in) :: td         ! temperature for saturation lookup
!
   real(r8) :: e       ! intermediate variable for es look-up
   real(r8) :: ai
   integer  :: i
!
   e = max(min(td,tmax),tmin)   ! partial pressure
   i = int(e-tmin)+1
   ai = aint(e-tmin)
   estblf = (tmin+ai-e+1.)* &
            estbl(i)-(tmin+ai-e)* &
            estbl(i+1)
   end function estblf


function findvalue(ix,n,ain,indxa)
!----------------------------------------------------------------------- 
! 
! purpose: 
! subroutine for finding ix-th smallest value in the array
! the elements are rearranged so that the ix-th smallest
! element is in the ix place and all smaller elements are
! moved to the elements up to ix (with random order).
!
! algorithm: based on the quicksort algorithm.
!
! author:       t. craig
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
   implicit none
!
! arguments
!
   integer, intent(in) :: ix                ! element to search for
   integer, intent(in) :: n                 ! total number of elements
   integer, intent(inout):: indxa(n)        ! array of integers
   real(r8), intent(in) :: ain(n)           ! array to search
!
   integer findvalue                        ! return value
!
! local variables
!
   integer i,j
   integer il,im,ir

   integer ia
   integer itmp
!
!---------------------------routine-----------------------------
!
   il=1
   ir=n
   do
      if (ir-il <= 1) then
         if (ir-il == 1) then
            if (ain(indxa(ir)) < ain(indxa(il))) then
               itmp=indxa(il)
               indxa(il)=indxa(ir)
               indxa(ir)=itmp
            endif
         endif
         findvalue=indxa(ix)
         return
      else
         im=(il+ir)/2
         itmp=indxa(im)
         indxa(im)=indxa(il+1)
         indxa(il+1)=itmp
         if (ain(indxa(il+1)) > ain(indxa(ir))) then
            itmp=indxa(il+1)
            indxa(il+1)=indxa(ir)
            indxa(ir)=itmp
         endif
         if (ain(indxa(il)) > ain(indxa(ir))) then
            itmp=indxa(il)
            indxa(il)=indxa(ir)
            indxa(ir)=itmp
         endif
         if (ain(indxa(il+1)) > ain(indxa(il))) then
            itmp=indxa(il+1)
            indxa(il+1)=indxa(il)
            indxa(il)=itmp
         endif
         i=il+1
         j=ir
         ia=indxa(il)
         do
            do
               i=i+1
               if (ain(indxa(i)) >= ain(ia)) exit
            end do
            do
               j=j-1
               if (ain(indxa(j)) <= ain(ia)) exit
            end do
            if (j < i) exit
            itmp=indxa(i)
            indxa(i)=indxa(j)
            indxa(j)=itmp
         end do
         indxa(il)=indxa(j)
         indxa(j)=ia
         if (j >= ix)ir=j-1
         if (j <= ix)il=i
      endif
   end do
end function findvalue


subroutine radini(gravx   ,cpairx  ,epsilox ,stebolx, pstdx )
!----------------------------------------------------------------------- 
! 
! purpose: 
! initialize various constants for radiation scheme; note that
! the radiation scheme uses cgs units.
! 
! method: 
! <describe the algorithm(s) used in the routine.> 
! <also include any applicable external references.> 
! 
! author: w. collins (h2o parameterization) and j. kiehl
! 
!-----------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid,       only: pver, pverp
!  use comozp,       only: cplos, cplol
!  use pmgrid,       only: masterproc, plev, plevp
!  use radae,        only: radaeini
!  use physconst,    only: mwdry, mwco2
#if ( defined spmd )
!   use mpishorthand
#endif
   implicit none

!------------------------------arguments--------------------------------
!
! input arguments
!
   real, intent(in) :: gravx      ! acceleration of gravity (mks)
   real, intent(in) :: cpairx     ! specific heat of dry air (mks)
   real, intent(in) :: epsilox    ! ratio of mol. wght of h2o to dry air
   real, intent(in) :: stebolx    ! stefan-boltzmann's constant (mks)
   real(r8), intent(in) :: pstdx      ! standard pressure (pascals)
!
!---------------------------local variables-----------------------------
!
   integer k       ! loop variable

   real(r8) v0         ! volume of a gas at stp (m**3/kmol)
   real(r8) p0         ! standard pressure (pascals)
   real(r8) amd        ! effective molecular weight of dry air (kg/kmol)
   real(r8) goz        ! acceleration of gravity (m/s**2)
!
!-----------------------------------------------------------------------
!
! set general radiation consts; convert to cgs units where appropriate:
!
   gravit  =  100.*gravx
   rga     =  1./gravit
   gravmks =  gravx
   cpair   =  1.e4*cpairx
   epsilo  =  epsilox
   sslp    =  1.013250e6
   stebol  =  1.e3*stebolx
   rgsslp  =  0.5/(gravit*sslp)
   dpfo3   =  2.5e-3
   dpfco2  =  5.0e-3
   dayspy  =  365.
   pie     =  4.*atan(1.)
!
! initialize ozone data.
!
   v0  = 22.4136         ! volume of a gas at stp (m**3/kmol)
   p0  = 0.1*sslp        ! standard pressure (pascals)
   amd = 28.9644         ! molecular weight of dry air (kg/kmol)
   goz = gravx           ! acceleration of gravity (m/s**2)
!
! constants for ozone path integrals (multiplication by 100 for unit
! conversion to cgs from mks):
!
   cplos = v0/(amd*goz)       *100.0
   cplol = v0/(amd*goz*p0)*0.5*100.0
!
! derived constants
! if the top model level is above ~90 km (0.1 pa), set the top level to compute
! longwave cooling to about 80 km (1 pa)
! wrf: assume top level > 0.1 mb
!  if (hypm(1) .lt. 0.1) then
!     do k = 1, pver
!        if (hypm(k) .lt. 1.) ntoplw  = k
!     end do
!  else
      ntoplw = 1
!  end if
!   if (masterproc) then
!     write (6,*) 'radini: ntoplw =',ntoplw, ' pressure:',hypm(ntoplw)
!   endif

   call radaeini( pstdx, mwdry, mwco2 )
   return
end subroutine radini

subroutine oznini(ozmixm,pin,levsiz,num_months,xlat,                &
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte)
!
! this subroutine assumes uniform distribution of ozone concentration.
! it should be replaced by monthly climatology that varies latitudinally and vertically
!

      implicit none

   integer,      intent(in   )    ::   ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte   

   integer,      intent(in   )    ::   levsiz, num_months

   real,  dimension( ims:ime, jms:jme ), intent(in   )  ::     xlat

   real,  dimension( ims:ime, levsiz, jms:jme, num_months ),      &
          intent(out   ) ::                                  ozmixm

   real,  dimension(levsiz), intent(out )  ::                   pin

! local
   integer, parameter :: latsiz = 64
   integer, parameter :: lonsiz = 1
   integer :: i, j, k, itf, jtf, ktf, m, pin_unit, lat_unit, oz_unit
   real    :: interp_pt
   character*256 :: message

   real,  dimension( lonsiz, levsiz, latsiz, num_months )    ::   &
                                                            ozmixin

   real,  dimension(latsiz)                ::             lat_ozone
#ifdef AMIPW_PHYSICS
   jtf=jte
   ktf=kte
   itf=ite
#else
   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)
#endif


!-- read in ozone pressure data

     write(message,*)'num_months = ',num_months
#ifndef AMIPW_PHYSICS
     call wrf_debug(50,message)
#endif

      pin_unit = 27
        open(pin_unit, file='ozone_plev.formatted',form='formatted',status='old')
        do k = 1,levsiz
        read (pin_unit,*)pin(k)
        end do
      close(27)

      do k=1,levsiz
        pin(k) = pin(k)*100.
      end do

!-- read in ozone lat data

      lat_unit = 28
        open(lat_unit, file='ozone_lat.formatted',form='formatted',status='old')
        do j = 1,latsiz
        read (lat_unit,*)lat_ozone(j)
        end do
      close(28)


!-- read in ozone data

      oz_unit = 29
      open(oz_unit, file='ozone.formatted',form='formatted',status='old')

      do m=2,num_months
      do j=1,latsiz ! latsiz=64
      do k=1,levsiz ! levsiz=59
      do i=1,lonsiz ! lonsiz=1
        read (oz_unit,*)ozmixin(i,k,j,m)
      enddo
      enddo
      enddo
      enddo
      close(29)


!-- latitudinally interpolate ozone data (and extend longitudinally)
!-- using function lin_interpol2(x, f, y) result(g)
! purpose:
!   interpolates f(x) to point y
!   assuming f(x) = f(x0) + a * (x - x0)
!   where a = ( f(x1) - f(x0) ) / (x1 - x0)
!   x0 <= x <= x1
!   assumes x is monotonically increasing
!    real, intent(in), dimension(:) :: x  ! grid points
!    real, intent(in), dimension(:) :: f  ! grid function values
!    real, intent(in) :: y                ! interpolation point
!    real :: g                            ! interpolated function value
!---------------------------------------------------------------------------

      do m=2,num_months
      do j=jts,jtf
      do k=1,levsiz
      do i=its,itf
         interp_pt=xlat(i,j)
         ozmixm(i,k,j,m)=lin_interpol2(lat_ozone(:),ozmixin(1,k,:,m),interp_pt)
      enddo
      enddo
      enddo
      enddo

! old code for fixed ozone

!     pin(1)=70.
!     do k=2,levsiz
!     pin(k)=pin(k-1)+16.
!     enddo

!     do k=1,levsiz
!         pin(k) = pin(k)*100.
!     end do

!     do m=1,num_months
!     do j=jts,jtf
!     do i=its,itf
!     do k=1,2
!      ozmixm(i,k,j,m)=1.e-6
!     enddo
!     do k=3,levsiz
!      ozmixm(i,k,j,m)=1.e-7
!     enddo
!     enddo
!     enddo
!     enddo

end subroutine oznini


subroutine aerosol_init(m_psp,m_psn,m_hybi,aerosolcp,aerosolcn,paerlev,naer_c,shalf,pptop,    &
#ifdef AMIPW_PHYSICS
                     hypm, &
#endif
                     ids, ide, jds, jde, kds, kde,                  &
                     ims, ime, jms, jme, kms, kme,                  &
                     its, ite, jts, jte, kts, kte)
!
!  this subroutine assumes a uniform aerosol distribution in both time and space.
!  it should be modified if aerosol data are available from wrf-chem or other sources
!
      implicit none

   integer,      intent(in   )    ::   ids,ide, jds,jde, kds,kde, &
                                       ims,ime, jms,jme, kms,kme, &
                                       its,ite, jts,jte, kts,kte

   integer,      intent(in   )    ::   paerlev,naer_c 

   real,     intent(in)                        :: pptop
   real,     dimension( kms:kme ), intent(in)  :: shalf
#ifdef AMIPW_PHYSICS
   real,     dimension( kms:kme ), intent(in)  :: hypm
#endif
   real,  dimension( ims:ime, paerlev, jms:jme, naer_c ),      &
          intent(inout   ) ::                                  aerosolcn , aerosolcp

   real,  dimension(paerlev), intent(out )  ::                m_hybi
   real,  dimension( ims:ime, jms:jme),  intent(out )  ::       m_psp,m_psn 

   real ::                                                      psurf
   real, dimension(29) :: hybi  
   integer k ! index through vertical levels

   integer :: i, j, itf, jtf, ktf,m

   data hybi/0, 0.0065700002014637, 0.0138600002974272, 0.023089999333024, &
    0.0346900001168251, 0.0491999983787537, 0.0672300010919571,      &
     0.0894500017166138, 0.116539999842644, 0.149159997701645,       &
    0.187830001115799, 0.232859998941422, 0.284209996461868,         &
    0.341369986534119, 0.403340011835098, 0.468600004911423,         &
    0.535290002822876, 0.601350009441376, 0.66482001543045,          &
    0.724009990692139, 0.777729988098145, 0.825269997119904,         & 
    0.866419970989227, 0.901350021362305, 0.930540025234222,         & 
    0.954590022563934, 0.974179983139038, 0.990000009536743, 1/
#ifdef AMIPW_PHYSICS
   jtf=jte
   ktf=kte
   itf=ite
#else
   jtf=min0(jte,jde-1)
   ktf=min0(kte,kde-1)
   itf=min0(ite,ide-1)
#endif

    do k=1,paerlev
      m_hybi(k)=hybi(k)
    enddo

!
! mxaerl = max number of levels (from bottom) for background aerosol
! limit background aerosol height to regions below 900 mb
!

   psurf = 1.e05
   mxaerl = 0
!  do k=pver,1,-1
   do k=kms,kme-1
#ifdef AMIPW_PHYSICS
     if (hypm(k) >= 9.e4) mxaerl = mxaerl + 1
#else
!     if (hypm(k) >= 9.e4) mxaerl = mxaerl + 1
      if (shalf(k)*psurf+pptop  >= 9.e4) mxaerl = mxaerl + 1
#endif
   end do
   mxaerl = max(mxaerl,1)
!  if (masterproc) then
      write(6,*)'aerosols:  background aerosol will be limited to ', &
                'bottom ',mxaerl,' model interfaces.'
!               'bottom ',mxaerl,' model interfaces. top interface is ', &
!               hypi(pverp-mxaerl),' pascals'
!  end if

     do j=jts,jtf
     do i=its,itf
      m_psp(i,j)=psurf
      m_psn(i,j)=psurf
     enddo
     enddo

     do j=jts,jtf
     do i=its,itf
     do k=1,paerlev
! aerosolc arrays are upward cumulative (kg/m2) at each level
! here we assume uniform vertical distribution (aerosolc linear with hybi)
      aerosolcp(i,k,j,idxsul)=1.e-7*(1.-hybi(k))
      aerosolcn(i,k,j,idxsul)=1.e-7*(1.-hybi(k))
      aerosolcp(i,k,j,idxsslt)=1.e-22*(1.-hybi(k))
      aerosolcn(i,k,j,idxsslt)=1.e-22*(1.-hybi(k))
      aerosolcp(i,k,j,idxdustfirst)=1.e-7*(1.-hybi(k))
      aerosolcn(i,k,j,idxdustfirst)=1.e-7*(1.-hybi(k))
      aerosolcp(i,k,j,idxdustfirst+1)=1.e-7*(1.-hybi(k))
      aerosolcn(i,k,j,idxdustfirst+1)=1.e-7*(1.-hybi(k))
      aerosolcp(i,k,j,idxdustfirst+2)=1.e-7*(1.-hybi(k))
      aerosolcn(i,k,j,idxdustfirst+2)=1.e-7*(1.-hybi(k))
      aerosolcp(i,k,j,idxdustfirst+3)=1.e-7*(1.-hybi(k))
      aerosolcn(i,k,j,idxdustfirst+3)=1.e-7*(1.-hybi(k))
      aerosolcp(i,k,j,idxocpho)=1.e-7*(1.-hybi(k))
      aerosolcn(i,k,j,idxocpho)=1.e-7*(1.-hybi(k))
      aerosolcp(i,k,j,idxbcpho)=1.e-9*(1.-hybi(k))
      aerosolcn(i,k,j,idxbcpho)=1.e-9*(1.-hybi(k))
      aerosolcp(i,k,j,idxocphi)=1.e-7*(1.-hybi(k))
      aerosolcn(i,k,j,idxocphi)=1.e-7*(1.-hybi(k))
      aerosolcp(i,k,j,idxbcphi)=1.e-8*(1.-hybi(k))
      aerosolcn(i,k,j,idxbcphi)=1.e-8*(1.-hybi(k))
     enddo
     enddo
     enddo

     call aer_optics_initialize
 

end subroutine aerosol_init

  subroutine aer_optics_initialize
#ifndef AMIPW_PHYSICS
use module_wrf_error
#endif

!   use shr_kind_mod, only: r8 => shr_kind_r8
!   use pmgrid  ! masterproc is here
!   use iofilemod, only: getfil

!#if ( defined spmd )
!    use mpishorthand
!#endif
    implicit none

!   include 'netcdf.inc'


    integer :: nrh_opac  ! number of relative humidity values for opac data
    integer :: nbnd      ! number of spectral bands, should be identical to nspint
    real(r8), parameter :: wgt_sscm = 6.0 / 7.0
    integer :: krh_opac  ! rh index for opac rh grid
    integer :: krh       ! another rh index
    integer :: ksz       ! dust size bin index
    integer :: kbnd      ! band index

    real(r8) :: rh   ! local relative humidity variable

    integer, parameter :: irh=8
    real(r8) :: rh_opac(irh)        ! opac relative humidity grid
    real(r8) :: ksul_opac(irh,nspint)    ! sulfate  extinction
    real(r8) :: wsul_opac(irh,nspint)    !          single scattering albedo
    real(r8) :: gsul_opac(irh,nspint)    !          asymmetry parameter
    real(r8) :: ksslt_opac(irh,nspint)   ! sea-salt
    real(r8) :: wsslt_opac(irh,nspint)
    real(r8) :: gsslt_opac(irh,nspint)
    real(r8) :: kssam_opac(irh,nspint)   ! sea-salt accumulation mode
    real(r8) :: wssam_opac(irh,nspint)
    real(r8) :: gssam_opac(irh,nspint)
    real(r8) :: ksscm_opac(irh,nspint)   ! sea-salt coarse mode
    real(r8) :: wsscm_opac(irh,nspint)
    real(r8) :: gsscm_opac(irh,nspint)
    real(r8) :: kcphil_opac(irh,nspint)  ! hydrophilic organic carbon
    real(r8) :: wcphil_opac(irh,nspint)
    real(r8) :: gcphil_opac(irh,nspint)
    real(r8) :: dummy(nspint)

      logical                 :: opened
      logical , external      :: wrf_dm_on_monitor

      character*80 errmess
      integer cam_aer_unit
      integer :: i
#ifdef AMIPW_PHYSICS
      integer :: comm, ierr
      comm = mpi_comm_world
#endif

!   read aerosol optics data
#ifndef AMIPW_PHYSICS
      if ( wrf_dm_on_monitor() ) then
        do i = 10,99
          inquire ( i , opened = opened )
          if ( .not. opened ) then
            cam_aer_unit = i
            goto 2010
          endif
        enddo
        cam_aer_unit = -1
 2010   continue
      endif
      call wrf_dm_bcast_bytes ( cam_aer_unit , iwordsize )
      if ( cam_aer_unit < 0 ) then
        call wrf_error_fatal ( 'module_ra_cam: aer_optics_initialize: can not find unused fortran unit to read in lookup table.' )
      endif
#endif
      cam_aer_unit = 999

#ifndef AMIPW_PHYSICS
        if ( wrf_dm_on_monitor() ) then
#else
       if(mpi_rank().eq.0)then
#endif
          open(cam_aer_unit,file='CAM_AEROPT_DATA',                  &
               form='unformatted',status='old',err=9010)
#ifndef AMIPW_PHYSICS
          call wrf_debug(50,'reading cam_aeropt_data')
#else
          print*,'reading cam_aeropt_data'
#endif
       endif
#ifndef AMIPW_PHYSICS
#define dm_bcast_macro(a) call wrf_dm_bcast_bytes ( a , size ( a ) * r8 )
#endif
#ifndef AMIPW_PHYSICS
         if ( wrf_dm_on_monitor() ) then
#else
         if ( mpi_rank().eq.0 ) then
#endif
         read (cam_aer_unit,err=9010) dummy
         read (cam_aer_unit,err=9010) rh_opac 
         read (cam_aer_unit,err=9010) ksul_opac 
         read (cam_aer_unit,err=9010) wsul_opac 
         read (cam_aer_unit,err=9010) gsul_opac 
         read (cam_aer_unit,err=9010) kssam_opac 
         read (cam_aer_unit,err=9010) wssam_opac 
         read (cam_aer_unit,err=9010) gssam_opac 
         read (cam_aer_unit,err=9010) ksscm_opac 
         read (cam_aer_unit,err=9010) wsscm_opac 
         read (cam_aer_unit,err=9010) gsscm_opac
         read (cam_aer_unit,err=9010) kcphil_opac 
         read (cam_aer_unit,err=9010) wcphil_opac 
         read (cam_aer_unit,err=9010) gcphil_opac 
         read (cam_aer_unit,err=9010) kcb 
         read (cam_aer_unit,err=9010) wcb 
         read (cam_aer_unit,err=9010) gcb 
         read (cam_aer_unit,err=9010) kdst 
         read (cam_aer_unit,err=9010) wdst 
         read (cam_aer_unit,err=9010) gdst 
         read (cam_aer_unit,err=9010) kbg 
         read (cam_aer_unit,err=9010) wbg 
         read (cam_aer_unit,err=9010) gbg
         read (cam_aer_unit,err=9010) kvolc 
         read (cam_aer_unit,err=9010) wvolc 
         read (cam_aer_unit,err=9010) gvolc
         endif

#ifndef AMIPW_PHYSICS
         dm_bcast_macro(rh_opac)
         dm_bcast_macro(ksul_opac)
         dm_bcast_macro(wsul_opac)
         dm_bcast_macro(gsul_opac)
         dm_bcast_macro(kssam_opac)
         dm_bcast_macro(wssam_opac)
         dm_bcast_macro(gssam_opac)
         dm_bcast_macro(ksscm_opac)
         dm_bcast_macro(wsscm_opac)
         dm_bcast_macro(gsscm_opac)
         dm_bcast_macro(kcphil_opac)
         dm_bcast_macro(wcphil_opac)
         dm_bcast_macro(gcphil_opac)
         dm_bcast_macro(kcb)
         dm_bcast_macro(wcb)
         dm_bcast_macro(gcb)
         dm_bcast_macro(kvolc)
         dm_bcast_macro(wvolc)
         dm_bcast_macro(gvolc)
         dm_bcast_macro(kdst)
         dm_bcast_macro(wdst)
         dm_bcast_macro(gdst)
         dm_bcast_macro(kbg)
         dm_bcast_macro(wbg)
         dm_bcast_macro(gbg)

         if ( wrf_dm_on_monitor() ) close (cam_aer_unit)
#else

         call mpi_bcast(rh_opac,    size(rh_opac)    ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(ksul_opac,  size(ksul_opac)  ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(wsul_opac,  size(wsul_opac)  ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(gsul_opac,  size(gsul_opac)  ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(kssam_opac, size(kssam_opac) ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(wssam_opac, size(wssam_opac) ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(gssam_opac, size(gssam_opac) ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(ksscm_opac, size(ksscm_opac) ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(wsscm_opac, size(wsscm_opac) ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(gsscm_opac, size(gsscm_opac) ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(kcphil_opac,size(kcphil_opac),mpi_real8, 0, comm,ierr)
         call mpi_bcast(wcphil_opac,size(wcphil_opac),mpi_real8, 0, comm,ierr)
         call mpi_bcast(gcphil_opac,size(gcphil_opac),mpi_real8, 0, comm,ierr)
         call mpi_bcast(kcb,        size(kcb)        ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(wcb,        size(wcb)        ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(gcb,        size(gcb)        ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(kvolc,      size(kvolc)      ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(wvolc,      size(wvolc)      ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(gvolc,      size(gvolc)      ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(kdst,       size(kdst)       ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(wdst,       size(wdst)       ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(gdst,       size(gdst)       ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(kbg,        size(kbg)        ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(wbg,        size(wbg)        ,mpi_real8, 0, comm,ierr)
         call mpi_bcast(gbg,        size(gbg)        ,mpi_real8, 0, comm,ierr)

         if ( mpi_rank().eq.0 ) close (cam_aer_unit)
#endif
!-------------------------------------------------

    ! map opac aerosol species onto cam aerosol species
    ! cam name             opac name
    ! sul   or so4         = suso                  sulfate soluble
    ! sslt  or sslt        = 1/7 ssam + 6/7 sscm   sea-salt accumulation/coagulation mode
    ! cphil or cphi        = waso                  water soluble (carbon)
    ! cphob or cpho        = waso @ rh = 0
    ! cb    or bcphi/bcpho = soot

    ksslt_opac(:,:) = (1.0 - wgt_sscm) * kssam_opac(:,:) + wgt_sscm * ksscm_opac(:,:)

    wsslt_opac(:,:) = ( (1.0 - wgt_sscm) * kssam_opac(:,:) * wssam_opac(:,:) &
                  + wgt_sscm * ksscm_opac(:,:) * wsscm_opac(:,:) ) &
                  / ksslt_opac(:,:)

    gsslt_opac(:,:) = ( (1.0 - wgt_sscm) * kssam_opac(:,:) * wssam_opac(:,:) * gssam_opac(:,:) &
                  + wgt_sscm * ksscm_opac(:,:) * wsscm_opac(:,:) * gsscm_opac(:,:) ) &
                   / ( ksslt_opac(:,:) * wsslt_opac(:,:) )

    do i=1,nspint
    kcphob(i) = kcphil_opac(1,i)
    wcphob(i) = wcphil_opac(1,i)
    gcphob(i) = gcphil_opac(1,i)
    end do

    ! interpolate optical properties of hygrospopic aerosol species
    !   onto a uniform relative humidity grid

    nbnd = nspint

    do krh = 1, nrh
      rh = 1.0_r8 / nrh * (krh - 1)
      do kbnd = 1, nbnd
        ksul(krh, kbnd) = exp_interpol( rh_opac, &
          ksul_opac(:, kbnd) / ksul_opac(1, kbnd), rh ) * ksul_opac(1, kbnd)
        wsul(krh, kbnd) = lin_interpol( rh_opac, &
          wsul_opac(:, kbnd) / wsul_opac(1, kbnd), rh ) * wsul_opac(1, kbnd)
        gsul(krh, kbnd) = lin_interpol( rh_opac, &
          gsul_opac(:, kbnd) / gsul_opac(1, kbnd), rh ) * gsul_opac(1, kbnd)
        ksslt(krh, kbnd) = exp_interpol( rh_opac, &
          ksslt_opac(:, kbnd) / ksslt_opac(1, kbnd), rh ) * ksslt_opac(1, kbnd)
        wsslt(krh, kbnd) = lin_interpol( rh_opac, &
          wsslt_opac(:, kbnd) / wsslt_opac(1, kbnd), rh ) * wsslt_opac(1, kbnd)
        gsslt(krh, kbnd) = lin_interpol( rh_opac, &
          gsslt_opac(:, kbnd) / gsslt_opac(1, kbnd), rh ) * gsslt_opac(1, kbnd)
        kcphil(krh, kbnd) = exp_interpol( rh_opac, &
          kcphil_opac(:, kbnd) / kcphil_opac(1, kbnd), rh ) * kcphil_opac(1, kbnd)
        wcphil(krh, kbnd) = lin_interpol( rh_opac, &
          wcphil_opac(:, kbnd) / wcphil_opac(1, kbnd), rh ) * wcphil_opac(1, kbnd)
        gcphil(krh, kbnd) = lin_interpol( rh_opac, &
          gcphil_opac(:, kbnd) / gcphil_opac(1, kbnd), rh )  * gcphil_opac(1, kbnd)
      end do
    end do

     return
9010 continue
     write( errmess , '(a35,i4)' ) 'module_ra_cam: error reading unit ',cam_aer_unit
#ifndef AMIPW_PHYSICS
     call wrf_error_fatal(errmess)
#endif

end subroutine aer_optics_initialize


subroutine radaeini( pstdx, mwdryx, mwco2x )
#ifndef AMIPW_PHYSICS
use module_wrf_error
#endif
!
! initialize radae module data
!
!
! input variables
!
   real(r8), intent(in) :: pstdx   ! standard pressure (dynes/cm^2)
   real(r8), intent(in) :: mwdryx  ! molecular weight of dry air 
   real(r8), intent(in) :: mwco2x  ! molecular weight of carbon dioxide
!
!      variables for loading absorptivity/emissivity
!
   integer ncid_ae                ! netcdf file id for abs/ems file

   integer pdimid                 ! pressure dimension id
   integer psize                  ! pressure dimension size

   integer tpdimid                ! path temperature dimension id
   integer tpsize                 ! path temperature size

   integer tedimid                ! emission temperature dimension id
   integer tesize                 ! emission temperature size

   integer udimid                 ! u (h2o path) dimension id
   integer usize                  ! u (h2o path) dimension size

   integer rhdimid                ! relative humidity dimension id
   integer rhsize                 ! relative humidity dimension size

   integer    ah2onwid            ! var. id for non-wndw abs.
   integer    eh2onwid            ! var. id for non-wndw ems.
   integer    ah2owid             ! var. id for wndw abs. (adjacent layers)
   integer cn_ah2owid             ! var. id for continuum trans. for wndw abs.
   integer cn_eh2owid             ! var. id for continuum trans. for wndw ems.
   integer ln_ah2owid             ! var. id for line trans. for wndw abs.
   integer ln_eh2owid             ! var. id for line trans. for wndw ems.
   
!  character*(nf_max_name) tmpname! dummy variable for var/dim names
   character(len=256) locfn       ! local filename
   integer tmptype                ! dummy variable for variable type
   integer ndims                  ! number of dimensions
!  integer dims(nf_max_var_dims)  ! vector of dimension ids
   integer natt                   ! number of attributes
!
! variables for setting up h2o table
!
   integer t                     ! path temperature
   integer tmin                  ! mininum path temperature
   integer tmax                  ! maximum path temperature
   integer itype                 ! type of sat. pressure (=0 -> h2o only)
   integer i
   real(r8) tdbl

      logical                 :: opened
      logical , external      :: wrf_dm_on_monitor

      character*80 errmess
      integer cam_abs_unit
#ifdef AMIPW_PHYSICS
      integer :: comm, ierr
      comm = mpi_comm_world
#endif

!
! constants to set
!
   p0     = pstdx
   amd    = mwdryx
   amco2  = mwco2x
!
! coefficients for h2o emissivity and absorptivity for overlap of h2o 
!    and trace gases.
!
   c16  = coefj(3,1)/coefj(2,1)
   c17  = coefk(3,1)/coefk(2,1)
   c26  = coefj(3,2)/coefj(2,2)
   c27  = coefk(3,2)/coefk(2,2)
   c28  = .5
   c29  = .002053
   c30  = .1
   c31  = 3.0e-5
!
! initialize further longwave constants referring to far wing
! correction for overlap of h2o and trace gases; r&d refers to:
!
!            ramanathan, v. and  p.downey, 1986: a nonisothermal
!            emissivity and absorptivity formulation for water vapor
!            journal of geophysical research, vol. 91., d8, pp 8649-8666
!
   fwcoef = .1           ! see eq(33) r&d
   fwc1   = .30          ! see eq(33) r&d
   fwc2   = 4.5          ! see eq(33) and eq(34) in r&d
   fc1    = 2.6          ! see eq(34) r&d

#ifndef AMIPW_PHYSICS
      if ( wrf_dm_on_monitor() ) then
        do i = 10,99
          inquire ( i , opened = opened )
          if ( .not. opened ) then
            cam_abs_unit = i
            goto 2010
          endif
        enddo
        cam_abs_unit = -1
 2010   continue
      endif
      call wrf_dm_bcast_bytes ( cam_abs_unit , iwordsize )
      if ( cam_abs_unit < 0 ) then
        call wrf_error_fatal ( 'module_ra_cam: radaeinit: can not find unused fortran unit to read in lookup table.' )
      endif
#endif

#ifdef AMIPW_PHYSICS
        cam_abs_unit = 888 
        if(mpi_rank().eq.0)then
#else
        if ( wrf_dm_on_monitor() ) then
#endif
          open(cam_abs_unit,file='CAM_ABS_DATA',                  &
               form='unformatted',status='old',err=9010)
#ifndef AMIPW_PHYSICS
          call wrf_debug(50,'reading cam_abs_data')
#else
          print*,'reading cam_abs_data'
#endif
        endif

#ifdef AMIPW_PHYSICS
#define dm_bcast_macro(a) call wrf_dm_bcast_bytes ( a , size ( a ) * r8 )
#endif
#ifdef AMIPW_PHYSICS
        if ( mpi_rank().eq.0 ) then
#else
         if ( wrf_dm_on_monitor() ) then
#endif
         read (cam_abs_unit,err=9010) ah2onw
         read (cam_abs_unit,err=9010) eh2onw 
         read (cam_abs_unit,err=9010) ah2ow 
         read (cam_abs_unit,err=9010) cn_ah2ow 
         read (cam_abs_unit,err=9010) cn_eh2ow 
         read (cam_abs_unit,err=9010) ln_ah2ow 
         read (cam_abs_unit,err=9010) ln_eh2ow 

         endif

#ifdef AMIPW_PHYSICS
        call mpi_bcast(ah2onw , size(ah2onw)  , mpi_real8, 0, comm,ierr)
        call mpi_bcast(eh2onw , size(eh2onw)  , mpi_real8, 0, comm,ierr)
        call mpi_bcast(ah2ow  , size(ah2ow)   , mpi_real8, 0, comm,ierr)
        call mpi_bcast(cn_ah2ow,size(cn_ah2ow), mpi_real8, 0, comm,ierr)
        call mpi_bcast(cn_eh2ow,size(cn_eh2ow), mpi_real8, 0, comm,ierr)
        call mpi_bcast(ln_ah2ow,size(ln_ah2ow), mpi_real8, 0, comm,ierr)
        call mpi_bcast(ln_eh2ow,size(ln_eh2ow), mpi_real8, 0, comm,ierr)
        if ( mpi_rank().eq.0 ) close (cam_abs_unit)
#else
         dm_bcast_macro(ah2onw)
         dm_bcast_macro(eh2onw)
         dm_bcast_macro(ah2ow)
         dm_bcast_macro(cn_ah2ow)
         dm_bcast_macro(cn_eh2ow)
         dm_bcast_macro(ln_ah2ow)
         dm_bcast_macro(ln_eh2ow)
         if ( wrf_dm_on_monitor() ) close (cam_abs_unit)
#endif

      
! set up table of h2o saturation vapor pressures for use in calculation
!     effective path rh.  need separate table from table in wv_saturation 
!     because:
!     (1. path temperatures can fall below minimum of that table; and
!     (2. abs/emissivity tables are derived with rh for water only.
!
      tmin = nint(min_tp_h2o)
      tmax = nint(max_tp_h2o)+1
      itype = 0
      do t = tmin, tmax
!        call gffgch(dble(t),estblh2o(t-tmin),itype)
         tdbl = t
         call gffgch(tdbl,estblh2o(t-tmin),itype)
      end do

     return
9010 continue
     write( errmess , '(a35,i4)' ) 'module_ra_cam: error reading unit ',cam_abs_unit
#ifndef AMIPW_PHYSICS
     call wrf_error_fatal(errmess)
#endif
end subroutine radaeini

 
end module module_ra_cam_support
