export year=2021
export mon=06

for day in 01 ;do
export hadv=33
export hnrk=3
export vadv=3
export vnrk=3
export glevel=g8b3
export glv=9
export ncell=5898242
export model_timestep=180
export dycore_timestep=15
export tracer_timestep=90
export nproc=1600
export nnode=50
export working_mode=dtp
export testcase=6
export advection_scheme=2
export pv_order=2
export pt_adv_flag=33
export isolate_adv=.false.
export test_dual_cell=.false.
export limiter=.false.
export mypath=`pwd`
export outdir=${mypath}/run-g8b3-Beg${year}${mon}${day}-hadv${hadv}-hnrk${hnrk}-vadv${vadv}-vnrk${vnrk}/
mkdir ${outdir}
cd ${outdir}

cp /THL8/home/zhangyi/grist/run/bin/ParGRIST-dtp-dcmip2016-tc.exe  ${outdir}/par.exe
#cp /THL8/home/zhangyi/grist/inputdata/wrf-data/RRTM_DATA          ${outdir}
#cp /THL8/home/zhangyi/grist/inputdata/wrf-data/RRTMG*             ${outdir}
#cp /THL8/home/zhangyi/grist/inputdata/wrf-data/camradData/CAM_A*  ${outdir}
#cp /THL8/home/zhangyi/grist/inputdata/wrf-data/camradData/ozone*  ${outdir}
#cp /THL8/home/zhangyi/grist/inputdata/noahmp_data/*               ${outdir}
#cp /THL8/home/zhangyi/grist/run/CPTP_T1/config/grist_lsm_noahmp.nml    ${outdir}
cp ../run.sh   ${outdir}
cp ../batch.sh ${outdir}

# platform dependant scripts
cat>grist.nml<<EOF
&ctl_para
 day_duration           = 10
 model_timestep         = ${model_timestep}
 h1_history_freq        = `expr 10800 / ${model_timestep}`
 h1_restart_freq        = `expr 86400 / ${model_timestep}`
 run_type               = 'init'
 working_mode           = '${working_mode}'
 write_history          = .true.
 write_history_h1       = .true.
 write_restart_h1       = .false.
 write_stepinfo         = .true.
 test_real_case         = .false.
 doAquaPlanet           = .false.
 start_ymd              = ${year}${mon}${day}
 start_tod              = 0
 comm_group_size        = 10
 grid_info              = "G8B3UR"
/
&data_para
 outdir                 = '${outdir}'
 gridFilePath           = '/THL8/home/zhangyi/grist/inputdata/mesh/G8B3_grid/'
 gridFileNameHead       = 'grist.grid_file.g8b3.ccvt'
 staticFilePath         = '/THL8/home/zhangyi/public/GRIST/data/static/static_uniform_g9.nc'
 initialAtmFilePath     = '/THL8/home/zhangyi/public/GRIST/data/GFS/geniniFromGFS/data/init/grist.gfs.initial.pl.G9UR_${year}${mon}${day}.00.nc'
 initialLndFilePath     = '/THL8/home/zhangyi/public/GRIST/data/GFS/geniniFromGFS/data/init/grist.gfs.initial.pl.G9UR_${year}${mon}${day}.00.nc'
 sstFilePath            = '/THL8/home/zhangyi/public/GRIST/data/GFS/geniniFromGFS/data/sst/'
 initialDataSorc        = 'ERAIP'
 numMonSST              = 12
 sstFile_year_beg       = ${year}
 real_sst_style         = 'CYCLE'
 sstFileNameHead        = 'realNoMissGFSSstSic${year}${mon}${day}.'
 sstFileNameTail        = '.GRIST.2621442.nc'
/
&swe_para
 testcase               = ${testcase}
 initialfield           = ${testcase}
 nsplit                 = 6
 time_scheme            = '${time_scheme}'
 advection_scheme       = ${advection_scheme}
 isolate_advection_test = ${isolate_adv}
 conserve_scheme        = 'te'
 test_dual_cell         = ${test_dual_cell}
 pv_order               = 2,2,3
 generate_cdo_file      = .true.
 use_limiter_rk3        = ${limiter}
 use_tr_mbs             = .false.
 use_tspas              = .false.
 ke_method              = 1
 ke_alpha               = 0.9
 mass_wtd_vor           = .true.
 modon_radius           = 500000.
/
&dycore_para
 gcm_testcase           = 'DCMIP2016-TC'
 mas_adv_flag           = 2
 pot_adv_flag           = 33,33,33
 phi_adv_flag           = 33,33,33
 www_adv_flag           = 33,33,33
 ver_adv_flag           = 3
 hor_pgf_flag           = 6
 nrk                    = 3
 nh_dynamics            = .true.
 www_damping_coef       = 0.0
 d3d_damping_coef       = 0.12
 zd                     = 12000.
 imbeta                 = 0.55
 imkesi                 = 0.55
 smg_coef               = 0.01
 ko4_coef               = 1e11
 ref_leng               = 10000
 do_dycore_laplacian_2nd= .true.
 do_dycore_laplacian_4th= .false.
 ad_dycore_laplacian_2nd= .true.
 ad_dycore_laplacian_4th= .false.
 eqs_vert_diff          = .true.
 nexpdif_level          = 30
 nsponge                = 5
 spcoef                 = 0
 use_expdif_tstep       = .false.
 profile_type           = 'tp1'
 dycore_timestep        = ${dycore_timestep}
 vr_mode                = .false.
 index_flag             = 'bfs'
 www_k2coef_type        = 'SMG'
 smooth_topo            = .false.
 nsmooth_topo           = 12
 smooth_type            = 'cellAvg'
 topo_type              = 'static'
 doNotDiagnose          = .true.
/
&tracer_para
 ntracer                = 1
 nmif                   = 1
 mif_index              = 1
 isolate_tracer_test    = .false.
 tracer_timestep        = ${tracer_timestep}
 tracer_hori_timestep   = ${tracer_timestep}
 tracer_vert_timestep   = ${tracer_timestep}
 do_tracer_laplacian_2nd= .false.
 tracer_hvsplit_method  = 1
 tracer_transport_hori_scheme = 'rk3o3'
 tracer_transport_vert_scheme = 'rk3o3-aimp'
/
&physics_para
 use_phys               = .true.
 physpkg                = 'DCMIP2016-TC'
 sub_physpkg            = 'none'
 ptend_wind_rk_on       = .true.
 ptend_heat_rk_on       = .true.
 ptend_dycore_f1_on     = .false.
 ptend_tracer_f1_on     = .false.
 ptend_f2_on            = .true.
 ptendSubDiabPhys       = .true.
 physics_coupling       = 'P3'
 use_som                = .false.
 nspecies               = 1
/
&mesh_para
 mesh_nv                = ${ncell}
 mesh_glv               = ${glv}
 nlev                   = 30
 nlevp                  = 31
 nlev_inidata           = 41
 levsoil                = 4
/
&ccvt_para
/
EOF
#cat>grist_amipw_phys.nml<<EOF
#&wrfphys_para
# wrfphys_cu_scheme      = 'NTDKV381'
# wrfphys_cf_scheme      = 'RANDALL'
# wrfphys_ra_scheme      = 'RRTMGV381'
# wrfphys_rasw_scheme    = 'RRTMGV381'
# wrfphys_ralw_scheme    = 'RRTMGV381'
# wrfphys_mp_scheme      = 'WSM6V381'
# wrfphys_bl_scheme      = 'YSUV381'
# wrfphys_sf_scheme      = 'SFCLAYV381'
# wrfphys_lm_scheme      = 'noahmp'
# wphys_has_req          = 1
# unuse_cu               = .false.
# step_cu                = 5
# step_ra                = 5
# use_gwdo               = .true.
#/
#EOF
#./batch.sh
cd ../
done
