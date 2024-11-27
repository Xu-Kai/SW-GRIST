export year=2012
export mon=05
export day=26
export hadv=33
export hnrk=3
export vadv=3
export vnrk=3
export glevel=g6
export glv=6
export ncell=40962
export model_timestep=1200
export dycore_timestep=300
export tracer_timestep=600
export nproc=64
export nnode=2
export working_mode=ampwrf
export testcase=6
export advection_scheme=2
export pv_order=3
export pt_adv_flag=33
export isolate_adv=.false.
export test_dual_cell=.false.
export limiter=.false.
export mypath=`pwd`
export outdir=${mypath}/Latest-sponge2-pv3-rad2h-rel105-${glevel}-hadv${hadv}-hnrk${hnrk}-vadv${vadv}-vnrk${vnrk}/
mkdir ${outdir}
cd ${outdir}

cp /THL8/home/zhangyi/grist/run/bin/ParGRIST-amp-wrf-real.exe    ${outdir}/par.exe
cp /THL8/home/zhangyi/grist/inputdata/wrf-data/RRTM_DATA         ${outdir}
cp /THL8/home/zhangyi/grist/inputdata/wrf-data/RRTMG*            ${outdir}
cp /THL8/home/zhangyi/grist/inputdata/wrf-data/camradData/CAM_A* ${outdir}
cp /THL8/home/zhangyi/grist/inputdata/wrf-data/camradData/ozone* ${outdir}
cp /THL8/home/zhangyi/grist/inputdata/noahmp_data/*              ${outdir}

cat>run.sh<<EOF
#!/bin/sh
export LD_LIBRARY_PATH=/THL8/home/zhangyi/software/netcdf-4.0.1/lib/:${LD_LIBRARY_PATH}
export LD_LIBRARY_PATH=/THL8/home/zhangyi/software/metis-5.1.0/build/Linux-x86_64/libmetis/:${LD_LIBRARY_PATH}
yhrun -N ${nnode} -n ${nproc} par.exe
EOF
cat>batch.sh<<EOF
yhbatch -N ${nnode} -n ${nproc} -p debug3N run.sh
EOF
chmod u+x run.sh batch.sh

cat>grist.nml<<EOF
&ctl_para
 day_duration           = 150
 model_timestep         = ${model_timestep}
 h1_history_freq        = `expr 86400 / ${model_timestep}`
 h1_restart_freq        = `expr 8640000 / ${model_timestep}`
 run_type               = 'init'
 working_mode           = '${working_mode}'
 write_history          = .true.
 write_history_h0       = .true.
 write_restart_h0       = .true.
 test_real_case         = .true.
 doAquaPlanet           = .false.
 start_ymd              = 20000526
 start_tod              = 0
 comm_group_size        = 10
 grid_info              = "G6"
/
&data_para
 outdir                 = '${outdir}'
 gridFilePath           = '/THL8/home/zhangyi/grist/inputdata/mesh/uniform-pol/'
 gridFileNameHead       = 'grist.grid_file.g6.ccvt'
 staticFilePath         = '/THL8/home/zhangyi/grist/inputdata/init/initialData/20120526/uniform-g6/static_g6.nc'
 sstFilePath            = '/THL8/home/zhangyi/grist/inputdata/bdy/new/g6/'
 initialAtmFilePath     = '/THL8/home/zhangyi/grist/inputdata/init/initialData/20120526/uniform-g6/grist.initial.pl.g6_20120526.nc.new.nc'
 initialLndFilePath     = '/THL8/home/zhangyi/grist/inputdata/init/initialData/20120526/uniform-g6/grist_gfs_20120526.g6.nc'
 initialDataSorc        = 'ERAIP'
 sstFile_year_beg       = 2000
 real_cycle_style       = 'CYCLE'
 numMonSST              = 12
 sstFileNameHead        = 'realNoMissingNewSstSic.' 
 sstFileNameTail        = '.grist.g6.nc'
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
 pv_order               = ${pv_order}
 generate_cdo_file      = .true.
 use_limiter_rk3        = ${limiter}
 use_tr_mbs             = .false.
 use_tspas              = .false.
 ke_method              = 1
 ke_alpha               = 0.9
 mass_wtd_vor           = .true.
 modon_radius           = 500000.
 ppm_coef               = 0.25
 ffsl_alpha             = 1.
/
&dycore_para
 gcm_testcase           = 'real-ERAIP'
 mas_adv_flag           = 2
 pot_adv_flag           = 33
 phi_adv_flag           = 33
 www_adv_flag           = 33
 ver_adv_flag           = 3
 hor_pgf_flag           = 6
 nrk                    = 3
 nh_dynamics            = .false.
 www_damping_coef       = 0.2
 d3d_damping_coef       = 0.12
 emf_damping_coef       = 0.
 zd                     = 12000.
 imbeta                 = 0.55
 imkesi                 = 0.55
 smg_coef               = 0.015
 ko4_coef               = 2e14
 ref_leng               = 7000
 do_dycore_laplacian_2nd= .true.
 do_dycore_laplacian_4th= .false.
 ad_dycore_laplacian_2nd= .true.
 ad_dycore_laplacian_4th= .false.
 eqs_vert_diff          = .true.
 nexpdif_level          = 30
 nsponge                = 5
 spcoef                 = 2
 use_expdif_tstep       = .true.
 profile_type           = 'tp1'
 dycore_timestep        = ${dycore_timestep}
 vr_mode                = .false.
/
&tracer_para
 ntracer                = 6
 nmif                   = 3
 mif_index              = 1,2,3
 nrk_hori               = ${hnrk}
 nrk_vert               = ${vnrk}
 tracer_hadv_flag       = ${hadv}
 tracer_vadv_flag       = ${vadv}
 tracer_hori_limiter    = .true.
 tracer_vert_limiter    = .true.
 isolate_tracer_test    = .false.
 tracer_timestep        = ${tracer_timestep}
 tracer_hori_timestep   = ${tracer_timestep}
 tracer_vert_timestep   = ${tracer_timestep}
 do_tracer_laplacian_2nd= .false.
 tracer_hvsplit_method  = 1
/
&physics_para
 use_phys               = .true.
 physpkg                = 'WRF2_PHYSICS'
 sub_physpkg            = 'none'
 ptend_wind_rk_on       = .true.
 ptend_heat_rk_on       = .true.
 ptend_dycore_f1_on     = .false.
 ptend_tracer_f1_on     = .false.
 ptend_f2_on            = .true.
 ptendSubDiabPhys       = .false.
 physics_coupling       = 'P3'
 TC_pbl_flag            = .false.
 nspecies               = 6
/
&mesh_para
 mesh_nv                = ${ncell}
 mesh_glv               = ${glv}
 nlev                   = 30
 nlevp                  = 31
 nlev_inidata           = 37
 levsoil                = 4
/
&ccvt_para
/
EOF
cat>grist_wrfphys.nml<<EOF
&wrfphys_para
 wrfphys_cu_scheme      = 'NTDKV381'
 wrfphys_cf_scheme      = 'CAM3'
 wrfphys_ra_scheme      = 'RRTMGV381'
 wrfphys_rasw_scheme    = 'RRTMGV381'
 wrfphys_ralw_scheme    = 'RRTMGV381'
 wrfphys_mp_scheme      = 'WSM6V381'
 wrfphys_bl_scheme      = 'YSUV381'
 wrfphys_sf_scheme      = 'SFCLAYV381'
 wrfphys_lm_scheme      = 'noahmp'
 wphys_has_req          = 0
 step_cu                = 1
 step_ra                = 6
 use_gwdo               = .false.
/
EOF
#./batch.sh
