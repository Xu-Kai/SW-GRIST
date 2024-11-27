export glevel=g4
export time_scheme=rk3
export testcase=6
export timestep=200
export advection_scheme=3
export pv_order=3
export isolate_adv=.false.
export test_dual_cell=.false.
export limiter=.false.
export mypath=`pwd`
export outdir=${mypath}/1-2021swe-kite3-${glevel}-tc${testcase}-${time_scheme}-adv${advection_scheme}-pv${pv_order}/
mkdir ${outdir}

cd ${outdir}
cp /Users/zhangyi/grist/run/par-grist-swe.exe ./

cat>grist.nml<<EOF
&ctl_para
 model_timestep         = ${timestep}
 day_duration           = 12.
 outdir                 = '${outdir}'
 grddir                 = '/Users/zhangyi/grist/run/MeshFiles/2017-1104-ccvt-3000/'
 h1_history_freq        = `expr 86400 / ${timestep}`
 restart_freq           = `expr 86400 / ${timestep}`
 write_history_h1       = .true.
 run_type               = 'init'
 gridFileNameHead       = "grist.grid_file.g6.ccvt"
/
&data_para
/
&swe_para
 swe_timestep           = ${timestep}
 testcase               = ${testcase}  ! swe
 initialfield           = ${testcase}  ! swe
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
 ke_method              = 0
 ke_alpha               = 0.375
 mass_wtd_vor           = .true.
 modon_radius           = 500000.
 ppm_coef               = 0.25
 ffsl_alpha             = 1.
 scalar_type            = 'point'
/
&dycore_para
/
&tracer_para
/
&physics_para
/
&mesh_para
 mesh_nv                = 2562
 mesh_glv               = 4
 nlev                   = 18
 nlevp                  = 19
/
&ccvt_para
/
EOF
t=$(date +"%H-%M-%S")
runnum=4
nohup mpirun -np ${runnum} ./par-grist-swe.exe > log.${t} &
