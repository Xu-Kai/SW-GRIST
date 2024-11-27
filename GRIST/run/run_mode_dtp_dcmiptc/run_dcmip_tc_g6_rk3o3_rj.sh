export hadv=33
export hnrk=3
export vadv=3
export vnrk=3
export glevel=g6
export glv=6
export ncell=40962
export model_timestep=1200
export dycore_timestep=150
export tracer_timestep=600
export nproc=160
export working_mode=dtp
export testcase=6
export advection_scheme=2
export pv_order=3
export pt_adv_flag=33
export isolate_adv=.false.
export test_dual_cell=.false.
export limiter=.false.
export mypath=`pwd`
export outdir=${mypath}/reg20210807-nhrrrimp-rk3o3-dbio-dcmiptc-${glevel}-hadv${hadv}-hnrk${hnrk}-vadv${vadv}-vnrk${vnrk}/
mkdir ${outdir}
cd ${outdir}
cp /g13/zhangyi/mac/run/bin/ParGRIST-dycore-btbc.exe ${outdir}/par.exe

cat>run.sbatch<<EOF
#!/usr/bin/bash
#SBATCH --comment=CAM
#SBATCH -J CAM
#SBATCH -n ${nproc}
#SBATCH -p normal 
#SBATCH -o gcm_%j.out
#SBATCH -e gcm_%j.err

##set runtime environment variables

ulimit -s unlimited
ulimit -c unlimited

module load compiler/intel/composer_xe_2017.2.174
module load mpi/intelmpi/2017.2.174
export I_MPI_PMI_LIBRARY=/opt/gridview/slurm17/lib/libpmi.so
export LD_LIBRARY_PATH=/g13/zhangyi/softwares/intel2017/metis-5.1.0/build/Linux-x86_64/libmetis/:${LD_LIBRARY_PATH}
srun ./par.exe

EOF

cat>grist.nml<<EOF
&ctl_para
 day_duration           = 30
 model_timestep         = ${model_timestep}
 h1_history_freq        = `expr 86400 / ${model_timestep}`
 h1_restart_freq        = `expr 86400 / ${model_timestep}`
 run_type               = 'init'
 working_mode           = '${working_mode}'
 write_history          = .false.
 write_history_h1       = .true.
 write_restart_h1       = .false.
 comm_group_size        = 1
/
&data_para
 outdir                 = '${outdir}'
 gridFilePath           = '/g13/zhangyi/mac/inputdata/mesh/uniform-pol/'
 gridFileNameHead       = 'grist.grid_file.g6.ccvt'
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
 pv_order               = ${pv_order},${pv_order},${pv_order}
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
 gcm_testcase           = 'DCMIP2016-TC'
 mas_adv_flag           = 2
 pot_adv_flag           = 33,33,33
 phi_adv_flag           = 33,33,33
 www_adv_flag           = 33,33,33
 ver_adv_flag           = 3
 hor_pgf_flag           = 6
 nrk                    = 3
 nh_dynamics            = .true.
 www_damping_coef       = 0.
 d3d_damping_coef       = 0.1
 emf_damping_coef       = 0.
 zd                     = 0.
 imbeta                 = 0.55
 imkesi                 = 0.55
 smg_coef               = 0.015
 do_dycore_laplacian_2nd= .true.
 ad_dycore_laplacian_2nd= .true.
 eqs_vert_diff          = .true.
 nexpdif_level          = 30
 profile_type           = 'tp1'
 dycore_timestep        = ${dycore_timestep}
/
&tracer_para
 ntracer                = 1
 nmif                   = 1
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
 do_tracer_laplacian_2nd= .true.
 tracer_hvsplit_method  = 1
 tracer_mxrt_fixer      = .false.
/
&physics_para
 use_phys               = .true.
 physpkg                = 'DCMIP2016-TC'
 ptend_wind_rk_on       = .true.
 ptend_heat_rk_on       = .true.
 ptend_dycore_f1_on     = .false.
 ptend_tracer_f1_on     = .true.
 ptend_f2_on            = .false.
 ptendSubDiabPhys       = .false.
 TC_pbl_flag            = .false.
/
&mesh_para
 mesh_nv                = ${ncell}
 mesh_glv               = ${glv}
 nlev                   = 30
 nlevp                  = 31
/
&ccvt_para
/
EOF
#sbatch run.sbatch
