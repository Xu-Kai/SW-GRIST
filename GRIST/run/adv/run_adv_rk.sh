for glv in  5 6 7 8  ;do
for adv in  5 ;do
for ini in 20 ;do
 
 if [[ ${glv} -eq 4 ]] ;then
    export timestep=400
    export mesh_nv=2562
    export nproc=16
 fi
 if [[ ${glv} -eq 5 ]] ;then
    export timestep=1800
    export mesh_nv=10242
    export nproc=32
 fi
 if [[ ${glv} -eq 6 ]] ;then
    export timestep=900
    export mesh_nv=40962
    export nproc=64
 fi
 if [[ ${glv} -eq 7 ]] ;then
    export timestep=450
    export mesh_nv=163842
    export nproc=64
 fi
 if [[ ${glv} -eq 8 ]] ;then
    export timestep=225
    export mesh_nv=655362
    export nproc=256
 fi
export glevel=g${glv}
export time_scheme=rk3
export testcase=20
export advection_scheme=${adv}
export pv_order=2
export isolate_adv=.true.
export test_dual_cell=.false.
export limiter=.false.
export mypath=`pwd`
export outdir=${mypath}/def-${glevel}-tc${testcase}-ini${ini}-${time_scheme}-adv${advection_scheme}-pv${pv_order}/
mkdir ${outdir}
cd ${outdir}
cp /g13/zhangyi/mac/run/bin/ParGRIST-swe.exe par.exe

cat>run.sbatch<<EOF
#!/usr/bin/bash
#SBATCH --comment=GRAPES
#SBATCH -J GRAPES
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
 day_duration           = 12 
 model_timestep         = ${timestep}
 h1_history_freq        = `expr 86400 / ${timestep}`
 h1_restart_freq        = 9999999
 run_type               = 'init'
/
&data_para
 outdir                 = '${outdir}'
 gridFilePath           = '/g13/zhangyi/mac/run/MPI_SCVT_UR/G${glv}/'
 gridFileNameHead       = 'grist.grid_file.g${glv}.ccvt'
/
&swe_para
 swe_timestep           = ${timestep}
 testcase               = ${testcase}  ! swe
 initialfield           = ${ini}       ! swe
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
 scalar_type            = 'cell'
 use_streamf            = .true.
 limiter_type           = 'fct'
/
&dycore_para
/
&tracer_para
/
&physics_para
/
&mesh_para
 mesh_nv                = ${mesh_nv} 
 mesh_glv               = ${glv}
 nlev                   = 18
 nlevp                  = 19
/
&ccvt_para
/
EOF
sbatch run.sbatch
cd ${mypath}
done
done
done
