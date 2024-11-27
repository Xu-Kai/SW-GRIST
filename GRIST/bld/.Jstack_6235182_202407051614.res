
==============================
  T0 -> 0.0-127.63
  T1 -> 2-31,33-39,41-43,46-47,50-51,65-67,70-71,73-75,77-95,97-111
  T2 -> 0,32,40,44,64,72,96,116
  T3 -> 1,45,48-49,52,54,56-57,59,68-69,76,113-114,117-118
  T4 -> 60,62,120,123-126
  T5 -> 127
  T6 -> 53,55,58,61,63,112,115,119,121-122
------------------------------
 -0- /home/export/online1/mdt00/shisuan/swgbcm/xukai/grist_bencmark/demo-g6-aqua/../GRIST-main/bld/ParGRIST-GCM_AMIPW_MIXED-mixed_sunway_xfort.exe 
   -1- slave__Waiting_For_Task ([0/8192])
     -2- slave_jobserver at jobserver-slave.S:83 ([8192/8192] T0 6.63/v84351 6.62/v84351 6.61/v84351...)
   -1- main at ../../src/atmosphere/gcm/grist_atmos.F90:24 ([0/128])
     -2- grist_atmos at ../../src/atmosphere/gcm/grist_atmos.F90:155 ([0/128])
       -3- grist_gcm_control_driver::grist_gcm_run at ../../src/atmosphere/gcm/grist_gcm_control_driver.F90:414 ([0/128])
         -4- grist_gcm_control_driver::grist_gcm_control_run_amipw at ../../src/atmosphere/gcm/grist_gcm_control_driver.F90:897 ([0/127])
           -5- grist_dycore_time_integration_2d_mixed::dycore_time_integration_run at ../../src/atmosphere/dynamics_mixed/dycore/grist_dycore_time_integration_2d_mixed.F90:701 ([0/110])
             -6- grist_nh_driver_module_mixed::grist_nh_dynamics_run at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:261 ([0/86])
               -7- grist_nh_driver_module_mixed::grist_nh_dynamics_run_diagvars at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:764 ([0/86])
                 -8- pmpi_bcast__ ([0/86])
                   -9- PMPI_Bcast ([86/86] T1 6/v84351 7/v84351 8/v84351...)
             -6- grist_nh_driver_module_mixed::grist_nh_dynamics_run at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:165 ([0/24])
               -7- grist_nh_driver_module_mixed::grist_nh_dynamics_run_explicit at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:377 ([0/8])
                 -8- grist_mpi::reduce_real8_1 at ../../src/grist_lib/src/grist_mpi.F90:126 ([0/8])
                   -9- pmpi_reduce__ ([0/8])
                     -10- PMPI_Reduce ([8/8] T2 0/v84350 116/v84369 44/v84357...)
               -7- grist_nh_driver_module_mixed::grist_nh_dynamics_run_explicit at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:355 ([0/16])
                 -8- grist_nh_explicit_tend_module_2d_mixed::grist_nh_et_adv_face at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_explicit_tend_module_2d_mixed.F90:135 ([0/16])
                   -9- grist_config_partition::exchange_data_2d_r4 at ../../src/infrastructure/grid/grist_config_partition.F90:3134 ([0/16])
                     -10- pmpi_waitall__ ([0/16])
                       -11- PMPI_Waitall ([16/16] T3 1/v84350 114/v84369 117/v84369...)
           -5- grist_dycore_time_integration_2d_mixed::dycore_time_integration_run at ../../src/atmosphere/dynamics_mixed/dycore/grist_dycore_time_integration_2d_mixed.F90:351 ([0/7])
             -6- grist_config_partition::exchange_data_1d at ../../src/infrastructure/grid/grist_config_partition.F90:2279 ([0/7])
               -7- pmpi_waitall__ ([0/7])
                 -8- PMPI_Waitall ([7/7] T4 126/v84371 60/v84360 62/v84360...)
           -5- grist_dycore_time_integration_2d_mixed::dycore_time_integration_run at ../../src/atmosphere/dynamics_mixed/dycore/grist_dycore_time_integration_2d_mixed.F90:398 ([0/10])
             -6- grist_config_partition::exchange_data_2d at ../../src/infrastructure/grid/grist_config_partition.F90:2869 ([0/10])
               -7- pmpi_waitall__ ([0/10])
                 -8- PMPI_Waitall ([10/10] T6 115/v84369 119/v84369 55/v84359...)
         -4- grist_gcm_control_driver::grist_gcm_control_run_amipw at ../../src/atmosphere/gcm/grist_gcm_control_driver.F90:870 ([0/1])
           -5- cesm_swlu_prof_init_ at ../../src/infrastructure/utils/swlu_prof_wrapper.c:42 ([0/1])
             -6- swlu_prof_init at samprof.c:1143 ([0/1])
               -7- memset at ../sysdeps/sw_64/sw_64sw6a/memset.S:202 ([0/1])
                 -8- signal handler called ([0/1])
                   -9- samprof at samprof.c:976 ([0/1])
                     -10- samprof_real at samprof.c:836 ([0/1])
                       -11- signal handler called ([0/1])
                         -12- swch_catchsig ([1/1] T5 127/v84371)
==============================

==========================================================
node(taskid):svrstart,wait,sout
vn084350(0       ):	0.55	10.56	10.56
vn084351(6       ):	0.53	10.54	10.54
vn084352(12      ):	0.56	10.57	10.58
vn084353(18      ):	0.46	10.46	10.47
vn084354(24      ):	0.59	10.59	10.59
vn084355(30      ):	0.55	10.56	10.56
vn084356(36      ):	0.54	10.55	10.57
vn084357(42      ):	0.52	10.52	10.53
vn084358(48      ):	0.55	10.58	10.58
vn084359(54      ):	0.49	10.50	10.51
vn084360(60      ):	0.57	10.57	10.58
vn084361(66      ):	0.56	10.57	10.57
vn084362(72      ):	0.54	10.55	10.56
vn084363(78      ):	0.51	10.52	10.52
vn084364(84      ):	0.53	10.53	10.54
vn084365(90      ):	0.44	10.45	10.45
vn084366(96      ):	0.55	10.55	10.55
vn084367(102     ):	0.54	10.54	10.55
vn084368(108     ):	0.54	10.55	10.57
vn084369(114     ):	0.50	10.50	10.51
vn084370(120     ):	0.58	10.59	10.59
vn084371(126     ):	0.44	10.45	10.45
before:0.005664,scan:31.008309,first_data:30.006974,process:1.001391,show:0.022218,total:31.036247

