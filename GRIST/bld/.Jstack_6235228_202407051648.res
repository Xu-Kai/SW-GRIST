
==============================
  T0 -> 0.0-127.63
  T1 -> 1-31,33-39,41,47-63,85,89,97-103,105-111,113-115,123-127
  T2 -> 0,32,40,42,84,86,90,96,112,120
  T3 -> 43-46,80,82,87-88,91,94-95,104,116,118-119,121-122
  T4 -> 64,67-70,76,78-79,81,83,92-93,117
  T5 -> 65-66,71-72,74-75,77
  T6 -> 73
------------------------------
 -0- /home/export/online1/mdt00/shisuan/swgbcm/xukai/grist_bencmark/demo-g6-aqua/../GRIST-main/bld/ParGRIST-GCM_AMIPW_MIXED-mixed_sunway_xfort.exe 
   -1- slave__Waiting_For_Task ([8192/8192] T0 126.63/v84735 126.62/v84735 126.61/v84735...)
   -1- main at ../../src/atmosphere/gcm/grist_atmos.F90:24 ([0/128])
     -2- grist_atmos at ../../src/atmosphere/gcm/grist_atmos.F90:155 ([0/128])
       -3- grist_gcm_control_driver::grist_gcm_run at ../../src/atmosphere/gcm/grist_gcm_control_driver.F90:414 ([0/128])
         -4- grist_gcm_control_driver::grist_gcm_control_run_amipw at ../../src/atmosphere/gcm/grist_gcm_control_driver.F90:897 ([0/127])
           -5- grist_dycore_time_integration_2d_mixed::dycore_time_integration_run at ../../src/atmosphere/dynamics_mixed/dycore/grist_dycore_time_integration_2d_mixed.F90:701 ([0/107])
             -6- grist_nh_driver_module_mixed::grist_nh_dynamics_run at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:261 ([0/80])
               -7- grist_nh_driver_module_mixed::grist_nh_dynamics_run_diagvars at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:764 ([0/80])
                 -8- pmpi_bcast__ ([0/80])
                   -9- PMPI_Bcast ([80/80] T1 126/v84735 127/v84735 30/v84719...)
             -6- grist_nh_driver_module_mixed::grist_nh_dynamics_run at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:165 ([0/27])
               -7- grist_nh_driver_module_mixed::grist_nh_dynamics_run_explicit at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:377 ([0/10])
                 -8- grist_mpi::reduce_real8_1 at ../../src/grist_lib/src/grist_mpi.F90:126 ([0/10])
                   -9- pmpi_reduce__ ([0/10])
                     -10- PMPI_Reduce ([10/10] T2 32/v84719 0/v84714 96/v84730...)
               -7- grist_nh_driver_module_mixed::grist_nh_dynamics_run_explicit at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_driver_module_mixed.F90:355 ([0/17])
                 -8- grist_nh_explicit_tend_module_2d_mixed::grist_nh_et_adv_face at ../../src/atmosphere/dynamics_mixed/dycore/nhd/grist_nh_explicit_tend_module_2d_mixed.F90:135 ([0/17])
                   -9- grist_config_partition::exchange_data_2d_r4 at ../../src/infrastructure/grid/grist_config_partition.F90:3134 ([0/17])
                     -10- pmpi_waitall__ ([0/17])
                       -11- PMPI_Waitall ([17/17] T3 104/v84731 121/v84734 122/v84734...)
           -5- grist_dycore_time_integration_2d_mixed::dycore_time_integration_run at ../../src/atmosphere/dynamics_mixed/dycore/grist_dycore_time_integration_2d_mixed.F90:398 ([0/13])
             -6- grist_config_partition::exchange_data_2d at ../../src/infrastructure/grid/grist_config_partition.F90:2869 ([0/13])
               -7- pmpi_waitall__ ([0/13])
                 -8- PMPI_Waitall ([13/13] T4 64/v84724 92/v84729 93/v84729...)
           -5- grist_dycore_time_integration_2d_mixed::dycore_time_integration_run at ../../src/atmosphere/dynamics_mixed/dycore/grist_dycore_time_integration_2d_mixed.F90:351 ([0/7])
             -6- grist_config_partition::exchange_data_1d at ../../src/infrastructure/grid/grist_config_partition.F90:2279 ([0/7])
               -7- pmpi_waitall__ ([0/7])
                 -8- PMPI_Waitall ([7/7] T5 65/v84724 72/v84726 74/v84726...)
         -4- grist_gcm_control_driver::grist_gcm_control_run_amipw at ../../src/atmosphere/gcm/grist_gcm_control_driver.F90:870 ([0/1])
           -5- cesm_swlu_prof_init_ at ../../src/infrastructure/utils/swlu_prof_wrapper.c:42 ([0/1])
             -6- swlu_prof_init at samprof.c:1142 ([0/1])
               -7- memset at ../sysdeps/sw_64/sw_64sw6a/memset.S:202 ([0/1])
                 -8- signal handler called ([0/1])
                   -9- samprof at samprof.c:976 ([0/1])
                     -10- samprof_real at samprof.c:805 ([0/1])
                       -11- signal handler called ([0/1])
                         -12- swch_catchsig ([1/1] T6 73/v84726)
==============================

==========================================================
node(taskid):svrstart,wait,sout
vn084714(0       ):	0.49	10.49	10.50
vn084715(6       ):	0.47	10.47	10.47
vn084716(12      ):	0.47	10.48	10.48
vn084717(18      ):	0.51	10.51	10.52
vn084718(24      ):	0.53	10.54	10.54
vn084719(30      ):	0.45	10.46	10.46
vn084720(36      ):	0.51	10.51	10.54
vn084721(42      ):	0.52	10.53	10.53
vn084722(48      ):	0.50	10.51	10.51
vn084723(54      ):	0.50	10.51	10.51
vn084724(60      ):	0.50	10.51	10.51
vn084725(66      ):	0.56	10.58	10.59
vn084726(72      ):	0.53	10.55	10.55
vn084727(78      ):	0.49	10.50	10.51
vn084728(84      ):	0.54	10.55	10.57
vn084729(90      ):	0.50	10.51	10.51
vn084730(96      ):	0.50	10.50	10.50
vn084731(102     ):	0.48	10.49	10.49
vn084732(108     ):	0.50	10.51	10.51
vn084733(114     ):	0.51	10.52	10.53
vn084734(120     ):	0.50	10.50	10.51
vn084735(126     ):	0.39	10.40	10.41
before:0.005284,scan:31.007767,first_data:30.007472,process:1.000383,show:0.014150,total:31.027289

