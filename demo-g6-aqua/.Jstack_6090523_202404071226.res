
==============================
  T0 -> 0.0-5.63
  T1 -> 0
  T2 -> 1-5
------------------------------
 -0- /home/export/online1/mdt00/shisuan/swgbcm/xduan/demo-g6-aqua/../grist-opt/bld/ParGRIST-AMIPW-sunway_xfort.exe 
   -1- slave__Waiting_For_Task ([0/384])
     -2- slave_jobserver at jobserver-slave.S:83 ([384/384] T0 0.63/v0 0.62/v0 0.61/v0...)
   -1- main at ../../src/atmosphere/gcm/grist_atmos.F90:24 ([0/6])
     -2- grist_atmos at ../../src/atmosphere/gcm/grist_atmos.F90:160 ([0/6])
       -3- grist_gcm_control_driver::grist_gcm_run at ../../src/atmosphere/gcm/grist_gcm_control_driver.F90:353 ([0/6])
         -4- grist_gcm_io_h1_module::gcm_output_history_h1 at ../../src/atmosphere/gcm/io_template/grist_gcm_io_h1_module.F90:412 ([0/6])
           -5- grist_gcm_io_h1_module::gcm_output_atm_h1_file at ../../src/atmosphere/gcm/io_template/grist_gcm_io_h1_module.F90:636 ([0/6])
             -6- grist_gcm_io_h1_module::write_atm_file1d at ../../src/atmosphere/gcm/io_template/grist_gcm_io_h1_module.F90:966 ([0/6])
               -7- grist_fileio_list_1d_module_par::wrap_output_1d_unordered_sp at ../../src/infrastructure/io/grist_fileio_list_1d_module_par.F90:708 ([0/1])
                 -8- grist_wrap_pf::wrap_create at ../../src/infrastructure/io/grist_wrap_pf.F90:35 ([0/1])
                   -9- nfmpi_create_ ([0/1])
                     -10- ncmpi_create ([0/1])
                       -11- ncmpiio_create ([0/1])
                         -12- PMPI_File_open ([1/1] T1 0/v0)
               -7- grist_fileio_list_1d_module_par::wrap_output_1d_unordered_sp at ../../src/infrastructure/io/grist_fileio_list_1d_module_par.F90:747 ([0/5])
                 -8- pmpi_comm_split__ ([0/5])
                   -9- PMPI_Comm_split ([5/5] T2 1/v0 2/v0 3/v0...)
==============================

==========================================================
node(taskid):svrstart,wait,sout
vn000000(0       ):	0.50	10.51	10.51
before:0.005228,scan:31.007893,first_data:30.007665,process:1.000289,show:0.021716,total:31.034898

