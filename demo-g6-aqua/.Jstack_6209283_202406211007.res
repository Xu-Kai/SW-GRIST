
==============================
  T0 -> 0.0-127.63
  T1 -> 0-127
------------------------------
 -0- /home/export/online1/mdt00/shisuan/swgbcm/xukai/grist_bencmark/demo-g6-aqua/../grist-opt/bld/ParGRIST-AMIPW_MIXED-opt_mixed_sunway_xfort.exe 
   -1- slave__Waiting_For_Task ([8192/8192] T0 6.63/v49993 6.62/v49993 6.61/v49993...)
   -1- main at ../../src_mixed_new/atmosphere/gcm/grist_atmos.F90:24 ([0/128])
     -2- grist_atmos at ../../src_mixed_new/atmosphere/gcm/grist_atmos.F90:68 ([0/128])
       -3- cesm_swlu_debug_init_ at ../../src_mixed_new/infrastructure/utils/swlu_prof_wrapper.c:49 ([0/128])
         -4- swlu_debug_init at samprof.c:1665 ([0/128])
           -5- swlu_read_elf at samprof.c:639 ([0/128])
             -6- sections_init at samprof.c:390 ([0/128])
               -7- strcmp at ../sysdeps/sw_64/strcmp.S:46 ([0/128])
                 -8- signal handler called ([0/128])
                   -9- swch_catchsig ([128/128] T1 6/v49993 7/v49993 8/v49993...)
==============================

==========================================================
node(taskid):svrstart,wait,sout
vn049992(0       ):	0.55	10.56	10.56
vn049993(6       ):	0.54	10.55	10.57
vn050058(12      ):	0.54	10.55	10.55
vn050059(18      ):	0.57	10.58	10.58
vn050060(24      ):	0.54	10.55	10.56
vn050061(30      ):	0.53	10.54	10.54
vn050064(36      ):	0.52	10.53	10.53
vn050065(42      ):	0.59	10.59	10.59
vn050066(48      ):	0.55	10.56	10.56
vn050067(54      ):	0.53	10.53	10.55
vn050068(60      ):	0.58	10.59	10.59
vn050069(66      ):	0.60	10.64	10.65
vn050070(72      ):	0.52	10.52	10.53
vn050071(78      ):	0.55	10.56	10.58
vn050072(84      ):	0.55	10.55	10.55
vn050073(90      ):	0.61	10.62	10.62
vn050074(96      ):	0.58	10.59	10.59
vn050075(102     ):	0.60	10.60	10.61
vn050076(108     ):	0.53	10.53	10.53
vn050077(114     ):	0.55	10.56	10.57
vn050078(120     ):	0.57	10.57	10.58
vn050079(126     ):	0.47	10.48	10.48
before:0.004934,scan:31.008000,first_data:30.006767,process:1.001294,show:0.027947,total:31.040942

