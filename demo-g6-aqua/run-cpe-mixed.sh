

cp grist_amipw_phys.nml.rad grist_amipw_phys.nml 
bsub  -I -q q_share -J grist-cpe-mixed -n 128 -share_size 12000 -cgsp 64 -host_stack 1024 -cache_size 128 -b -o grist-cpe-mixed.log ../GRIST-OPT/bld/ParGRIST-GCM_AMIPW_MIXED-mixed_sunway_xfort.exe
