

cp grist_amipw_phys.nml.org grist_amipw_phys.nml 
bsub  -I -q q_share -J grist-mpe-dp -n 128 -share_size 12000 -cgsp 64 -host_stack 1024 -cache_size 128 -b -o grist-mpe-dp.log ../GRIST/bld/ParGRIST-GCM_AMIPW-sunway_xfort.exe
