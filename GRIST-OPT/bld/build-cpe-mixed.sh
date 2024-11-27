
export PATH=$PWD/../../SWGOMP/script/:${PATH}
cd .. 
rm src 
ln -sf src_mixed_new src

cd bld

make GCM_AMIPW_MIXED-mixed_sunway_xfort
