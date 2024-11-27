export PATH=$PWD/../../SWGOMP/script/:${PATH}
cd .. 
rm src
ln -sf src-main src
cd bld
make GCM_AMIPW-sunway_xfort
