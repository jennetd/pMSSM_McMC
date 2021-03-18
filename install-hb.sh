cd packages 
tar -zxvf higgsbounds.tar.gz
cd higgsbounds
export CMAKE_CURRENT_SOURCE_DIR=$PWD
mkdir build && cd build 
cmake .. -DFeynHiggs_ROOT=/pMSSM_McMC/packages/FeynHiggs-2.16.1 && make
