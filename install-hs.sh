cd packages/
tar -zxvf higgssignals.tar.gz
cd higgssignals
mkdir build && cd build 
cmake .. -DFeynHiggs_ROOT=/pMSSM_McMC/packages/FeynHiggs-2.18.0 && make
