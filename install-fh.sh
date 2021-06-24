cd packages 
tar -zxvf FeynHiggs-2.18.0.tar.gz
cd FeynHiggs-2.18.0
export CMAKE_CURRENT_SOURCE_DIR=$PWD
make && make install
