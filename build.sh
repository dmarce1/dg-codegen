mkdir $1
cd $1
cmake  -DCMAKE_INSTALL_PREFIX=$HOME/local/$1/ -DCMAKE_BUILD_TYPE=$1 ..
make -j12


