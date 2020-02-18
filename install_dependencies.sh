#!/bin/sh
set -ex
mkdir Installs
cd Installs
export INSTALLATION_PATH=$PWD

sudo apt-get -qq update #-qq for less output
#sudo apt-get -y upgrade #do not run this on travis, takes too much time, only if necessary

# Helpful utils
sudo apt-get -y install wget

# For cloning
# sudo apt-get -y install git

# Development
sudo apt-get -y install build-essential cmake
sudo apt-get -y install pkg-config

# Boost
sudo apt-get -y install libboost-all-dev

# OpenCV dependencies
# File I/O
sudo apt-get -y install libjpeg8-dev libpng12-dev libtiff5-dev libjasper-dev
#Video I/O
sudo apt-get -y install libavcodec-dev libavformat-dev libswscale-dev libavresample-dev
sudo apt-get -y install libv4l-dev
sudo apt-get -y install libx264-dev
#sudo apt-get -y install ffmpeg
#sudo apt-get -y install libxvidcore-dev
# GUI
sudo apt-get -y install libgtk-3-dev
sudo apt-get -y install libx11-dev

# libxml
wget https://github.com/GNOME/libxml2/archive/v2.9.9.tar.gz
tar -xzf v2.9.9.tar.gz
cd libxml2-2.9.9
./autogen.sh
make
sudo make install
sudo ldconfig
cd $INSTALLATION_PATH

# Optimization
sudo apt-get -y install libatlas-base-dev gfortran liblapacke-dev
sudo apt-get -y install libpthread-stubs0-dev

# Python 3...?
#sudo apt-get -y install python3-dev
#sudo apt-get -y install python3-numpy

# OpenCV
#git clone https://github.com/opencv/opencv.git
wget https://github.com/opencv/opencv/archive/3.3.0.tar.gz
tar -zxf 3.3.0.tar.gz -C $INSTALLATION_PATH
cd opencv-3.3.0
mkdir build
cd build
CC=gcc \
CXX=g++ \
FC=gfortran \
CFLAGS="-O3 -fPIC" \
CPPFLAGS="-O3 -fPIC" \
CXXFLAGS="-O3 -fPIC" \
LDFLAGS="-Wl,-rpath, -ldl -lrt" \
cmake \
    -DCMAKE_BUILD_TYPE=RELEASE \
    -DWITH_CUDA=OFF \
    -DWITH_DOCS=OFF \
    -DBUILD_TESTS=ON \
    -DBUILD_EXAMPLES=OFF \
    -DINSTALL_C_EXAMPLES=OFF \
    -DBUILD_SHARED_LIBS=ON \
    ../
make -j 8
sudo make install
sudo ldconfig
cd $INSTALLATION_PATH

#dlib
wget https://github.com/davisking/dlib/archive/v19.7.tar.gz
tar -zxf v19.7.tar.gz -C $INSTALLATION_PATH
mkdir dlib_build
cd dlib_build
cmake \
    -DCMAKE_BUILD_TYPE=RELEASE \
    ../dlib-19.7
make -j 8
sudo make install
sudo ldconfig

