## Debian/Ubuntu
(enable universe repository)

sudo apt-get install git g++ autoconf libtool zlib1g zlib1g-dev libnifti2 libnifti-dev libgsl0-dev libxerces-c-dev libglib2.0-dev
git clone https://github.com/janeisklar/RSTools.git rstools
cd rstools
./autogen.sh
./configure CFLAGS=-I/usr/include/nifti CPPFLAGS=-I/usr/include/nifti
make
sudo make install

export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH

## CentOS

sudo yum install cmake git gcc-c++ autoconf libtool zlib zlib-devel gsl gsl-devel xerces-c xerces-c-devel glib2 glib2-devel

wget http://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz\?r\= -O nifticlib.tgz
tar xzf nifticlib.tgz
cd nifticlib-2.0.0
mkdir build
cd build
cmake -D BUILD_SHARED_LIBS=On ..
make
sudo make install

git clone https://github.com/janeisklar/RSTools.git rstools
cd rstools
./autogen.sh
./configure CFLAGS=-I/usr/local/include/nifti CPPFLAGS=-I/usr/local/include/nifti
make
sudo make install