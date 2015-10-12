# Installation instructions
## General
These instructions are intended to be executed following the OS-specific installation steps that are outlined in the subsequent chapters.
```bash
# Download job files
cd <anywhere you like>
git clone git@github.com:janeisklar/RSTools-preprocessing-jobs.git jobs
# if you are _NOT_ using the CMRR multiband EPI sequence:
git checkout tags/v2

# Download additional resources (templates, masks, etc.)
git clone git@github.com:janeisklar/RSTools-preprocessing-resources.git resources
```

To configure paths to additional tools and resources open /usr/local/etc/rstools.conf and adjust the paths for: fslPath, ANTSPATH, afniPath and ppsPath. The latter one (ppsPath) should point to the 'resources' directory you've downloaded in the previous step.

Also see the usage instructions at the end of the document.


## Debian/Ubuntu
(enable universe repository)

```bash
# compile rstools core package (contains resting-state tools only)
sudo apt-get install git g++ autoconf libtool zlib1g zlib1g-dev libnifti2 libnifti-dev libgsl0-dev libxerces-c-dev libglib2.0-dev

git clone https://github.com/janeisklar/RSTools.git
cd RSTools
./autogen.sh
./configure CFLAGS=-I/usr/include/nifti CPPFLAGS=-I/usr/include/nifti
make
sudo make install

echo 'export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH' >>~/.bashrc 

# compile rstools for the preprocessing of fMRI data
cd ..
git clone https://github.com/janeisklar/RSTools-preprocessing.git
cd RSTools-preprocessing
./autogen.sh
./configure
make
sudo make install

# compile rstools jobeditor (for editing/creating pipelines for batch preprocessing)
sudo apt-get install qt5-default qttools5-dev-tools

cd ..
git clone https://github.com/janeisklar/RSTools-jobeditor.git
cd RSTools-jobeditor
./autogen.sh
./configure
make
sudo make install

```

## CentOS

```bash
# compile rstools core package (contains resting-state tools only)
sudo yum install cmake git gcc-c++ autoconf libtool zlib zlib-devel gsl gsl-devel xerces-c xerces-c-devel glib2 glib2-devel

wget http://downloads.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz\?r\= -O nifticlib.tgz
tar xzf nifticlib.tgz
cd nifticlib-2.0.0
mkdir build
cd build
cmake -D BUILD_SHARED_LIBS=On ..
make
sudo make install

git clone https://github.com/janeisklar/RSTools.git
cd RSTools
./autogen.sh
./configure CFLAGS=-I/usr/local/include/nifti CPPFLAGS=-I/usr/local/include/nifti
make
sudo make install

# compile rstools for the preprocessing of fMRI data
cd ..
git clone https://github.com/janeisklar/RSTools-preprocessing.git
cd RSTools-preprocessing
./autogen.sh
./configure
make
sudo make install

# compile rstools jobeditor (for editing/creating pipelines for batch preprocessing)
sudo yum groupinstall "C Development Tools and Libraries"
sudo rpm -Uvh http://dl.fedoraproject.org/pub/epel/6/x86_64/epel-release-6-8.noarch.rpm
sudo yum install mesa-libGL-devel qt5-qtbase qt5-qtbase-devel

cd ..
git clone https://github.com/janeisklar/RSTools-jobeditor.git
cd RSTools-jobeditor
./autogen.sh
./configure
make
sudo make install
```

## OS X (Yosemite)

```bash
# install homebrew if necesary

brew update
brew install xerces-c
brew install gsl
brew install gcc --without-multilib
brew install libtool
brew install automake
brew install autoconf
brew install pkg-config
brew install glib
brew install wget

# nifticlib
wget http://vorboss.dl.sourceforge.net/project/niftilib/nifticlib/nifticlib_2_0_0/nifticlib-2.0.0.tar.gz
tar xzf nifticlib-2.0.0.tar.gz
cd nifticlib-2.0.0
mkdir build
cd build
cmake -D BUILD_SHARED_LIBS=1 -D CMAKE_INSTALL_PREFIX=/usr/local/nifticlib ..
make

cd ../..

#rstools core
git clone https://github.com/janeisklar/RSTools.git
cd RSTools
./autogen.sh

export COMMONFLAGS="-I. -I/usr/local/nifticlib/include/nifti -I/usr/local/Cellar/gcc/5.2.0/include/c++/5.2.0 -I/usr/local/Cellar/gcc/5.2.0/include/c++/5.2.0/x86_64-apple-darwin14.5.0 -I/usr/local/include -D_DARWIN_C_SOURCE"
export CFLAGS=$COMMONFLAGS
export CPPFLAGS=$COMMONFLAGS
export CXXFLAGS=$CPPFLAGS 
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig 
export LDFLAGS="-L/usr/local/lib -L/usr/local/nifticlib/lib" 
export LIBRARY_PATH="/usr/local/nifticlib/lib;/usr/local/lib" 
export CC="/usr/local/bin/gcc-5 -nostdinc++" 
export CPP="/usr/local/bin/cpp-5 -nostdinc++"
export CXX="/usr/local/bin/g++-5 -nostdinc++"
./configure
make
make install
cd ..

ln -v -s /usr/local/nifticlib/lib/libznz.2.dylib /usr/local/lib/libznz.2.dylib
ln -v -s /usr/local/nifticlib/lib/libniftiio.2.dylib /usr/local/lib/libniftiio.2.dylib

# rstools preprocessing
git clone https://github.com/janeisklar/RSTools-preprocessing.git
cd RSTools-preprocessing
./autogen.sh
./configure
make
make install
cd ..

# rstools jobeditor
brew install qt
git clone https://github.com/janeisklar/RSTools-jobeditor.git

####################################################################################
# if you run into problems change the following:
#   QMAKE_CFLAGS_X86_64 += -Xarch_ppc64 -mmacosx-version-min=10.5
# to
#   QMAKE_CFLAGS_X86_64 += -mmacosx-version-min=10.5
# in: /usr/local/Cellar/qt/4.8.7/mkspecs/common/g++-macx.conf
####################################################################################

cd RSTools-jobeditor
./autogen.sh
./configure
make
make install

```

# Usage instructions
## General
All binaries resinde in /usr/local/bin if not specified differently when compiling and are prefixed with rs*. Each binary provides a documentation when run without any arguments (with the exception of rsjobeditor).

The job files available for download assume the following file structure which however is not fixed in any way and can be changed by modifying the parameters in the job file:
```
study
 ├── jobs
 │   └── 7t_rs_multiband_distortion_correction.job
 ├── preproc
 │   ├── m1
 │   │   ├── hc
 │   │   │   ├── my_subject_id_1
 │   │   │   │   ├── dvols.nii
 │   │   │   │   ├── advols.nii
 │   │   │   │   ├── uadvols.nii
 │   │   │   │   ├── buadvols.nii
 │   │   │   │   ├── rbuadvols.nii
 │   │   │   │   ├── wrbuadvols.nii
 │   │   │   │   ├── swrbuadvols.nii
 │   │   │   │   └── ...
 │   │   │   └── my_subject_id_2
 │   │   └── patient
 │   └── m2
 └── subjects
     ├── m1
     │   └── hc
     │       ├── my_subject_id_1
     │       │   ├── rfieldmap_LPI.nii
     │       │   └── vols.nii
     │       └── my_subject_id_2
     └── m2
```

The subjects folder hereby contains the unprocessed nifti files for every subject and will not be written to. All created files go into the preproc folder which can thus be easily deleted when changing parameters and starting over.

## Creating/editing a batch job
```bash
rsjobeditor </path/to/jobfile.job>
```

## Running a batch job
```bash
rsbatch -j 7t_rs_multiband_distortion_correction.job -A subject=my_subject_id_1 -A group=m1/patient -t 8
```
Any parameter that is used in the job file but was not specified needs to be specified using the -A command. In the same vein, parameters that were specified can be overwritten. The "-o" option is usfull to see what the pipeline would execute without actually running the job. For other options like skipping steps see the in-program help.

## Further information
Additional information about the history of a pre-processed file can be retrieved using rsinfo. It will also list dicom header information if a dicom directory for the scan was specified in the job file.
```bash
rsinfo -i swrbuadvols.nii
````

## Parallelization
The pre-processing of multiple subjects can be parallelized by using GNU parallels.
For example:
```bash
CMD="rsbatch -A subject={} -j 7t_rs_multiband_distortion_correction.job --threads=8 --quiet"
ls subjects/m1/hc | parallel -j+0 --eta -S10/: $CMD
```

## Known issues
Currently the folders in the 'preproc' directory will not be created automatically and need to exist beforehand. This will be fixed soon.