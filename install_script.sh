#!/usr/bin/env bash
set -e

echo "This will attempt to download all the required libraries (potentially outdated and likely redundant) and compile QTLtools. You MUST have patience and basic utilities like C/C++ and fortran compilers, tar, gzip, bzip, and wget. It will create two new directories downloads and install, and when it's done will consume ~1.6G of space. It has been tested with several linux distros (ubuntu, arch, centos, manjaro) but is NOT guaranteed to work."
echo -e "\n    Usage ./install_script.sh [NoOfThreads]\n"


while true; do
    read -p "Do you want to continue (y/n)? " yn
    case $yn in
        y ) break;;
        n ) exit;;
        * ) echo "Please answer y or n.";;
    esac
done

FOO=$1
THREADS=${FOO:=1}
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]
then
        echo "Give a positive integer for the number of threads!"
        exit 1
fi

exit_fail(){ 
    echo "Script failed please check the previous output for the reason" 
}

trap exit_fail EXIT


#CREATE DIR STRUCTURE
mkdir downloads
mkdir -p install/lib/
mkdir install/include/

LD_LIBRARY_PATH_BCK=$LD_LIBRARY_PATH
LIBRARY_PATH_BCK=$LIBRARY_PATH
CPATH_BCK=$CPATH
PATH_BCK=$PATH
PKG_CONFIG_PATH_BCK=$PKG_CONFIG_PATH

QTLTOOLS_DIR=$PWD

export LIBRARY_PATH=$PWD/install/lib/:$LIBRARY_PATH
export LD_LIBRARY_PATH=$PWD/install/lib/:$LD_LIBRARY_PATH
export CPATH=$PWD/install/include/:$CPATH
export PATH=$PWD/install/bin/:$PATH
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:/$PWD/install/lib/pkgconfig/

cd downloads

#ZLIB
wget --no-check-certificate https://zlib.net/zlib-1.2.11.tar.gz
tar zxvf zlib-1.2.11.tar.gz
rm zlib-1.2.11.tar.gz
cd zlib-1.2.11
./configure --prefix=$PWD/../../install/
make -j $THREADS
make install
cd ..

#BZIP2
wget --no-check-certificate https://sourceware.org/pub/bzip2/bzip2-1.0.8.tar.gz
tar zxvf bzip2-1.0.8.tar.gz
rm bzip2-1.0.8.tar.gz
cd bzip2-1.0.8
make install PREFIX=$PWD/../../install/ CFLAGS="-Wall -Winline -O2 -g -D_FILE_OFFSET_BITS=64 -fPIC" -j $THREADS
cd ..

#LZMA
wget --no-check-certificate https://tukaani.org/xz/xz-5.2.5.tar.gz
tar zxvf xz-5.2.5.tar.gz
rm xz-5.2.5.tar.gz
cd xz-5.2.5
./configure --prefix=$PWD/../../install/
make -j $THREADS
make install
cd ..

#PCRE2
wget --no-check-certificate https://github.com/PhilipHazel/pcre2/releases/download/pcre2-10.37/pcre2-10.37.tar.gz
tar zxvf pcre2-10.37.tar.gz
rm pcre2-10.37.tar.gz
cd pcre2-10.37
./configure --prefix=$PWD/../../install/ --enable-jit
make -j $THREADS
make install
cd ..

#CURL
wget --no-check-certificate https://curl.se/download/curl-7.79.1.tar.gz
tar zxvf curl-7.79.1.tar.gz
rm curl-7.79.1.tar.gz
cd curl-7.79.1
./configure --prefix=$PWD/../../install/ --with-openssl
make -j $THREADS
make install
cd ..


#GSL
wget --no-check-certificate https://ftp.gnu.org/gnu/gsl/gsl-2.7.tar.gz
tar zxvf gsl-2.7.tar.gz
rm gsl-2.7.tar.gz
cd gsl-2.7
./configure --prefix=$PWD/../../install/
make -j $THREADS
make install
cd ..

#HTSlib
wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
tar xvfj htslib-1.13.tar.bz2
rm htslib-1.13.tar.bz2
cd htslib-1.13
./configure --prefix=$PWD/../../install/ --disable-libcurl
make -j $THREADS
make install
cd ..

#BOOST
wget --no-check-certificate https://boostorg.jfrog.io/artifactory/main/release/1.77.0/source/boost_1_77_0.tar.bz2
tar xvfj boost_1_77_0.tar.bz2
rm boost_1_77_0.tar.bz2
cd boost_1_77_0
./bootstrap.sh  --with-libraries=iostreams,program_options --prefix=$PWD/../../install/
./b2 install
cd ..

#R
wget --no-check-certificate https://cran.r-project.org/src/base/R-4/R-4.1.1.tar.gz
tar zxvf R-4.1.1.tar.gz
rm R-4.1.1.tar.gz
cd R-4.1.1
./configure --prefix=$PWD/../../install/ --libdir=$PWD/../../install/lib/ --without-readline --without-x
cd src/nmath/standalone/
make -j $THREADS
make install
cd $QTLTOOLS_DIR

#QTLTOOLS
make -f Makefile4InstallScript -j $THREADS

#CLEAN UP
#rm -rf install/ download/

#RESTORE ENV VARS
export LIBRARY_PATH=$LIBRARY_PATH_BCK
export CPATH=$CPATH_BCK
export PATH=$PATH_BCK
export PKG_CONFIG_PATH=$PKG_CONFIG_PATH_BCK
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH_BCK

trap - EXIT
echo "Script succesfully completed. If you want to install QTLtools please type (as root): make install. You may want to delete the install and downloads directories as well."
