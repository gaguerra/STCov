# STCov

STCov is a program for calculating statistics on pairwise coalescence times given a rooted, bifurcating species tree

##Installation instructions

STCov is implemented in C++ (a precompiled binary is not yet added to this github, but will be added soon). 

If you choose to build the software from scratch, follow the build instructions below

##Build Instructions

STCov requires the following libraries and executables in order to compile and run: 

* A C++-11 compiler (gcc 4.8 or later, for example)
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) - Used for linear algebra and matrix calculations

On OSX, from inside the STCov folder, build STCov by running: 

### clang++ -w -std=c++11 -stdlib=libc++ -I /usr/local/include/ main.cpp -o STCov


