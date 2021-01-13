#RAZOR
Collection of functions and drivers to perform computational reduction and/or feature identification. Approaches implemented: Reduced Basis, Manifold Learning.

Requirements
------------

To install this project, please ensure that you have installed the following (installation guides are provided on the respective websites or github repositories):

SOFTWARES:
  - [Git](http://git-scm.com)
  - A C++ compiler, e.g., [GCC](https://gcc.gnu.org/), [clang](http://clang.llvm.org/) [Verified on GCC v6.4.0 or newer]
  - [CMake](http://www.cmake.org)

LIBRARIES:
  - [smart-math](https://github.com/strath-ace-labs/smart-math)
  - [smart-uq](https://github.com/strath-ace-labs/smart-uq)
  - [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), v3.3.5 or newer
	- Optional: [pagmo](https://esa.github.io/pagmo2/install.html) 
	- Optional: [GSL](https://www.gnu.org/software/gsl/)

Edit CMAkeLists.txt file to switch ON-OFF the optional libraries

Compilation
------------

Run the following commands to download, build, and compile this project.

    cd path_to_compilation/
    git clone https://github.com/strath-ace-labs/RAZOR.git
    (Enter your username and password)
    cd RAZOR
    mkdir build && cd build
    cmake .. 
    make


Build options
-------------

To be able to use optional capabilities you can change the options in the CMakelists.txt in RAZOR main folder before running cmake or run `ccmake ..` to bring up the interface that can be used to toggle options.


Running RAZOR
-------------
For instruction on how to run the code, please refer to the tutorial folder


To do
------
- Implement Artifical Neural Network approach
- Regression tests 
- Adding new tutorials

