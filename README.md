# pose_terrier

## Dependencies
* C++11 or later
* [CMake](https://cmake.org/) (3.10 or later)
* [Ceres Solver](http://ceres-solver.org/) (tested with 1.13.0)
* [Eigen](http://eigen.tuxfamily.org/) (tested with 3.3.4)

## Build instructions

Clone the repository:
```
git clone [repository url]
```

CMake is employed for building. From the cloned directory, execute:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```
