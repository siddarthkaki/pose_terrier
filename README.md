# seeker_pose_estimator

## Dependencies
* C++11 or later
* [CMake](https://cmake.org/) (3.10 or later)
* [Ceres Solver](http://ceres-solver.org/) (tested with 1.13.0)
* [Eigen](http://eigen.tuxfamily.org/) (tested with 3.3.4)
* [Protocol Buffers](https://developers.google.com/protocol-buffers/) (tested with 3.0.0)

## Build instructions

Clone the repository:
```
git clone https://github.com/autognc/seeker_pose_estimator
```

CMake is employed for building. From the cloned directory, execute:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```
