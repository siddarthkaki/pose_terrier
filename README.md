# pose_terrier

<img src="https://github.com/siddarthkaki/pose_terrier/blob/master/logo/pose_terrier_logo.png" width="300" height="300">

## Dependencies
* C++11 or later
* [CMake](https://cmake.org/) (3.10 or later)
* [Ceres Solver](http://ceres-solver.org/) (tested with 1.13.0)
* [Eigen](http://eigen.tuxfamily.org/) (tested with 3.3.4)
* [Protocol Buffers](https://developers.google.com/protocol-buffers/) (tested with 3.6.1.3)

## Build instructions
Clone the repository:
```
git clone [repository url]
```

Link Eigen libraries to proper directory (if needed). For example:
```
sudo ln -s /usr/include/eigen3/Eigen /usr/include/Eigen
sudo ln -s /usr/include/eigen3/unsupported/ /usr/include/unsupported
```

Compile protobufs for employed protobuf version. From cloned directory, execute:
```
protoc pose.proto --cpp_out=$PWD
protoc measurement.proto --cpp_out=$PWD
```

CMake is employed for building. From the cloned directory, execute:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make -j4
```

## Usage
Before the first usage, link the `params.json` file to the build directory. Execute from the `build` directory:
```
ln -s ../params.json params.json
```

To run the non-linear least-squares pose estimator + Kalman filter, execute from the `build` directory:
```
./pose_estimator_main
```
This executable expects measurements (packaged as a protocol buffer of the form specified in `measurement.proto`) over a named pipe (the path for which is specified by the `pipe_path_input` field of `params.json`) at each time-step.

As an example, measurements can be simulated and published over the input named pipe by executing the following from the `build` directory:
```
./simulate_measurements
```


The `pose_estimator_main` executable has two forms of outputs: 1) publishing the solved-for pose (packaged as a protocol buffer of the form specified in `pose.proto`) over a named pipe (the path for which is specified by the `pipe_path_output` field of `params.json`) at each time-step, and 2) logging data to csv files on disk. The options for enabling and configuring both forms of outputs are located in `params.json`.

As an example, poses published by `pose_estimator_main` over the output named pipe can be read and outputted to console by executing the following from the `build` directory:
```
./pose_reader_example
```
To analyse data logged to csv files on disk, some MATLAB scripts are provided. From the `scripts` directory, run `post_analysis_main.m` when true pose data is unavailable, or run `post_analysis_main_with_truth.m` when true pose data is available and is logged to csv file.
