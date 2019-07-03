add_executable(pose_estimator pose_estimator.cc Utilities.cc PoseSolver.cc)
target_link_libraries(pose_estimator ceres glog)
