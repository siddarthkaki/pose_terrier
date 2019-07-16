set (CMAKE_CXX_STANDARD 11)

add_executable(main main.cc Utilities.cc PoseSolver.cc)
target_link_libraries(main ceres glog)

add_executable(pose_estimator pose_estimator.cc Utilities.cc PoseSolver.cc)
target_link_libraries(pose_estimator ceres glog)