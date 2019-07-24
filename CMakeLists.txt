set (CMAKE_CXX_STANDARD 11)

add_executable(main main.cc Utilities.cc PoseSolver.cc)
target_link_libraries(main ceres glog)

add_executable(multiple_pose_estimator_test multiple_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
target_link_libraries(multiple_pose_estimator_test ceres glog)