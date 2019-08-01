cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 11)

#add_executable(single_pose_estimator_test single_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
#target_link_libraries(single_pose_estimator_test ceres glog)

#add_executable(multiple_pose_estimator_test multiple_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
#target_link_libraries(multiple_pose_estimator_test ceres glog)

add_executable(pose_tracking_test pose_tracking_test.cc Utilities.cc PoseSolver.cc KalmanFilter.cc)
target_link_libraries(pose_tracking_test ceres glog)

add_executable(ekf_pose_tracking_test ekf_pose_tracking_test.cc Utilities.cc PoseSolver.cc KalmanFilter.cc)
target_link_libraries(ekf_pose_tracking_test ceres glog)