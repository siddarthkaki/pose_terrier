cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 11)

find_package(PythonLibs 2.7 REQUIRED)

#add_executable(single_pose_estimator_test single_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
#target_link_libraries(single_pose_estimator_test ceres glog)

#add_executable(multiple_pose_estimator_test multiple_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
#target_link_libraries(multiple_pose_estimator_test ceres glog)

add_executable(pose_tracking_test pose_tracking_test.cc Utilities.cc PoseSolver.cc third_party/KalmanFilter/KalmanFilter.cc)
target_include_directories(pose_tracking_test PUBLIC third_party/KalmanFilter ${PYTHON_INCLUDE_DIRS})
target_link_libraries(pose_tracking_test ceres glog ${PYTHON_LIBRARIES})