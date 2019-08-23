cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 11)

include(FindProtobuf)
find_package(Protobuf REQUIRED)
find_package(OpenCV REQUIRED)

#add_executable(single_pose_estimator_test single_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
#target_link_libraries(single_pose_estimator_test ceres glog)

#add_executable(multiple_pose_estimator_test multiple_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
#target_link_libraries(multiple_pose_estimator_test ceres glog)

#add_executable(pnp_multiple_pose_estimator_test pnp_multiple_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
#target_link_libraries(pnp_multiple_pose_estimator_test ceres glog ${OpenCV_LIBS})

add_executable(pose_tracking_test pose_tracking_test.cc Utilities.cc PoseSolver.cc KalmanFilter.cc)
target_link_libraries(pose_tracking_test ceres glog)

add_executable(pnp_pose_tracking_test pnp_pose_tracking_test.cc Utilities.cc PoseSolver.cc KalmanFilter.cc)
target_link_libraries(pnp_pose_tracking_test ceres glog ${OpenCV_LIBS})

#add_executable(ekf_pose_tracking_test ekf_pose_tracking_test.cc Utilities.cc PoseSolver.cc KalmanFilter.cc)
#target_link_libraries(ekf_pose_tracking_test ceres glog)

add_executable(pose_estimator_main pose_estimator_main.cc Utilities.cc PoseSolver.cc KalmanFilter.cc pose.pb.cc measurement.pb.cc)
target_include_directories(pose_estimator_main PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pose_estimator_main ceres glog ${PROTOBUF_LIBRARY})

add_executable(pnp_pose_estimator_main pnp_pose_estimator_main.cc Utilities.cc PoseSolver.cc KalmanFilter.cc pose.pb.cc measurement.pb.cc)
target_include_directories(pnp_pose_estimator_main PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pnp_pose_estimator_main ceres glog ${PROTOBUF_LIBRARY} ${OpenCV_LIBS})

add_executable(simulate_measurements simulate_measurements.cc Utilities.cc pose.pb.cc measurement.pb.cc)
target_include_directories(simulate_measurements PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(simulate_measurements glog ${PROTOBUF_LIBRARY})

#add_executable(pose_reader_example pose_reader_example.cc pose.pb.cc)
#target_include_directories(pose_reader_example PUBLIC ${PROTOBUF_INCLUDE_DIR})
#target_link_libraries(pose_reader_example glog ${PROTOBUF_LIBRARY})

