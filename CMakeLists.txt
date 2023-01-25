project(pose_terrier)
cmake_minimum_required(VERSION 3.10)

set (CMAKE_CXX_STANDARD 14)
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

include(FindProtobuf)
find_package(Protobuf REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Ceres REQUIRED)
find_package(glog REQUIRED)
find_package(Boost 1.74.0) 

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

#add_executable(single_pose_estimator_test single_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
#target_link_libraries(single_pose_estimator_test ceres glog)

add_executable(multiple_pose_estimator_test multiple_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
target_link_libraries(multiple_pose_estimator_test ceres glog::glog pthread)

add_executable(pnp_multiple_pose_estimator_test pnp_multiple_pose_estimator_test.cc Utilities.cc)
target_link_libraries(pnp_multiple_pose_estimator_test ceres ${OpenCV_LIBS})

add_executable(pose_tracking_test pose_tracking_test.cc Utilities.cc PoseSolver.cc KalmanFilter.cc)
target_link_libraries(pose_tracking_test ceres glog::glog)

add_executable(pnp_pose_tracking_test pnp_pose_tracking_test.cc Utilities.cc KalmanFilter.cc)
target_link_libraries(pnp_pose_tracking_test ceres ${OpenCV_LIBS})

add_executable(mekf_pose_tracking_test mekf_pose_tracking_test.cc Utilities.cc PoseSolver.cc KalmanFilter.cc MEKF.cc MEKF2.cc)
target_link_libraries(mekf_pose_tracking_test ceres glog::glog)

add_executable(pnp_mekf_pose_tracking_test pnp_mekf_pose_tracking_test.cc Utilities.cc KalmanFilter.cc MEKF.cc MEKF2.cc)
target_link_libraries(pnp_mekf_pose_tracking_test ceres glog::glog ${OpenCV_LIBS})

add_executable(pose_estimator_main pose_estimator_main.cc Utilities.cc PoseSolver.cc KalmanFilter.cc MEKF.cc MEKF2.cc pose.pb.cc measurement.pb.cc)
target_include_directories(pose_estimator_main PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pose_estimator_main ceres glog::glog ${PROTOBUF_LIBRARY} pthread)

add_executable(pose_estimator_main_with_quatera pose_estimator_main_with_quatera.cc Utilities.cc PoseSolver.cc KalmanFilter.cc MEKF.cc MEKF2.cc QuateRA.cc pose.pb.cc measurement.pb.cc)
target_include_directories(pose_estimator_main_with_quatera PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pose_estimator_main_with_quatera ceres glog::glog ${PROTOBUF_LIBRARY} pthread)

add_executable(pose_estimator_pose_only pose_estimator_pose_only.cc Utilities.cc PoseSolver.cc KalmanFilter.cc MEKF.cc MEKF2.cc pose.pb.cc measurement.pb.cc)
target_include_directories(pose_estimator_pose_only PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pose_estimator_pose_only ceres glog::glog ${PROTOBUF_LIBRARY} pthread)

add_executable(pose_estimator_x pose_estimator_x.cc Utilities.cc PoseSolver.cc KalmanFilter.cc MEKF.cc MEKF2.cc pose.pb.cc measurementx.pb.cc)
target_include_directories(pose_estimator_x PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pose_estimator_x ceres glog::glog ${PROTOBUF_LIBRARY} pthread)

add_executable(pnp_pose_estimator_main pnp_pose_estimator_main.cc Utilities.cc KalmanFilter.cc MEKF.cc MEKF2.cc pose.pb.cc measurement.pb.cc)
target_include_directories(pnp_pose_estimator_main PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pnp_pose_estimator_main ceres glog::glog ${PROTOBUF_LIBRARY} ${OpenCV_LIBS})

add_executable(simulate_measurements simulate_measurements.cc Utilities.cc pose.pb.cc measurement.pb.cc)
target_include_directories(simulate_measurements PUBLIC ${PROTOBUF_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
target_link_libraries(simulate_measurements glog::glog ${PROTOBUF_LIBRARY} ${Boost_LIBRARIES})
set_target_properties(simulate_measurements PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

add_executable(read_and_publish_measurements read_and_publish_measurements.cc Utilities.cc pose.pb.cc measurement.pb.cc)
target_include_directories(read_and_publish_measurements PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(read_and_publish_measurements glog::glog ${PROTOBUF_LIBRARY})
set_target_properties(read_and_publish_measurements PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

add_executable(pose_reader_example pose_reader_example.cc pose.pb.cc)
target_include_directories(pose_reader_example PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pose_reader_example glog::glog ${PROTOBUF_LIBRARY})
