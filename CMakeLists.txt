project(pose_terrier)
cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 14)
set(Boost_USE_STATIC_LIBS OFF)      # uncomment when performing Euler rotational dynamics propagation in simulate_measurements
set(Boost_USE_MULTITHREADED ON)     # uncomment when performing Euler rotational dynamics propagation in simulate_measurements
set(Boost_USE_STATIC_RUNTIME OFF)   # uncomment when performing Euler rotational dynamics propagation in simulate_measurements

include(FindProtobuf)
find_package(Protobuf REQUIRED)
find_package(Ceres REQUIRED)
find_package(glog REQUIRED)
find_package(OpenCV REQUIRED)   # uncomment when employing PnP
find_package(Boost 1.74.0)      # uncomment when performing Euler rotational dynamics propagation in simulate_measurements

SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# the primary real-time executable with both EPnP and NLS for static pose estimation
add_executable(pose_estimator_main pose_estimator_main.cc Utilities.cc PoseSolver.cc MEKF2.cc pose.pb.cc measurement.pb.cc)
target_include_directories(pose_estimator_main PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pose_estimator_main ceres glog::glog ${PROTOBUF_LIBRARY} ${OpenCV_LIBS} pthread)

# the primary real-time executable with only NLS for static pose estimation
add_executable(nls_pose_estimator_main nls_pose_estimator_main.cc Utilities.cc PoseSolver.cc MEKF2.cc pose.pb.cc measurement.pb.cc)
target_include_directories(nls_pose_estimator_main PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(nls_pose_estimator_main ceres glog::glog ${PROTOBUF_LIBRARY} pthread)

# the primary real-time executable with the quatera module built-in
add_executable(pose_estimator_main_quatera pose_estimator_main_quatera.cc Utilities.cc PoseSolver.cc MEKF2.cc QuateRA.cc pose.pb.cc measurement.pb.cc)
target_include_directories(pose_estimator_main_quatera PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pose_estimator_main_quatera ceres glog::glog ${PROTOBUF_LIBRARY} pthread)

# real-time executable with only pose as measurement inputs; deprecated
#add_executable(pose_estimator_pose_only pose_estimator_pose_only.cc Utilities.cc PoseSolver.cc MEKF2.cc pose.pb.cc measurement.pb.cc)
#target_include_directories(pose_estimator_pose_only PUBLIC ${PROTOBUF_INCLUDE_DIR})
#target_link_libraries(pose_estimator_pose_only ceres glog::glog ${PROTOBUF_LIBRARY} pthread)

# the primary real-time executable with relative bearings to keypoints and also pose as measurement inputs
#add_executable(pose_estimator_x pose_estimator_x.cc Utilities.cc PoseSolver.cc MEKF2.cc pose.pb.cc measurementx.pb.cc)
#target_include_directories(pose_estimator_x PUBLIC ${PROTOBUF_INCLUDE_DIR})
#target_link_libraries(pose_estimator_x ceres glog::glog ${PROTOBUF_LIBRARY} pthread)

# the primary real-time executable with OpenCV's PnP instead of NLS for static pose estimation; deprecated
#add_executable(pnp_pose_estimator_main pnp_pose_estimator_main.cc Utilities.cc MEKF2.cc pose.pb.cc measurement.pb.cc)
#target_include_directories(pnp_pose_estimator_main PUBLIC ${PROTOBUF_INCLUDE_DIR})
#target_link_libraries(pnp_pose_estimator_main ceres glog::glog ${PROTOBUF_LIBRARY} ${OpenCV_LIBS})

# test executable for multiple random static pose estimates
add_executable(multiple_pose_estimator_test multiple_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
target_link_libraries(multiple_pose_estimator_test ceres glog::glog pthread ${OpenCV_LIBS})

# test executable for multiple random static pose estimates with only NLS
add_executable(nls_multiple_pose_estimator_test nls_multiple_pose_estimator_test.cc Utilities.cc PoseSolver.cc)
target_link_libraries(nls_multiple_pose_estimator_test ceres glog::glog pthread)

# test executable for multiple random static pose estimates with OpenCV's PnP instead of NLS
add_executable(pnp_multiple_pose_estimator_test pnp_multiple_pose_estimator_test.cc Utilities.cc)
target_link_libraries(pnp_multiple_pose_estimator_test ceres ${OpenCV_LIBS})

# test executable to simulate keypoint relative bearing estimates
add_executable(simulate_measurements simulate_measurements.cc Utilities.cc pose.pb.cc measurement.pb.cc)
target_include_directories(simulate_measurements PUBLIC ${PROTOBUF_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
target_link_libraries(simulate_measurements glog::glog ${PROTOBUF_LIBRARY} ${Boost_LIBRARIES})
set_target_properties(simulate_measurements PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

# test executable to read and re-publish previously-saved keypoint relative bearing estimates
add_executable(read_and_publish_measurements read_and_publish_measurements.cc Utilities.cc pose.pb.cc measurement.pb.cc)
target_include_directories(read_and_publish_measurements PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(read_and_publish_measurements glog::glog ${PROTOBUF_LIBRARY})
set_target_properties(read_and_publish_measurements PROPERTIES LINK_FLAGS "-Wl,--no-as-needed")

# test executable to demonstrate reading ProtoBuf-packaged poses from a named pipe
add_executable(pose_reader_example pose_reader_example.cc pose.pb.cc)
target_include_directories(pose_reader_example PUBLIC ${PROTOBUF_INCLUDE_DIR})
target_link_libraries(pose_reader_example glog::glog ${PROTOBUF_LIBRARY})