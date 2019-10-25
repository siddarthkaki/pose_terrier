/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <Eigen/Core>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <fcntl.h>
#include <sys/stat.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>

#include "glog/logging.h"

#include "Utilities.h"
#include "KalmanFilter.h"
#include "MEKF.h"
#include "pose.pb.h"
#include "measurement.pb.h"

#include "third_party/json.hpp"
#include "third_party/CppRot/cpprot.h"

using Eigen::AngleAxisd;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXd;
using nlohmann::json;

#define GET_VARIABLE_NAME(Variable) (#Variable)

sig_atomic_t volatile finished = 0;

/**
 * @function exit_handler
 * @brief exit_handler function
 */
void exit_handler(int signal)
{
    printf("Caught signal %d\n", signal);
    finished = 1;
}

bool CheckValidMeasurement(const ProtoMeas::Measurements& measurements);

/**
 * @function main
 * @brief main function
 */
int main(int argc, char **argv)
{
    //google::InitGoogleLogging(argv[0]);
    google::InstallFailureSignalHandler();

    //-- Read-in problem geometry and params ---------------------------------/

    // read params from JSON file
    std::ifstream input_stream("params.json");
    json json_params;
    try
    {
        input_stream >> json_params;
    }
    catch (json::parse_error &e)
    {
        // output exception information
        std::cout << "message: " << e.what() << '\n'
                  << "exception id: " << e.id << '\n'
                  << "byte position of error: " << e.byte << std::endl;
    }

    const std::string pipe_path_input = json_params["pipe_path_input"];
    const std::string pipe_path_output = json_params["pipe_path_output"];

    // specify rigid position vector of camera wrt chaser in chaser frame
    Vector3d rCamVec;
    for (unsigned int idx = 0; idx < 3; idx++)
    {
        rCamVec(idx) = json_params["rCamVec"].at(idx);
    }

    // specify camera focal length
    double focal_length = json_params["focal_length"]; //5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    //double meas_std = double(json_params["meas_std_deg"]) * Utilities::DEG2RAD;

    // specify expected number of time-steps for memory pre-allocation
    const unsigned int num_poses_test = json_params["num_poses_test"];

    const bool output_to_pipe = json_params["output_to_pipe"];
    const bool log_to_file = json_params["log_to_file"];
    const bool log_periodically = json_params["log_periodically"];
    const unsigned int vector_reserve_size = json_params["vector_reserve_size"];

    const double kf_dt = json_params["kf_dt"];
    const double mekf_dt = kf_dt;

    //------------------------------------------------------------------------/

    //-- Init Filters --------------------------------------------------------/
    
    double kf_process_noise_std = 0.01;
    double kf_measurement_noise_std = 0.05;

    double mekf_process_noise_std = 0.01;
    double mekf_measurement_noise_std = 0.05;

    KF::KalmanFilter kf;
    kf.InitLinearPositionTracking(kf_process_noise_std, kf_measurement_noise_std, kf_dt);

    MEKF::MEKF mekf(mekf_dt);
    mekf.Init(mekf_process_noise_std, mekf_measurement_noise_std, mekf_dt);

    //-- Init sequence -------------------------------------------------------/

    // log path name prefixing and postfixing
    std::string init_time_str = std::to_string(std::time(nullptr));
    std::string prefix = "../data/" + init_time_str + "_";
    std::string postfix = ".csv";

    // declare vectors for storage
    std::vector<Pose> solved_poses, filtered_poses;
    //std::vector<VectorXd> kf_states;
    //std::vector<MatrixXd> kf_covars;
    std::vector<double> timestamps; // [s]

    // pre-allocate memory
    solved_poses.reserve(vector_reserve_size);
    filtered_poses.reserve(vector_reserve_size);
    //kf_states.reserve(vector_reserve_size);
    //kf_covars.reserve(vector_reserve_size);
    timestamps.reserve(vector_reserve_size);

    // clock object
    auto init_t = std::chrono::high_resolution_clock::now();
    double curr_elapsed_t = 0.0;

    // TODO: read in initial guess from json or other program
    // initial pose guess
    Pose pose0;
    pose0.pos << 0.0, 0.0, 25.0;
    pose0.quat.w() = 1.0;
    pose0.quat.vec() = Vector3d::Zero();

    std::cout << "Waiting for first measurement..." << std::endl;

    // set pipe names
    int fd_in, fd_out, rd_in = 0;
    const char *fifo_path_input = pipe_path_input.c_str();
    const char *fifo_path_output = pipe_path_output.c_str();

    // if pipes do not exist, create them
    struct stat buf;
    if (stat(fifo_path_input, &buf) != 0)
    {
        mkfifo(fifo_path_input, 0666); 
    }
    if (stat(fifo_path_output, &buf) != 0)
    {
        mkfifo(fifo_path_output, 0666); 
    }

    // set-up camera instrinic params
    MatrixXd PMat(3,3);
    PMat << focal_length,            0, 0,
                       0, focal_length, 0,
                       0,            0, 1;
    cv::Mat camera_matrix;
    cv::eigen2cv(PMat, camera_matrix);
    cv::Mat dist_coeffs = cv::Mat::zeros(4, 1, cv::DataType<double>::type); // Assuming no lens distortion

    // open FIFO for read only, with blocking to receive first measurement
    fd_in = open(fifo_path_input, O_RDONLY);

    // loop for waiting for first measurement to initialise filter
    bool received_first_meas = false;
    while (!received_first_meas)
    {
        // read byte size of measurement object from pipe
        size_t size;
        rd_in = read(fd_in, &size, sizeof(size));

        if (rd_in == sizeof(size)) // if successfully received size
        {
            // allocate sufficient buffer space
            void *buffer = malloc(size);

            // read serialised measurement object from pipe
            rd_in = read(fd_in, buffer, size);

            // close pipe
            close(fd_in);

            // deserialise from buffer array
            ProtoMeas::Measurements measurements;
            const bool protomeas_ok = measurements.ParseFromArray(buffer, size);

            // free memory
            free(buffer);

            if (protomeas_ok) // if measurement parsed properly from pipe
            {
                // sanitise measurement inputs
                if (CheckValidMeasurement(measurements))
                {
                    std::cout << "Received first measurement." << std::endl;
                    received_first_meas = true;

                    unsigned int num_feature_points = measurements.num_feature_points();

                    // 3D model feature point locations
                    std::vector<cv::Point3d> model_points;
                    // 2D image point locations
                    std::vector<cv::Point2d> image_points;

                    // fill in 2D, 3D points
                    for (unsigned int idx = 0; idx < num_feature_points; idx++)
                    {
                        model_points.push_back(cv::Point3d( measurements.feature_points(idx).x(),
                                                            measurements.feature_points(idx).y(),
                                                            measurements.feature_points(idx).z() ));
                    
                        double az = measurements.bearings(idx).az();
                        double el = measurements.bearings(idx).el();

                        image_points.push_back(cv::Point2d(tan(az)*focal_length, tan(el)*focal_length));
                    }

                    // Output rotation and translation
                    cv::Mat rotation_vector; // Rotation in axis-angle form
                    cv::Mat translation_vector;

                    // Solve for pose
                    cv::solvePnPRansac(model_points, image_points, camera_matrix, dist_coeffs, rotation_vector, translation_vector, true);

                    PoseSolution pose_sol;
                    cv::cv2eigen(translation_vector, pose_sol.pose.pos);
                    cv::Mat R;
                    cv::Rodrigues(rotation_vector, R); // R is 3x3
                    Eigen::Matrix3d mat;
                    cv::cv2eigen(R, mat);
                    Eigen::Quaterniond EigenQuat(mat);
                    pose_sol.pose.quat = EigenQuat;

                    Pose conj_pose = Utilities::ConjugatePose(pose_sol.pose);
                    
                    // initialise filter priors
                    // KF priors
                    VectorXd state0 = VectorXd::Zero(kf.num_states_);
                    state0.head(3) = pose_sol.pose.pos;
                    MatrixXd covar0 = 10.0 * MatrixXd::Identity(kf.num_states_, kf.num_states_);
                    covar0(0, 0) = 1.0;
                    covar0(1, 1) = 1.0;
                    covar0(2, 2) = 3.0;
                    kf.SetInitialStateAndCovar(state0, covar0);
                    kf.R_(0, 0) = 1.0;
                    kf.R_(1, 1) = 1.0;
                    kf.R_(2, 2) = 3.0;

                    std::cout << "KF Model: " << std::endl;
                    kf.PrintModelMatrices();
                    std::cout << std::endl;

                    // MEKF priors
                    Quaterniond init_quat = pose_sol.pose.quat;
                    Vector3d init_omega = 0.01 * Vector3d::Random();
                    Vector3d init_alpha = 0.1 * Vector3d::Random();
                    MatrixXd init_covar = MatrixXd::Identity(mekf.num_states_, mekf.num_states_);
                    init_covar(0, 0) = 0.1;
                    init_covar(1, 1) = 0.1;
                    init_covar(2, 2) = 0.1;
                    mekf.SetInitialStateAndCovar(init_quat, init_omega, init_alpha, init_covar);

                    std::cout << "MEKF Model: " << std::endl;
                    mekf.PrintModelMatrices();
                    std::cout << std::endl;

                    solved_poses.push_back(pose_sol.pose);
                    filtered_poses.push_back(pose_sol.pose);
                    //kf_states.push_back(state0);
                    //kf_covars.push_back(covar0);
                    timestamps.push_back(0.0);
                    init_t = std::chrono::high_resolution_clock::now();
                }
                else // if bad measurement received, wait for new one
                {
                    std::cout << "Received bad first measurement, waiting for valid first measurement." << std::endl;
                }
            }
            else // if message not parsed properly, wait for new one
            {
                std::cout << "Received unparsable first measurement, waiting for valid first measurement." << std::endl;
            }
        }
        else
        { // if no size message found
            //std::cout << "Awaiting first measurement..." << std::endl;
        }
    }

    //------------------------------------------------------------------------/

    // continue once first measurement has been received

    // set-up exit handler
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = exit_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);

    // Output rotation and translation
    cv::Mat rotation_vector; // Rotation in axis-angle form
    cv::Mat translation_vector;

    // open FIFO for measurement read only, without blocking
    fd_in = open(fifo_path_input, O_RDONLY | O_NONBLOCK);

    auto last_t = std::chrono::high_resolution_clock::now();
    auto curr_t = std::chrono::high_resolution_clock::now();

    //-- Main Loop -----------------------------------------------------------/
    while (true)
    {
        // TIMING : run dynamics at specified kf_dt rate
        double curr_delta_t = 0.0;
        while (curr_delta_t < kf_dt)
        {
            curr_t = std::chrono::high_resolution_clock::now();
            curr_delta_t = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(curr_t - last_t).count();
            curr_delta_t *= pow(10.0, -9.0);
        }
        //std::cout << "Predict." << std::endl;
        last_t = curr_t;

        // KF prediction step
        kf.Predict(VectorXd::Zero(kf.num_inputs_));

        // MEKF prediction step (state propagation in terms of quaternions, covariance propagation in terms of gibbs vector)
        mekf.Predict();

        PoseSolution pose_sol;

        //-- Check for new measurement ---------------------------------------/
        // read byte size of measurement object from pipe
        size_t size;
        rd_in = read(fd_in, &size, sizeof(size));

        if (rd_in == sizeof(size)) // if successfully received size
        {
            // allocate sufficient buffer space
            void *buffer = malloc(size);

            // read serialised measurement object from pipe
            rd_in = read(fd_in, buffer, size);

            // deserialise from buffer array
            ProtoMeas::Measurements measurements;
            const bool protomeas_ok = measurements.ParseFromArray(buffer, size);

            // free memory
            free(buffer);

            if (protomeas_ok) // if measurement parsed properly from pipe
            {
                // sanitise measurement inputs
                if (CheckValidMeasurement(measurements))
                {
                    unsigned int num_feature_points = measurements.num_feature_points();

                    // 3D model feature point locations
                    std::vector<cv::Point3d> model_points;
                    // 2D image point locations
                    std::vector<cv::Point2d> image_points;

                    // fill in 2D, 3D points
                    for (unsigned int idx = 0; idx < num_feature_points; idx++)
                    {
                        model_points.push_back(cv::Point3d( measurements.feature_points(idx).x(),
                                                            measurements.feature_points(idx).y(),
                                                            measurements.feature_points(idx).z() ));
                    
                        double az = measurements.bearings(idx).az();
                        double el = measurements.bearings(idx).el();

                        image_points.push_back(cv::Point2d(tan(az)*focal_length, tan(el)*focal_length));
                    }

                    //-- Measurement Update ----------------------------------/

                    // Solve for pose
                    cv::solvePnPRansac(model_points, image_points, camera_matrix, dist_coeffs, rotation_vector, translation_vector, true);

                    cv::cv2eigen(translation_vector, pose_sol.pose.pos);
                    cv::Mat R;
                    cv::Rodrigues(rotation_vector, R); // R is 3x3
                    Eigen::Matrix3d mat;
                    cv::cv2eigen(R, mat);
                    Eigen::Quaterniond EigenQuat(mat);
                    pose_sol.pose.quat = EigenQuat;

                    Pose conj_pose = Utilities::ConjugatePose(pose_sol.pose);
                    
                    // wrap NLS position solution as KF measurement
                    Vector3d pos_meas_wrapper = pose_sol.pose.pos;

                    // KF measurement update step
                    kf.Update(pos_meas_wrapper);

                    // wrap NLS attitude solution as MEKF measurement
                    VectorXd att_meas_wrapper(4);
                    att_meas_wrapper(0) = pose_sol.pose.quat.normalized().w();
                    att_meas_wrapper(1) = pose_sol.pose.quat.normalized().x();
                    att_meas_wrapper(2) = pose_sol.pose.quat.normalized().y();
                    att_meas_wrapper(3) = pose_sol.pose.quat.normalized().z();
                    
                    // MEKF measurement update step
                    mekf.Update(att_meas_wrapper);

                    // MEKF reset step
                    mekf.Reset();
                }
                else // if bad measurement received, skip measurement update
                {
                    std::cout << "Received bad measurement, skipping measurement update." << std::endl;
                }
            }
            else // if message not parsed properly, skip measurement update
            {
                std::cout << "Received unparsable measurement, skipping measurement update." << std::endl;
            }
        }
        else // if no size message found, skip measurement update
        {
        }

        kf.StoreAndClean();
        mekf.StoreAndClean();

        Pose pose_filtered;
        pose_filtered.pos = kf.last_state_estimate.head(3);
        pose_filtered.quat = mekf.quat_est_.normalized();

        curr_t = std::chrono::high_resolution_clock::now();
        curr_elapsed_t = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(curr_t - init_t).count();
        curr_elapsed_t *= pow(10.0, -9.0);

        //--------------------------------------------------------------------/

        //-- Data Storage ----------------------------------------------------/        
        solved_poses.push_back(pose_sol.pose);
        filtered_poses.push_back(pose_filtered);
        //kf_states.push_back(kf.last_state_estimate);
        //kf_covars.push_back(kf.last_covar_estimate);
        timestamps.push_back(curr_elapsed_t);

        // set NLS initial guess for next time-step to latest filtered estimate
        pose0 = filtered_poses.back();

        //-- Write Pose to Pipe ----------------------------------------------/
        if (output_to_pipe)
        {
            ProtoPose::Pose proto_pose;
            ProtoPose::Position *pos = proto_pose.mutable_pos();
            ProtoPose::Attitude *att = proto_pose.mutable_att();
            pos->set_x(pose_filtered.pos(0));
            pos->set_y(pose_filtered.pos(1));
            pos->set_z(pose_filtered.pos(2));
            att->set_qw(pose_filtered.quat.w());
            att->set_qx(pose_filtered.quat.x());
            att->set_qy(pose_filtered.quat.y());
            att->set_qz(pose_filtered.quat.z());
            proto_pose.set_time_stamp(curr_elapsed_t);

            // store byte size of pose object
            size_t size_out = proto_pose.ByteSize();

            // allocate sufficient buffer space
            void *buffer_out = malloc(size_out);

            // serialise to the buffer array
            proto_pose.SerializeToArray(buffer_out, size_out);

            // Open FIFO for write only, without blocking
            fd_out = open(fifo_path_output, O_WRONLY | O_NONBLOCK);

            // write size of pose object to pipe
            write(fd_out, &size_out, sizeof(size_out));

            // write serialised pose object to pipe
            write(fd_out, buffer_out, size_out);

            // close FIFO
            close(fd_out);

            // free memory
            free(buffer_out);
        }

        //-- Handling for Periodic Logging -----------------------------------/
        if (log_to_file && log_periodically && filtered_poses.size() >= vector_reserve_size)
        {
            // write to csv files
            bool append_mode = true;
            Utilities::WritePosesToCSV(solved_poses, prefix + "solved_poses" + postfix, append_mode);
            Utilities::WritePosesToCSV(filtered_poses, prefix + "filtered_poses" + postfix, append_mode);
            //Utilities::WriteKFStatesToCSV(kf_states, prefix + "kf_states" + postfix, append_mode);
            //Utilities::WriteKFCovarsToCSV(kf_covars, prefix + "kf_covars" + postfix, append_mode);
            Utilities::WriteTimestampsToFile(timestamps, prefix + "timestamps" + postfix, append_mode);

            // clear vectors
            solved_poses.clear();
            filtered_poses.clear();
            //kf_states.clear();
            //kf_covars.clear();
            timestamps.clear();
        }

        //-- Handling for Program Exit ---------------------------------------/
        if (finished)
        {
            if (log_to_file)
            {
                bool append_mode;

                if (log_periodically)
                {
                    append_mode = true;
                }
                else
                {
                    append_mode = false;
                }

                // write to csv files
                Utilities::WritePosesToCSV(solved_poses, prefix + "solved_poses" + postfix, append_mode);
                Utilities::WritePosesToCSV(filtered_poses, prefix + "filtered_poses" + postfix, append_mode);
                //Utilities::WriteKFStatesToCSV(kf_states, prefix + "kf_states" + postfix, append_mode);
                //Utilities::WriteKFCovarsToCSV(kf_covars, prefix + "kf_covars" + postfix, append_mode);
                Utilities::WriteTimestampsToFile(timestamps, prefix + "timestamps" + postfix, append_mode);
                printf("Logged data to file.\n");
            }

            // close pipe
            close(fd_in);

            printf("Exiting....\n");
            exit(1);
        }
        //--------------------------------------------------------------------/
    }

    return 0;
}

/**
 * @function CheckValidMeasurement
 * @brief checks whether measurement input sizing is consistent
 * @return true if consistent, false if not
 */
bool CheckValidMeasurement(const ProtoMeas::Measurements& measurements)
{
    unsigned int num_feature_points = measurements.num_feature_points();
    unsigned int verify_num_ft_pt_1 = measurements.feature_points_size();
    unsigned int verify_num_ft_pt_2 = measurements.bearings_size();

    if (num_feature_points != verify_num_ft_pt_1 || num_feature_points != verify_num_ft_pt_2)
    {
        return false;
    }
    else if ((num_feature_points == 0) || (verify_num_ft_pt_1 == 0) || (verify_num_ft_pt_2 == 0))
    {
        return false;
    }
    else
    {
        return true;
    }
}