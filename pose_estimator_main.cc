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

#include "ceres/ceres.h"
#include "glog/logging.h"

#include "cost_functor.h"
#include "Utilities.h"
#include "PoseSolver.h"
#include "KalmanFilter.h"
#include "pose.pb.h"
#include "measurement.pb.h"

#include "third_party/json.hpp"

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
    for (unsigned int idx = 0; idx < 2; idx++)
    {
        rCamVec(idx) = json_params["rCamVec"].at(idx);
    }

    // specify camera focal length
    //double focal_length = json_params["focal_length"]; //5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    //double meas_std = double(json_params["meas_std_deg"]) * Utilities::DEG2RAD;

    // specify expected number of time-steps for memory pre-allocation
    const unsigned int num_poses_test = json_params["num_poses_test"];

    const bool output_to_pipe = json_params["output_to_pipe"];
    const bool log_to_file = json_params["log_to_file"];
    const bool log_periodically = json_params["log_periodically"];
    const unsigned int vector_reserve_size = json_params["vector_reserve_size"];

    const double kf_dt = json_params["kf_dt"];

    //------------------------------------------------------------------------/

    //-- Init sequence -------------------------------------------------------/

    // declare vectors for storage
    std::vector<Pose> solved_poses, filtered_poses;
    std::vector<VectorXd> kf_states;
    std::vector<MatrixXd> kf_covars;
    std::vector<double> solution_times, timestamps; // [ms]

    // pre-allocate memory
    solved_poses.reserve(vector_reserve_size);
    filtered_poses.reserve(vector_reserve_size);
    kf_states.reserve(vector_reserve_size);
    kf_covars.reserve(vector_reserve_size);
    solution_times.reserve(vector_reserve_size);
    timestamps.reserve(vector_reserve_size);

    // Kalman Filter object
    KF::KalmanFilter kf;

    // TODO: read in initial guess from json or other program
    // initial pose guess
    Pose pose0;
    pose0.pos << 0.0, 0.0, 25.0;
    pose0.quat.w() = 1.0;
    pose0.quat.vec() = Vector3d::Zero();

    std::cout << "Waiting for first measurement." << std::endl;

    // set pipe names
    int fd_in, fd_out, rd_in = 0;
    const char *fifo_path_input = pipe_path_input.c_str();
    const char *fifo_path_output = pipe_path_output.c_str();
    
    // Open FIFO for read only, with blocking to receive first measurement
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

                    // Construct Eigen::MatrixXd out of feature point locations
                    MatrixXd rFeaMat(num_feature_points, 3);
                    for (unsigned int idx = 0; idx < num_feature_points; idx++)
                    {
                        rFeaMat(idx, 0) = measurements.feature_points(idx).x();
                        rFeaMat(idx, 1) = measurements.feature_points(idx).y();
                        rFeaMat(idx, 2) = measurements.feature_points(idx).z();
                    }

                    // Construct Eigen::VectorXd out of measurements
                    VectorXd yVec(2 * num_feature_points);
                    for (unsigned int idx = 0; idx < num_feature_points; idx++)
                    {
                        yVec(2 * idx + 0) = measurements.bearings(idx).az();
                        yVec(2 * idx + 1) = measurements.bearings(idx).el();
                    }

                    // solve for pose with ceres (via wrapper)
                    PoseSolution pose_sol = PoseSolver::SolvePoseReinit(pose0, yVec, rCamVec, rFeaMat);

                    Pose conj_pose_temp = Utilities::ConjugatePose(pose_sol.pose);
                    Pose conj_pose = PoseSolver::SolvePose(conj_pose_temp, yVec, rCamVec, rFeaMat).pose;

                    // initialise KF

                    double kf_process_noise_std = 0.01;
                    double kf_measurement_noise_std = 0.05;

                    kf.InitLinearPoseTracking(kf_process_noise_std, kf_measurement_noise_std, kf_dt);
                    VectorXd state0 = VectorXd::Zero(kf.num_states_);
                    state0.head(3) = pose_sol.pose.pos;
                    state0.segment(3, 3) = pose_sol.pose.quat.toRotationMatrix().eulerAngles(0, 1, 2);
                    MatrixXd covar0 = 10.0 * MatrixXd::Identity(kf.num_states_, kf.num_states_);
                    covar0( 0,  0) = 1.0;
                    covar0( 1,  1) = 1.0;
                    covar0( 2,  2) = 3.0;
                    covar0( 9,  9) = 10.0 * Utilities::DEG2RAD;
                    covar0(10, 10) = 10.0 * Utilities::DEG2RAD;
                    covar0(11, 11) = 10.0 * Utilities::DEG2RAD;
                    kf.SetInitialStateAndCovar(state0, covar0);

                    kf.R_(0, 0) = 1.0;
                    kf.R_(1, 1) = 1.0;
                    kf.R_(2, 2) = 3.0;
                    kf.R_(3, 3) = 10.0 * Utilities::DEG2RAD;
                    kf.R_(4, 4) = 10.0 * Utilities::DEG2RAD;
                    kf.R_(5, 5) = 10.0 * Utilities::DEG2RAD;

                    kf.PrintModelMatrices();

                    solved_poses.push_back(pose_sol.pose);
                    filtered_poses.push_back(pose_sol.pose);
                    kf_states.push_back(state0);
                    kf_covars.push_back(covar0);
                    timestamps.push_back(0.0);
                }
                else // if bad measurement received, wait for new one
                {
                    std::cout << "Received bad first measurement, waiting for valid first measurement." << std::endl;
                }
            }
            else // if message not parsed properly, wait for new one
            {
                std::cout << "Received bad first measurement, waiting for valid first measurement." << std::endl;
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

                    // Construct Eigen::MatrixXd out of feature point locations
                    MatrixXd rFeaMat(num_feature_points, 3);
                    for (unsigned int idx = 0; idx < num_feature_points; idx++)
                    {
                        rFeaMat(idx, 0) = measurements.feature_points(idx).x();
                        rFeaMat(idx, 1) = measurements.feature_points(idx).y();
                        rFeaMat(idx, 2) = measurements.feature_points(idx).z();
                    }

                    // Construct Eigen::VectorXd out of measurements
                    VectorXd yVec(2 * num_feature_points);
                    for (unsigned int idx = 0; idx < num_feature_points; idx++)
                    {
                        yVec(2 * idx + 0) = measurements.bearings(idx).az();
                        yVec(2 * idx + 1) = measurements.bearings(idx).el();
                    }

                    //-- Measurement update ----------------------------------/

                    // Note: NLS initial guess for this time-step (pose0) is set to
                    //       the last filtered estimate by this point

                    // solve for pose with ceres (via wrapper)
                    pose_sol = PoseSolver::SolvePoseReinit(pose0, yVec, rCamVec, rFeaMat);

                    Pose conj_pose_temp = Utilities::ConjugatePose(pose_sol.pose);
                    Pose conj_pose = PoseSolver::SolvePose(conj_pose_temp, yVec, rCamVec, rFeaMat).pose;

                    // wrap NLS pose solution as KF measurement
                    VectorXd pose_meas_wrapper(6);
                    pose_meas_wrapper.head(3) = pose_sol.pose.pos;
                    pose_meas_wrapper.tail(3) = pose_sol.pose.quat.toRotationMatrix().eulerAngles(0, 1, 2);

                    // wrap NLS conjugate pose solution as KF measurement
                    VectorXd conj_pose_meas_wrapper(6);
                    conj_pose_meas_wrapper.head(3) = conj_pose.pos;
                    conj_pose_meas_wrapper.tail(3) = conj_pose.quat.toRotationMatrix().eulerAngles(0, 1, 2);

                    // choose as measurement whichever pose produces the smallest measurement residual norm
                    double pose_meas_norm = (pose_meas_wrapper - kf.H_ * kf.statek1k_).norm();
                    double conj_pose_meas_norm = (conj_pose_meas_wrapper - kf.H_ * kf.statek1k_).norm();
                    if (pose_meas_norm < conj_pose_meas_norm)
                    {
                        kf.Update(pose_meas_wrapper);
                    }
                    else
                    {
                        kf.Update(conj_pose_meas_wrapper);
                    }
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

        Pose pose_filtered;
        VectorXd pose_filt_wrapper = kf.last_state_estimate;
        pose_filtered.pos = pose_filt_wrapper.head(3);
        pose_filtered.quat = (AngleAxisd(pose_filt_wrapper( 9), Vector3d::UnitX()) *
                              AngleAxisd(pose_filt_wrapper(10), Vector3d::UnitY()) *
                              AngleAxisd(pose_filt_wrapper(11), Vector3d::UnitZ()));

        //--------------------------------------------------------------------/

        //-- Data Storage ----------------------------------------------------/
        solved_poses.push_back(pose_sol.pose);
        filtered_poses.push_back(pose_filtered);
        kf_states.push_back(kf.last_state_estimate);
        kf_covars.push_back(kf.last_covar_estimate);

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
            proto_pose.set_time_stamp(0.0); // TODO: TIME-STAMP

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
            // TODO TIMESTAMP FILENAME
            Utilities::WritePosesToCSV(solved_poses, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(solved_poses))), append_mode);
            Utilities::WritePosesToCSV(filtered_poses, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(filtered_poses))), append_mode);
            Utilities::WriteKFStatesToCSV(kf_states, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(kf_states))), append_mode);
            Utilities::WriteKFCovarsToCSV(kf_covars, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(kf_covars))), append_mode);

            // clear vectors
            solved_poses.clear();
            filtered_poses.clear();
            kf_states.clear();
            kf_covars.clear();
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
                // TODO TIMESTAMP FILENAME
                Utilities::WritePosesToCSV(solved_poses, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(solved_poses))), append_mode);
                Utilities::WritePosesToCSV(filtered_poses, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(filtered_poses))), append_mode);
                Utilities::WriteKFStatesToCSV(kf_states, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(kf_states))), append_mode);
                Utilities::WriteKFCovarsToCSV(kf_covars, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(kf_covars))), append_mode);
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