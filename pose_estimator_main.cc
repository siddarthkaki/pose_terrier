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
//#include "glog/logging.h"

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

/**
 * @function main
 * @brief main function
 */
int main(int argc, char **argv)
{
    //google::InitGoogleLogging(argv[0]);

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

    std::string pipe_path_input  = json_params["pipe_path_input"];
    std::string pipe_path_output = json_params["pipe_path_output"];

    // specify rigid position vector of camera wrt chaser in chaser frame
    Vector3d rCamVec;
    for (unsigned int idx = 0; idx < 2; idx++)
    {
        rCamVec(idx) = json_params["rCamVec"].at(idx);
    }

    // specify camera focal length
    double focal_length = json_params["focal_length"]; //5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    double meas_std = double(json_params["meas_std_deg"]) * Utilities::DEG2RAD;

    // specify expected number of time-steps for memory pre-allocation
    unsigned int num_poses_test = json_params["num_poses_test"];

    double kf_dt = json_params["kf_dt"];

    //------------------------------------------------------------------------/

    //-- Init sequence -------------------------------------------------------/

    // declare vectors for storage
    std::vector<Pose> solved_poses, solved_poses_conj, filtered_poses;
    std::vector<double> solution_times, timestamps; // [ms]

    // pre-allocate memory
    solved_poses.reserve(num_poses_test);
    solved_poses_conj.reserve(num_poses_test);
    filtered_poses.reserve(num_poses_test);
    solution_times.reserve(num_poses_test);
    timestamps.reserve(num_poses_test);

    // Kalman Filter object
    KF::KalmanFilter kf;

    // TODO: read in initial guess from json or other program
    // initial pose guess
    Pose pose0;
    pose0.pos << 0.0, 0.0, 25.0;
    pose0.quat.w() = 1.0;
    pose0.quat.vec() = Vector3d::Zero();

    // initialise pipe
    int fd_in, rd_in = 0;
    const char *fifo_path_input = pipe_path_input.c_str();
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

            // deserialise from buffer array
            Proto::Measurements measurements;
            measurements.ParseFromArray(buffer, size);

            // close pipe
            close(fd_in);

            std::cout << "Received first measurement." << std::endl;
            received_first_meas = true;

            unsigned int num_feature_points = measurements.num_feature_points();

            // Construct Eigen::MatrixXd out of feature point locations
            MatrixXd rFeaMat(num_feature_points, 3);
            for (unsigned int idx = 0; idx < num_feature_points; idx++)
            {
                rFeaMat(idx,0) = measurements.feature_points(idx).x();
                rFeaMat(idx,1) = measurements.feature_points(idx).y();
                rFeaMat(idx,2) = measurements.feature_points(idx).z();
            }
            
            // Construct Eigen::VectorXd out of measurements
            VectorXd yVec(2*num_feature_points);
            for (unsigned int idx = 0; idx < num_feature_points; idx++)
            {
                yVec(2*idx+0) = measurements.bearings(idx).az();
                yVec(2*idx+1) = measurements.bearings(idx).el();
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

            filtered_poses.push_back(pose_sol.pose);
            timestamps.push_back(0.0);
        }
        else
        { // if no size message found
            std::cout << "Awaiting first measurement..." << std::endl;
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


    // Open FIFO for read only, without blocking
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
            curr_delta_t = (double)std::chrono::duration_cast<std::chrono::seconds>(curr_t - last_t).count();
        }
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

            // read serialised pose object from pipe
            rd_in= read(fd_in, buffer, size);

            // deserialise from buffer array
            Proto::Measurements measurements;
            measurements.ParseFromArray(buffer, size);

            unsigned int num_feature_points = measurements.num_feature_points();

            // Construct Eigen::MatrixXd out of feature point locations
            MatrixXd rFeaMat(num_feature_points, 3);
            for (unsigned int idx = 0; idx < num_feature_points; idx++)
            {
                rFeaMat(idx,0) = measurements.feature_points(idx).x();
                rFeaMat(idx,1) = measurements.feature_points(idx).y();
                rFeaMat(idx,2) = measurements.feature_points(idx).z();
            }
            
            // Construct Eigen::VectorXd out of measurements
            VectorXd yVec(2*num_feature_points);
            for (unsigned int idx = 0; idx < num_feature_points; idx++)
            {
                yVec(2*idx+0) = measurements.bearings(idx).az();
                yVec(2*idx+1) = measurements.bearings(idx).el();
            }

            //-- Measurement update ------------------------------------------/

            // set NLS initial guess to last filtered estimate
            pose0 = filtered_poses.back();

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
        else
        { } // if no size message found, skip measurement update

        kf.StoreAndClean();

        Pose pose_filtered;
        VectorXd pose_filt_wrapper = kf.states.back();
        pose_filtered.pos = pose_filt_wrapper.head(3);
        pose_filtered.quat = AngleAxisd(pose_filt_wrapper( 9), Vector3d::UnitX()) *
                             AngleAxisd(pose_filt_wrapper(10), Vector3d::UnitY()) *
                             AngleAxisd(pose_filt_wrapper(11), Vector3d::UnitZ());

        //--------------------------------------------------------------------/

        //-- Performance Metrics & Storage -----------------------------------/

        /*
        // store info from ith run
        true_poses.push_back(pose_true);
        solved_poses.push_back(pose_sol.pose);
        solved_poses_conj.push_back(conj_pose);
        filtered_poses.push_back(pose_filtered);
        solution_times.push_back((double)duration);
        pos_scores.push_back(pos_score);
        att_scores.push_back(att_score); //std::min(att_score,conj_att_score) );
        */
        solved_poses.push_back(pose_sol.pose);
        filtered_poses.push_back(pose_filtered);

        //-- Handling for program exit ----------------------------------------/
        if (finished)
        {
            // write to csv file
            Utilities::WritePosesToCSV(solved_poses, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(solved_poses))));
            Utilities::WritePosesToCSV(filtered_poses, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(filtered_poses))));
            Utilities::WriteKFStatesToCSV(kf.states, Utilities::WrapVarToPath(std::string("kf_states")));
            Utilities::WriteKFCovarsToCSV(kf.covars, Utilities::WrapVarToPath(std::string("kf_covars")));
            printf("Logged data to file.\nExiting....\n");
            exit(1);
        }
        //--------------------------------------------------------------------/
    }

    return 0;
}