/* Copyright (c) 2022 Siddarth Kaki
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
#include <signal.h>

#include "ceres/ceres.h"
#include "glog/logging.h"

#include "Utilities.h"
#include "PoseSolver.h"
#include "MEKF2.h"
#include "QuateRA.h"
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

    struct sigaction sa;
    sa.sa_handler = SIG_IGN;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    if (sigaction(SIGPIPE, &sa, 0) == -1)
    {
        std::cout << "SIGPIPE ERROR" << std::endl;
    }
    
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

    // specify initial guess of relative position vector of target wrt chaser in chaser frame
    Vector3d rPos0;
    for (unsigned int idx = 0; idx < 3; idx++)
    {
        rPos0(idx) = json_params["rPos0"].at(idx);
    }

    // specify rigid position vector of camera wrt chaser in chaser frame
    Vector3d rCamVec;
    for (unsigned int idx = 0; idx < 3; idx++)
    {
        rCamVec(idx) = json_params["rCamVec"].at(idx);
    }

    // specify camera focal length
    //double focal_length = json_params["focal_length"]; //5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    double bearing_meas_std = double(json_params["bearing_meas_std_deg"]) * Utilities::DEG2RAD;

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
    
    // double kf_process_noise_std = 0.01;
    // double kf_measurement_noise_std = 0.05;

    double mekf_process_noise_std = double(json_params["process_noise_std"]);//0.01;
    double mekf_measurement_noise_std = double(json_params["measurement_noise_std"]);//0.05;
    double tau = double(json_params["tau"]);
    double max_flip_thresh_deg = double(json_params["max_flip_thresh_deg"]);
    double qpsd = double(json_params["qpsd"]);

    MEKF2::MEKF2 mekf(mekf_dt);
    mekf.Init(
        mekf_process_noise_std, 
        mekf_measurement_noise_std, 
        mekf_dt, 
        tau,
        qpsd,  
        max_flip_thresh_deg
    );

    //-- Init QuateRA --------------------------------------------------------/

    bool adapt_window_size = (json_params["adapt_window_size"]);
    unsigned int L = int(json_params["window_size"]);
    unsigned int Lmin = int(json_params["min_window_size"]);
    unsigned int Lmax = int(json_params["max_window_size"]);
    double angle_noise_std_deg =  double(json_params["quatera_angle_noise_std_deg"]);
    double threshold_n = double(json_params["eps_threshold_n"]);

    QuateRA::QuateRA quatera(mekf_dt);
    quatera.Init(
        adapt_window_size,
        L,
        Lmin,
        Lmax,
        angle_noise_std_deg*Utilities::DEG2RAD,
        threshold_n
    );

    //-- Init sequence -------------------------------------------------------/

    // log path name prefixing and postfixing
    std::string init_time_str = std::to_string(std::time(nullptr));
    std::string prefix = "../data/" + init_time_str + "_";
    //std::string prefix = "../data/";
    std::string postfix = ".csv";

    // declare vectors for storage
    std::vector<Pose> solved_poses, filtered_poses;
    std::vector<VectorXd> solved_omegas, filtered_omegas, filtered_alphas, filtered_pos_states;
    std::vector<VectorXd> filtered_covar_diag;
    std::vector<double> timestamps; // [s]

    // pre-allocate memory
    solved_poses.reserve(vector_reserve_size);
    filtered_poses.reserve(vector_reserve_size);
    solved_omegas.reserve(vector_reserve_size);
    filtered_omegas.reserve(vector_reserve_size);
    filtered_alphas.reserve(vector_reserve_size);
    filtered_pos_states.reserve(vector_reserve_size);
    filtered_covar_diag.reserve(vector_reserve_size);
    timestamps.reserve(vector_reserve_size);

    // clock object
    auto init_t = std::chrono::high_resolution_clock::now();
    double curr_elapsed_t = 0.0;

    // initial pose guess
    Pose pose0;
    pose0.pos = rPos0;
    pose0.quat.w() = 1.0;
    pose0.quat.vec() = Vector3d::Zero();

    std::cout << "Initial pos guess: " << pose0.pos.transpose() << std::endl;
    std::cout << "Initial att guess: " << pose0.quat.w() << " " << pose0.quat.vec().transpose() << std::endl << std::endl;

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
    
    // open FIFO for read only, with blocking to receive first measurement
    fd_in = open(fifo_path_input, O_RDONLY);

    // loop for waiting for first measurement to initialise filters
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
                    bool valid_pose = true;
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
                    PoseSolution pose_sol = PoseSolver::SolvePoseReinit(pose0, yVec, rCamVec, rFeaMat, bearing_meas_std);

                    // check if pose solution is valid
                    if (abs(pose_sol.pose.quat.norm() - 1.0) < 1e-4)
                    {
                        // Pose conj_pose_temp = Utilities::ConjugatePose(pose_sol.pose);
                        // Pose conj_pose = PoseSolver::SolvePose(conj_pose_temp, yVec, rCamVec, rFeaMat, bearing_meas_std).pose;

                        // MEKF priors
                        Quaterniond init_quat = pose_sol.pose.quat;
                        Vector3d init_omega = 0.005 * Vector3d::Random();
                        Vector3d init_alpha = 0.01 * Vector3d::Random();
                        MatrixXd init_covar = MatrixXd::Identity(mekf.num_states_, mekf.num_states_);
                        VectorXd x0 = VectorXd::Zero(mekf.num_pos_states_);
                        x0.head(3) = pose_sol.pose.pos;
                        x0.segment(3,3) = 0.005 * Vector3d::Random();
                        x0.tail(3) = 0.001 * Vector3d::Random();
                        mekf.SetInitialStateAndCovar(init_quat, init_omega, init_alpha, x0, init_covar);

                        std::cout << "init_omega: " << init_omega << std::endl << std::endl;
                        std::cout << "Q_att: " << std::endl << mekf.Q_.topLeftCorner(9, 9) << std::endl << std::endl;

                        // QuateRA init measurement
                        quatera.InitMeasurement(Utilities::QuatToVec4(init_quat), pose_sol.cov_pose.topLeftCorner(3,3));

                        solved_poses.push_back(pose_sol.pose);
                        filtered_poses.push_back(pose_sol.pose);
                        solved_omegas.push_back(Vector3d::Zero());
                        filtered_omegas.push_back(init_omega);
                        filtered_alphas.push_back(init_alpha);
                        filtered_pos_states.push_back(mekf.state_est_.tail(9));
                        filtered_covar_diag.push_back(pose_sol.cov_pose.diagonal());
                        timestamps.push_back(0.0);
                        init_t = std::chrono::high_resolution_clock::now();
                    }
                    else // reject the measurement
                    {
                        std::cout << std::fixed << std::setprecision(9);
                        std::cout << "Invalid pose solution; skipping measurement. Quat norm: " << pose_sol.pose.quat.norm() << std::endl;
                        received_first_meas = false;
                        valid_pose = false;
                    }
                    
                    // write pose to pipe 
                    if (output_to_pipe)
                    {
                        ProtoPose::Pose proto_pose;
                        ProtoPose::Position *pos = proto_pose.mutable_pos();
                        ProtoPose::Attitude *att = proto_pose.mutable_att();
                        pos->set_x(pose_sol.pose.pos(0));
                        pos->set_y(pose_sol.pose.pos(1));
                        pos->set_z(pose_sol.pose.pos(2));
                        att->set_qw(pose_sol.pose.quat.w());
                        att->set_qx(pose_sol.pose.quat.x());
                        att->set_qy(pose_sol.pose.quat.y());
                        att->set_qz(pose_sol.pose.quat.z());
                        proto_pose.set_time_stamp(0.0);
                        proto_pose.set_valid_pose(valid_pose);

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
     // close pipe
    close(fd_in);
    // set-up exit handler
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = exit_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);

    // open FIFO for measurement read only, without blocking
    fd_in = open(fifo_path_input, O_RDONLY | O_NONBLOCK);

    // set NLS initial guess for next time-step to initial filtered estimate
    pose0 = filtered_poses.back();

    std::cout << "Running estimator..." << std::endl;

    auto last_t = std::chrono::high_resolution_clock::now();
    auto curr_t = std::chrono::high_resolution_clock::now();

    unsigned int meas_counter = 1;

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

        last_t = curr_t;

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
                    pose_sol = PoseSolver::SolvePoseReinit(pose0, yVec, rCamVec, rFeaMat, bearing_meas_std);

                    // check if pose solution is valid
                    if (abs(pose_sol.pose.quat.norm() - 1.0) < 1e-4)
                    {

                        meas_counter++;

                        if (meas_counter % 2 != 0)
                        {
                            // wrap NLS pose solution as MEKF measurement
                            VectorXd meas_wrapper(7);
                            meas_wrapper(0) = pose_sol.pose.quat.normalized().w();
                            meas_wrapper(1) = pose_sol.pose.quat.normalized().x();
                            meas_wrapper(2) = pose_sol.pose.quat.normalized().y();
                            meas_wrapper(3) = pose_sol.pose.quat.normalized().z();
                            meas_wrapper.tail(3) = pose_sol.pose.pos;
                            
                            // MEKF measurement update step
                            mekf.R_ = pose_sol.cov_pose;
                            mekf.Update(meas_wrapper);

                            // MEKF reset step
                            mekf.Reset();
                        }
                        else
                        {
                            // QuateRA measurement
                            curr_t = std::chrono::high_resolution_clock::now();
                            curr_elapsed_t = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(curr_t - init_t).count();
                            curr_elapsed_t *= pow(10.0, -9.0);
                            
                            quatera.angle_noise_std_ = sqrt(mekf.covar_est_.diagonal().head(3).mean());

                            quatera.Update(Utilities::QuatToVec4(pose_sol.pose.quat), pose_sol.cov_pose.topLeftCorner(3,3), curr_elapsed_t);

                            mekf.AngVelUpdate(quatera.ang_vel_est_, quatera.covar_est_.bottomRightCorner(3,3));
                        }

                        //quatera.eps_mean_ = eps_mean * pow(euler_noise_std*Utilities::RAD2DEG/2.0, 2);
                        //quatera.eps_std_ = eps_std * pow(euler_noise_std*Utilities::RAD2DEG/2.0, 2);

                        //std::cout << quatera.eps_mean_ << " eps mean" << std::endl;
                        //std::cout << quatera.eps_std_ << " eps std" << std::endl << std::endl;

                        /*
                        std::cout << Utilities::RAD2DEG*mekf.omega_est_.transpose() << std::endl;
                        std::cout << Utilities::RAD2DEG*quatera.ang_vel_est_.transpose() << std::endl << std::endl;

                        std::cout << Utilities::RAD2DEG*mekf.omega_est_.norm() << std::endl;
                        std::cout << Utilities::RAD2DEG*quatera.ang_vel_est_.norm() << std::endl << std::endl;
                        */
                    
                        std::cout << curr_elapsed_t << " sec" << std::endl << std::endl;
                    }
                    else // reject the measurement
                    {
                        std::cout << std::fixed << std::setprecision(9);
                        std::cout << "Invalid pose solution; skipping measurement. Quat norm: " << pose_sol.pose.quat.norm() << std::endl;
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

        mekf.StoreAndClean();

        Pose pose_filtered;
        pose_filtered.pos = mekf.pos_est_;
        pose_filtered.quat = mekf.quat_est_.normalized();

        Vector3d alpha_filtered = mekf.state_est_.segment(6, 3);

        // VectorXd covar_filtered_diag = mekf.covar_est_.diagonal();
        // Vector6d pose_covar_filtered_diag;
        // pose_covar_filtered_diag << pose_covar_filtered_diag.segment(0,3), pose_covar_filtered_diag.segment(9,3);

        curr_t = std::chrono::high_resolution_clock::now();
        curr_elapsed_t = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(curr_t - init_t).count();
        curr_elapsed_t *= pow(10.0, -9.0);

        //--------------------------------------------------------------------/

        //-- Data Storage ----------------------------------------------------/        
        solved_poses.push_back(pose_sol.pose);
        filtered_poses.push_back(pose_filtered);
        solved_omegas.push_back(quatera.ang_vel_est_);
        filtered_omegas.push_back(mekf.omega_est_);
        filtered_alphas.push_back(alpha_filtered);
        filtered_pos_states.push_back(mekf.state_est_.tail(9));
        filtered_covar_diag.push_back(mekf.covar_est_.diagonal());
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
            proto_pose.set_valid_pose(true);
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
            // write to csv files
            Utilities::WritePosesToCSV(solved_poses, prefix + "solved_poses" + postfix, append_mode);
            Utilities::WritePosesToCSV(filtered_poses, prefix + "filtered_poses" + postfix, append_mode);
            Utilities::WriteKFStatesToCSV(solved_omegas, prefix + "solved_omegas" + postfix, append_mode);
            Utilities::WriteKFStatesToCSV(filtered_omegas, prefix + "filtered_omegas" + postfix, append_mode);
            Utilities::WriteKFStatesToCSV(filtered_alphas, prefix + "filtered_alphas" + postfix, append_mode);
            Utilities::WriteKFStatesToCSV(filtered_pos_states, prefix + "filtered_pos_states" + postfix, append_mode);
            Utilities::WriteKFStatesToCSV(filtered_covar_diag, prefix + "filtered_covar_diag" + postfix, append_mode);
            Utilities::WriteTimestampsToFile(timestamps, prefix + "timestamps" + postfix, append_mode);
            Utilities::WriteTimestampsToFile(quatera.L_history, prefix + "window_size" + postfix, append_mode);
            Utilities::WriteTimestampsToFile(quatera.eps_history, prefix + "eps_history" + postfix, append_mode);
                
            // clear vectors
            solved_poses.clear();
            filtered_poses.clear();
            solved_omegas.clear();
            filtered_omegas.clear();
            filtered_alphas.clear();
            filtered_pos_states.clear();
            filtered_covar_diag.clear();
            timestamps.clear();
            quatera.L_history.clear();
            quatera.eps_history.clear();
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
                Utilities::WriteKFStatesToCSV(solved_omegas, prefix + "solved_omegas" + postfix, append_mode);
                Utilities::WriteKFStatesToCSV(filtered_omegas, prefix + "filtered_omegas" + postfix, append_mode);
                Utilities::WriteKFStatesToCSV(filtered_alphas, prefix + "filtered_alphas" + postfix, append_mode);
                Utilities::WriteKFStatesToCSV(filtered_pos_states, prefix + "filtered_pos_states" + postfix, append_mode);
                Utilities::WriteKFStatesToCSV(filtered_covar_diag, prefix + "filtered_covar_diag" + postfix, append_mode);
                Utilities::WriteTimestampsToFile(timestamps, prefix + "timestamps" + postfix, append_mode);
                Utilities::WriteTimestampsToFile(quatera.L_history, prefix + "window_size" + postfix, append_mode);
                Utilities::WriteTimestampsToFile(quatera.eps_history, prefix + "eps_history" + postfix, append_mode);
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