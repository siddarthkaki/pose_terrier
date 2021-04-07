/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <Eigen/Core>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <math.h>

#include "ceres/ceres.h"
//#include "glog/logging.h"

#include "Utilities.h"
#include "PoseSolver.h"
#include "MEKF2.h"

#include "third_party/json.hpp"
#include "third_party/CppRot/cpprot.h"

using Eigen::AngleAxisd;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXd;
using nlohmann::json;

#define GET_VAR_NAME(Variable) (#Variable)

/**
 * @function main
 * @brief main function
 */
int main(int argc, char **argv)
{

    // timing
    std::chrono::high_resolution_clock::time_point t_program_1 = std::chrono::high_resolution_clock::now();

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

    // specify rigid position vector of camera wrt chaser in chaser frame
    Vector3d rCamVec;
    for (unsigned int idx = 0; idx < 3; idx++)
    {
        rCamVec(idx) = json_params["rCamVec"].at(idx);
    }

    // specify camera focal length
    double focal_length = json_params["focal_length"]; //5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    double bearing_meas_std = double(json_params["bearing_meas_std_deg"]) * Utilities::DEG2RAD;

    const double kf_dt = json_params["kf_dt"];

    // specify rigid position vector of feature points wrt target in target frame
    unsigned int num_features = json_params["rFeaMat"].size();
    MatrixXd rFeaMat(num_features, 3);
    for (unsigned int idx = 0; idx < num_features; idx++)
    {
        for (unsigned int jdx = 0; jdx < 3; jdx++)
        {
            rFeaMat(idx, jdx) = json_params["rFeaMat"][idx]["fea" + std::to_string(idx + 1)][jdx];
        }
    }
    
    /**
    num_features = 11;
    rFeaMat = 2.5 * MatrixXd::Random(num_features, 3);
    **/

    unsigned int num_poses_test = json_params["num_poses_test"];

    //------------------------------------------------------------------------/

    //-- Init KFs ------------------------------------------------------------/

    double mekf_process_noise_std = 0.01;
    double mekf_measurement_noise_std = 0.05;
    double mekf_dt = kf_dt;

    MEKF2::MEKF2 mekf(mekf_dt);
    mekf.Init(mekf_process_noise_std, mekf_measurement_noise_std, mekf_dt);

    //------------------------------------------------------------------------/

    //-- Loop ----------------------------------------------------------------/

    std::vector<Pose> true_poses, solved_poses, solved_poses_conj, filtered_poses;
    std::vector<VectorXd> kf_states;
    std::vector<MatrixXd> kf_covars;
    std::vector<double> solution_times; // [ms]
    std::vector<double> pos_scores;
    std::vector<double> att_scores;

    true_poses.reserve(num_poses_test);
    solved_poses.reserve(num_poses_test);
    solved_poses_conj.reserve(num_poses_test);
    filtered_poses.reserve(num_poses_test);
    kf_states.reserve(num_poses_test);
    kf_covars.reserve(num_poses_test);
    solution_times.reserve(num_poses_test);
    pos_scores.reserve(num_poses_test);
    att_scores.reserve(num_poses_test);

    // initial pose guess
    Pose pose0;

    // true pose
    Pose pose_true;
    //pose_true.quat = Quaterniond::UnitRandom();
    pose_true.quat.w() = 1.0;
    pose_true.quat.vec() = Vector3d::Zero();
    pose_true.pos << 0.5, -0.25, 30.0;

    bool first_run = true;

    for (unsigned int pose_idx = 0; pose_idx < num_poses_test; pose_idx++)
    {
        //-- Simulate Measurements -------------------------------------------/

        // generate true pose values for ith run
        Quaterniond quat_step = AngleAxisd( 0.001, Vector3d::UnitX()) *
                                AngleAxisd(-0.001, Vector3d::UnitY()) *
                                AngleAxisd( 0.001, Vector3d::UnitZ());
        pose_true.quat = pose_true.quat * quat_step;

        pose_true.pos(0) += 0.01;
        pose_true.pos(1) -= 0.01;
        pose_true.pos(2) += 0.05;

        // express feature points in chaser frame at the specified pose
        MatrixXd rMat = Utilities::FeaPointsTargetToChaser(pose_true, rCamVec, rFeaMat);

        // generate simulated measurements
        VectorXd yVec = Utilities::SimulateMeasurements(rMat, focal_length);

        // add Gaussian noise to simulated measurements
        VectorXd yVecNoise = Utilities::AddGaussianNoiseToVector(yVec, bearing_meas_std);

        //--------------------------------------------------------------------/

        //-- Solve for pose --------------------------------------------------/

        // timing
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // if first timestep, set NLS initial guess to default
        if (first_run)
        {
            pose0.pos << 0.0, 0.0, 25.0;
            pose0.quat.w() = 1.0;
            pose0.quat.vec() = Vector3d::Zero();
        }
        else // else, set NLS initial guess to last filtered estimate
        {
            pose0 = filtered_poses.back();
        }

        // solve for pose with ceres (via wrapper)
        PoseSolution pose_sol = PoseSolver::SolvePoseReinit(pose0, yVecNoise, rCamVec, rFeaMat, bearing_meas_std);

        //Pose conj_pose_temp = Utilities::ConjugatePose(pose_sol.pose);
        //Pose conj_pose = PoseSolver::SolvePose(conj_pose_temp, yVecNoise, rCamVec, rFeaMat).pose;

        Pose pose_filtered;

        // if first time-step, then set KF priors
        if (first_run)
        {
            // MEKF priors
            Quaterniond init_quat = pose_sol.pose.quat;
            Vector3d init_omega = 0.01*Vector3d::Ones();
            Vector3d init_alpha = 0.1*Vector3d::Ones();
            MatrixXd init_covar = MatrixXd::Identity(mekf.num_states_, mekf.num_states_);
            VectorXd x0 = VectorXd::Zero(mekf.num_pos_states_);
            x0.head(3) = pose_sol.pose.pos;

            mekf.SetInitialStateAndCovar(init_quat, init_omega, init_alpha, x0, init_covar);

            mekf.PrintModelMatrices();

            pose_filtered.pos = pose_sol.pose.pos;
            pose_filtered.quat = pose_sol.pose.quat;
        }
        else // else, perform KF tracking
        {

            // MEKF prediction step (state propagation in terms of quaternions, covariance propagation in terms of gibbs vector)
            mekf.Predict();

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
            mekf.StoreAndClean();
            
            // storage
            pose_filtered.pos = mekf.pos_est_;//.head(3);
            pose_filtered.quat = mekf.quat_est_.normalized();
        }

        // timing
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

        // time taken to perform pose solution
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
        //--------------------------------------------------------------------/

        //-- Performance Metrics & Storage -----------------------------------/

        // compute position and attitude scores
        double pos_score = Utilities::PositionScore(pose_true.pos, pose_filtered.pos);
        double att_score = Utilities::AttitudeScore(pose_true.quat, pose_filtered.quat);
        //double conj_att_score = Utilities::AttitudeScore(pose_true.quat,     conj_pose.quat);

        // store info from ith run
        true_poses.push_back(pose_true);
        solved_poses.push_back(pose_sol.pose);
        filtered_poses.push_back(pose_filtered);
        solution_times.push_back((double)duration);
        pos_scores.push_back(pos_score);
        att_scores.push_back(att_score); //std::min(att_score,conj_att_score) );

        if (first_run)
        {
            first_run = false;
        }
    }

    // timing
    std::chrono::high_resolution_clock::time_point t_program_2 = std::chrono::high_resolution_clock::now();

    // time taken to perform pose solution
    auto program_duration = std::chrono::duration_cast<std::chrono::milliseconds>(t_program_2 - t_program_1).count();

    //-- Performance Metric Stats & Output -----------------------------------/

    // write to csv file
    Utilities::WritePosesToCSV(true_poses, Utilities::WrapVarToPath(std::string(GET_VAR_NAME(true_poses))), false);
    Utilities::WritePosesToCSV(solved_poses, Utilities::WrapVarToPath(std::string(GET_VAR_NAME(solved_poses))), false);
    Utilities::WritePosesToCSV(filtered_poses, Utilities::WrapVarToPath(std::string(GET_VAR_NAME(filtered_poses))), false);
    
    double pos_score_mean = Utilities::StdVectorMean(pos_scores);
    double att_score_mean = Utilities::StdVectorMean(att_scores);

    double pos_score_std = sqrt(Utilities::StdVectorVar(pos_scores));
    double att_score_std = sqrt(Utilities::StdVectorVar(att_scores));

    double solution_times_mean = Utilities::StdVectorMean(solution_times);

    std::cout << "num_runs :\t" << num_poses_test << std::endl
              << std::endl;

    std::cout << "pos_score_mean :\t" << pos_score_mean << " [m]" << std::endl;
    std::cout << "pos_score_std  :\t" << pos_score_std << " [m]" << std::endl
              << std::endl;

    std::cout << "att_score_mean :\t" << att_score_mean * Utilities::RAD2DEG << " [deg]" << std::endl;
    std::cout << "att_score_std  :\t" << att_score_std * Utilities::RAD2DEG << " [deg]" << std::endl
              << std::endl;

    std::cout << "mean_time  :\t" << solution_times_mean << " [ms]" << std::endl;
    std::cout << "total_time :\t" << solution_times_mean * num_poses_test << " [ms]" << std::endl;

    std::cout << "total_program_time :\t" << (double)program_duration << " [ms]" << std::endl;
    
    //------------------------------------------------------------------------/

    return 0;
}