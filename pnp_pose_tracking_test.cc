/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <Eigen/Core>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <math.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>

#include "Utilities.h"
#include "KalmanFilter.h"

#include "third_party/json.hpp"

using Eigen::AngleAxisd;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXd;
using nlohmann::json;

#define GET_VARIABLE_NAME(Variable) (#Variable)

/**
 * @function main
 * @brief main function
 */
int main(int argc, char **argv)
{
    std::srand((unsigned int) time(NULL));

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
    for (unsigned int idx = 0; idx < 2; idx++)
    {
        rCamVec(idx) = json_params["rCamVec"].at(idx);
    }

    // specify camera focal length
    double focal_length = json_params["focal_length"]; //5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    double meas_std = double(json_params["meas_std_deg"]) * Utilities::DEG2RAD;

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
    /*
    // TEMPORARY
    num_features = 11;
    rFeaMat = 2.5 * MatrixXd::Random(num_features, 3);
    std::cout << rFeaMat << std::endl;
    */

    unsigned int num_poses_test = json_params["num_poses_test"];

    //------------------------------------------------------------------------/

    // 3D model points
    std::vector<cv::Point3d> model_points;
    for (unsigned int idx = 0; idx < num_features; idx++)
    {
        model_points.push_back(cv::Point3d(rFeaMat(idx, 0), rFeaMat(idx, 1), rFeaMat(idx, 2)));
    }

    double f = focal_length;
        MatrixXd PMat(3,3);
        PMat << f, 0, 0,
                0, f, 0,
                0, 0, 1;

    cv::Mat camera_matrix;
    cv::eigen2cv(PMat, camera_matrix);

    cv::Mat dist_coeffs = cv::Mat::zeros(4, 1, cv::DataType<double>::type); // Assuming no lens distortion


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

    // Kalman Filter object
    KF::KalmanFilter kf;

    // initial pose guess
    Pose pose0;

    // true pose
    Pose pose_true;
    pose_true.pos << 0.5, -0.25, 30.0;
    //pose_true.quat = Quaterniond::UnitRandom();
    pose_true.quat.w() = 1.0;
    pose_true.quat.vec() = Vector3d::Zero();

    // Output rotation and translation
    cv::Mat rotation_vector; // Rotation in axis-angle form
    cv::Mat translation_vector;

    bool first_run = true;

    for (unsigned int pose_idx = 0; pose_idx < num_poses_test; pose_idx++)
    {
        //-- Simulate Measurements -------------------------------------------/

        // generate true pose values for ith run
        pose_true.pos(0) += 0.001;
        pose_true.pos(1) -= 0.001;
        pose_true.pos(2) += 0.01;
        Quaterniond quat_step = AngleAxisd(0.001, Vector3d::UnitX()) *
                                AngleAxisd(-0.001, Vector3d::UnitY()) *
                                AngleAxisd(0.001, Vector3d::UnitZ());
        pose_true.quat = pose_true.quat * quat_step;

        // express feature points in chaser frame at the specified pose
        MatrixXd rMat = Utilities::FeaPointsTargetToChaser(pose_true, rCamVec, rFeaMat);

        // 2D image points
        std::vector<cv::Point2d> image_points;

        // project feature points to image plane
        for (unsigned int idx = 0; idx < num_features; idx++)
        {
            Vector3d rVeci = rMat.row(idx);
            Vector2d imgPti = Utilities::CameraProjection(rVeci, focal_length);

            image_points.push_back(cv::Point2d(imgPti(0), imgPti(1)));
        }
        // TODO: ADD NOISE!

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

        Pose pose_filtered;

        // if first time-step, then initialise KF model, and set KF prior to first NLS solution
        if (first_run)
        {
            double kf_process_noise_std = 0.01;
            double kf_measurement_noise_std = 0.05;
            double kf_dt = 0.5;

            kf.InitLinearPoseTracking(kf_process_noise_std, kf_measurement_noise_std, kf_dt);
            VectorXd state0 = VectorXd::Zero(kf.num_states_);
            state0.head(3) = pose_sol.pose.pos;
            state0.segment(3, 3) = pose_sol.pose.quat.toRotationMatrix().eulerAngles(0, 1, 2);
            MatrixXd covar0 = 10.0 * MatrixXd::Identity(kf.num_states_, kf.num_states_);
            covar0(0, 0) = 1.0;
            covar0(1, 1) = 1.0;
            covar0(2, 2) = 3.0;
            covar0(9, 9) = 10.0 * Utilities::DEG2RAD;
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

            pose_filtered.pos = pose_sol.pose.pos;
            pose_filtered.quat = pose_sol.pose.quat;

            //solved_poses.push_back( pose_sol.pose );
        }
        else // else, perform KF tracking
        {
            // prediction step
            kf.Predict(VectorXd::Zero(kf.num_inputs_));

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

            /*
            double      meas_att_score = Utilities::AttitudeScore(pose_true.quat, pose_sol.pose.quat);
            double conj_meas_att_score = Utilities::AttitudeScore(pose_true.quat, conj_pose.quat);
            if ( meas_att_score < conj_meas_att_score )
            { kf.Update(pose_meas_wrapper); solved_poses.push_back( pose_sol.pose ); }
            else
            { kf.Update(conj_pose_meas_wrapper); solved_poses.push_back( conj_pose ); }
            */

            kf.StoreAndClean();

            VectorXd pose_filt_wrapper = kf.last_state_estimate;
            pose_filtered.pos = pose_filt_wrapper.head(3);
            pose_filtered.quat = AngleAxisd(pose_filt_wrapper(9), Vector3d::UnitX()) *
                                 AngleAxisd(pose_filt_wrapper(10), Vector3d::UnitY()) *
                                 AngleAxisd(pose_filt_wrapper(11), Vector3d::UnitZ());
            //pose_filtered.quat = pose_filtered.quat.conjugate(); // TODO KEEP OR REMOVE ?
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
        solved_poses_conj.push_back(conj_pose);
        filtered_poses.push_back(pose_filtered);
        kf_states.push_back(kf.last_state_estimate);
        kf_covars.push_back(kf.last_covar_estimate);
        solution_times.push_back((double)duration);
        pos_scores.push_back(pos_score);
        att_scores.push_back(att_score); //std::min(att_score,conj_att_score) );

        if (first_run)
        {
            first_run = false;
        }
    }

    //-- Performance Metric Stats & Output -----------------------------------/

    // write to csv file
    Utilities::WritePosesToCSV(true_poses, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(true_poses))), false);
    Utilities::WritePosesToCSV(solved_poses, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(solved_poses))), false);
    Utilities::WritePosesToCSV(filtered_poses, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(filtered_poses))), false);
    Utilities::WriteKFStatesToCSV(kf_states, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(kf_states))), false);
    Utilities::WriteKFCovarsToCSV(kf_covars, Utilities::WrapVarToPath(std::string(GET_VARIABLE_NAME(kf_covars))), false);

    double pos_score_mean = Utilities::StdVectorMean(pos_scores);
    double att_score_mean = Utilities::StdVectorMean(att_scores);

    double pos_score_std = sqrt(Utilities::StdVectorVar(pos_scores));
    double att_score_std = sqrt(Utilities::StdVectorVar(att_scores));

    double solution_times_mean = Utilities::StdVectorMean(solution_times);

    std::cout << "num_runs :\t" << num_poses_test << std::endl
              << std::endl;

    std::cout << "pos_score_mean :\t" << pos_score_mean << /*" [m]" <<*/ std::endl;
    std::cout << "pos_score_std  :\t" << pos_score_std << /*" [m]" <<*/ std::endl
              << std::endl;

    std::cout << "att_score_mean :\t" << att_score_mean * Utilities::RAD2DEG << " [deg]" << std::endl;
    std::cout << "att_score_std  :\t" << att_score_std * Utilities::RAD2DEG << " [deg]" << std::endl
              << std::endl;

    std::cout << "mean_time  :\t" << solution_times_mean << " [ms]" << std::endl;
    std::cout << "total_time :\t" << solution_times_mean * num_poses_test << " [ms]" << std::endl;
    //------------------------------------------------------------------------/

    return 0;
}