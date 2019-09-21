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

#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>

#include "cost_functor.h"
#include "mekf_f_functor.h"
#include "Utilities.h"
#include "KalmanFilter.h"
#include "MEKF.h"

#include "third_party/json.hpp"
#include "third_party/SmartBuffer.h"

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
    
    //num_features = 11;
    //rFeaMat = 2.5 * MatrixXd::Random(num_features, 3);

    unsigned int num_poses_test = json_params["num_poses_test"];

    double kf_process_noise_std = 0.01;
    double kf_measurement_noise_std = 0.05;
    double kf_dt = 0.01;

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

    // KF object
    //KalmanFilter::KalmanFilter kf;

    // MEKF object
    MEKF::MEKF mekf(kf_dt);

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

    // set-up for Jacobian computation with ceres
    //constexpr unsigned int num_states_F = 9; // eul (3), omega (3), alpha (3)
    constexpr unsigned int num_states_F = 6; // eul (3), omega (3)

    bool first_run = true;

    for (unsigned int pose_idx = 0; pose_idx < num_poses_test; pose_idx++)
    {
        //-- Simulate Measurements -------------------------------------------/

        // generate true pose values for ith run
        pose_true.pos(0) += 0.01;
        pose_true.pos(1) -= 0.01;
        pose_true.pos(2) += 0.05;
        Quaterniond quat_step = AngleAxisd( 0.001, Vector3d::UnitX()) *
                                AngleAxisd(-0.001, Vector3d::UnitY()) *
                                AngleAxisd( 0.001, Vector3d::UnitZ());
        pose_true.quat = pose_true.quat * quat_step;

        // express feature points in chaser frame at the specified pose
        MatrixXd rMat = Utilities::FeaPointsTargetToChaser(pose_true, rCamVec, rFeaMat);

        // generate simulated measurements
        VectorXd yVec = Utilities::SimulateMeasurements(rMat, focal_length);

        // add Gaussian noise to simulated measurements
        VectorXd yVecNoise = Utilities::AddGaussianNoiseToVector(yVec, meas_std);

        // 2D image points
        std::vector<cv::Point2d> image_points;

        // project feature points to image plane
        for (unsigned int idx = 0; idx < num_features; idx++)
        {
            double az = yVecNoise(idx * 2 + 0);
            double el = yVecNoise(idx * 2 + 1);

            image_points.push_back(cv::Point2d(tan(az)*focal_length, tan(el)*focal_length));
        }

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
            mekf.Init(kf_process_noise_std, kf_measurement_noise_std, kf_dt);
            VectorXd att_state0 = VectorXd::Zero(mekf.num_states_);
            //state0.head(3) = pose_sol.pose.pos;
            att_state0(0) = pose_sol.pose.quat.w();
            att_state0(1) = pose_sol.pose.quat.x();
            att_state0(2) = pose_sol.pose.quat.y();
            att_state0(3) = pose_sol.pose.quat.z();
            att_state0.tail(3) = Vector3d::Random();
            MatrixXd att_covar0 = 1.0 * MatrixXd::Identity(mekf.num_states_covar_, mekf.num_states_covar_);
            att_covar0(0,0) = 0.1;
            att_covar0(1,1) = 0.1;
            att_covar0(2,2) = 0.1;
            /*
            covar0(0, 0) = 1.0;
            covar0(1, 1) = 1.0;
            covar0(2, 2) = 3.0;
            covar0(9, 9) = 10.0 * Utilities::DEG2RAD;
            covar0(10, 10) = 10.0 * Utilities::DEG2RAD;
            covar0(11, 11) = 10.0 * Utilities::DEG2RAD;
            */
            mekf.SetInitialStateAndCovar(att_state0, att_covar0);
            
            // NOTE: F not accurate in print-out
            mekf.PrintModelMatrices();

            pose_filtered.pos = pose_sol.pose.pos;
            pose_filtered.quat = pose_sol.pose.quat;

        }
        else // else, perform KF tracking
        {
            // set-up for Jacobian computation with ceres
            ceres::CostFunction *mekf_f_wrapper = new ceres::AutoDiffCostFunction<MEKF_f_Functor, num_states_F, num_states_F>(
                new MEKF_f_Functor(mekf.dt_));

            // set-up for computing and storing F Jacobian matrix for current time-step
            
            // prepare statekk_ accessor; need pointer to pointer
            //Eigen::Matrix<double, 9, 1> parameters_state_buffer; // eul (3), omega (3), alpha (3)
            Eigen::Matrix<double, 6, 1> parameters_state_buffer; // eul (3), omega (3)
            Quaterniond quatkk_buffer;
            quatkk_buffer.w() = mekf.statekk_(0);
            quatkk_buffer.x() = mekf.statekk_(1);
            quatkk_buffer.y() = mekf.statekk_(2);
            quatkk_buffer.z() = mekf.statekk_(3);
            Vector3d eulkk_buffer = quatkk_buffer.toRotationMatrix().eulerAngles(0, 1, 2);
            parameters_state_buffer.head(3) = eulkk_buffer;
            //parameters_state_buffer.tail(6) = mekf.statekk_.tail(6);
            parameters_state_buffer.tail(3) = mekf.statekk_.tail(3);
            const double *parameters = parameters_state_buffer.data();

            // structures for autodiff evaluation
            SmartBuffer1D residuals(num_states_F);
            SmartBuffer2D jacobian(1, num_states_F*num_states_F);

            // Evaluate jacobian
            bool success = mekf_f_wrapper->Evaluate(&parameters, residuals.Get(), jacobian.Get());

            if (success == true)
            {
                mekf.F_ = Utilities::ConvertToEigenMatrix(jacobian.Get(), num_states_F, num_states_F);
            }
            else
            {
                std::cout << "Jacobian computation failed!" << std::endl;
            }

            // prediction step (state propagation in terms of quaternions, covariance propagation in terms of euler angles)
            mekf.Predict(VectorXd::Zero(mekf.num_inputs_), MEKF_f_Functor::MEKF_f_quat);

            // wrap NLS pose solution as KF measurement
            VectorXd att_meas_wrapper(4);
            att_meas_wrapper(0) = pose_sol.pose.quat.normalized().w();
            att_meas_wrapper(1) = pose_sol.pose.quat.normalized().x();
            att_meas_wrapper(2) = pose_sol.pose.quat.normalized().y();
            att_meas_wrapper(3) = pose_sol.pose.quat.normalized().z();

            // wrap NLS conjugate pose solution as KF measurement
            VectorXd conj_att_meas_wrapper(4);
            conj_att_meas_wrapper(0) = conj_pose.quat.normalized().w();
            conj_att_meas_wrapper(1) = conj_pose.quat.normalized().x();
            conj_att_meas_wrapper(2) = conj_pose.quat.normalized().y();
            conj_att_meas_wrapper(3) = conj_pose.quat.normalized().z();

            
            // choose as measurement whichever attitude produces the smallest measurement residual norm
            //double pose_meas_norm = (pose_meas_wrapper - mekf.H_ * mekf.statek1k_).norm();
            //double conj_pose_meas_norm = (conj_pose_meas_wrapper - mekf.H_ * mekf.statek1k_).norm();
            //if (pose_meas_norm < conj_pose_meas_norm)
            //{
            mekf.Update(att_meas_wrapper, MEKF::MEKF::MeasResidFunction);
            //}
            //else
            {
                //mekf.Update(conj_pose_meas_wrapper);
            }

            mekf.Reset();

            mekf.StoreAndClean();

            VectorXd att_filt_wrapper = mekf.last_state_estimate;
            pose_filtered.pos = Vector3d::Zero();

            pose_filtered.quat.w() = att_filt_wrapper(0);
            pose_filtered.quat.x() = att_filt_wrapper(1);
            pose_filtered.quat.y() = att_filt_wrapper(2);
            pose_filtered.quat.z() = att_filt_wrapper(3);
            pose_filtered.quat.normalize();

            // delete pointers after use
            /*
            for (unsigned int idx = 0; idx < num_states; idx++)
            {
                delete jacobians[idx];
            }
            delete parameters[0];
            delete jacobians;
            delete parameters;
            delete nonlinear_propagation_auto_diff_wrapper;
            */
           delete mekf_f_wrapper;
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
        kf_states.push_back(mekf.last_state_estimate);
        kf_covars.push_back(mekf.last_covar_estimate);
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
    Utilities::WritePosesToCSV(true_poses, Utilities::WrapVarToPath(std::string(GET_VAR_NAME(true_poses))), false);
    Utilities::WritePosesToCSV(solved_poses, Utilities::WrapVarToPath(std::string(GET_VAR_NAME(solved_poses))), false);
    Utilities::WritePosesToCSV(filtered_poses, Utilities::WrapVarToPath(std::string(GET_VAR_NAME(filtered_poses))), false);
    Utilities::WriteKFStatesToCSV(kf_states, Utilities::WrapVarToPath(std::string(GET_VAR_NAME(kf_states))), false);
    Utilities::WriteKFCovarsToCSV(kf_covars, Utilities::WrapVarToPath(std::string(GET_VAR_NAME(kf_covars))), false);

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
    //------------------------------------------------------------------------/

    return 0;
}