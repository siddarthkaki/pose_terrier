/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <Eigen/Core>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <math.h>

#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>

#include "Utilities.h"

#include "third_party/json.hpp"

using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXd;
using nlohmann::json;

/**
 * @function main
 * @brief main function
 */
int main(int argc, char **argv)
{
    std::srand((unsigned int) time(NULL));

    //-- Read-in problem geometry and params ---------------------------------/

    // read params from JSON file
    std::ifstream input_stream("../params.json");
    json json_params;
    input_stream >> json_params;

    // specify rigid position vector of camera wrt chaser in chaser frame
    Vector3d rCamVec;
    for (unsigned int idx = 0; idx < 3; idx++)
    {
        rCamVec(idx) = json_params["rCamVec"].at(idx);
    }

    // specify camera focal length
    double focal_length = json_params["focal_length"]; //5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    double meas_std = double(json_params["meas_std_deg"]) * Utilities::DEG2RAD;

    // specify rigid position vector of feature points wrt target in target frame
    unsigned int num_features = json_params["rFeaMat"].size();
    MatrixXd rFeaMat(num_features,3);
    for (unsigned int idx = 0; idx < num_features; idx++)
    {   for (unsigned int jdx = 0; jdx < 3; jdx++)
        { rFeaMat(idx,jdx) = json_params["rFeaMat"][idx]["fea" + std::to_string(idx+1)][jdx]; } }

    /*
    // TEMPORARY
    int num_features = 11;
    rFeaMat = 2.5 * MatrixXd::Random(num_features, 3);
    std::cout << rFeaMat << std::endl << std::endl;
    */

    unsigned int num_poses_test = json_params["num_poses_test"];

    //------------------------------------------------------------------------/

    // 3D model points
    std::vector<cv::Point3d> model_points;
    for (unsigned int idx = 0; idx < num_features; idx++)
    {
        model_points.push_back(cv::Point3d(rFeaMat(idx, 0), rFeaMat(idx, 1), rFeaMat(idx, 2)));
    }

    //-- Loop ----------------------------------------------------------------/

    std::vector<Pose> solved_poses;
    std::vector<Pose> solved_poses_conj;
    std::vector<double> solution_times; // [ms]
    std::vector<double> pos_scores;
    std::vector<double> att_scores;

    solved_poses.reserve(num_poses_test);
    solved_poses_conj.reserve(num_poses_test);
    solution_times.reserve(num_poses_test);
    pos_scores.reserve(num_poses_test);
    att_scores.reserve(num_poses_test);

    for (unsigned int pose_idx = 0; pose_idx < num_poses_test; pose_idx++)
    {

        //-- Simulate Measurements -------------------------------------------/

        // generate true state values for ith run
        Pose poseTrue;
        poseTrue.pos << 0.0, 0.0, 25.0;
        poseTrue.pos.head(2) = Utilities::AddGaussianNoiseToVector(poseTrue.pos.head(2), 1);
        poseTrue.pos.tail(1) = Utilities::AddGaussianNoiseToVector(poseTrue.pos, 3).tail(1);
        poseTrue.quat = Quaterniond::UnitRandom();

        //std::cout << poseTrue.pos.transpose() << "\t\t" << poseTrue.quat.w() << " " << poseTrue.quat.vec().transpose() << std::endl;

        // express feature points in chaser frame at the specified pose
        MatrixXd rMat = Utilities::FeaPointsTargetToChaser(poseTrue, rCamVec, rFeaMat);

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

        double f = focal_length;
        MatrixXd PMat(3,3);
        PMat << f, 0, 0,
                0, f, 0,
                0, 0, 1;
        
        cv::Mat camera_matrix;
        cv::eigen2cv(PMat, camera_matrix);

        cv::Mat dist_coeffs = cv::Mat::zeros(4,1,cv::DataType<double>::type); // Assuming no lens distortion

        // timing
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        // Output rotation and translation
        cv::Mat rotation_vector; // Rotation in axis-angle form
        cv::Mat translation_vector;
        
        // Solve for pose
        cv::solvePnPRansac(model_points, image_points, camera_matrix, dist_coeffs, rotation_vector, translation_vector);
    
        Pose pose;
        cv::cv2eigen(translation_vector, pose.pos);
        cv::Mat R;
        cv::Rodrigues(rotation_vector, R); // R is 3x3
        Eigen::Matrix3d mat;
        cv::cv2eigen(R, mat);
        Eigen::Quaterniond EigenQuat(mat);
        pose.quat = EigenQuat;

        Pose conj_pose = Utilities::ConjugatePose(pose);

        //std::cout << pose.pos.transpose() << "\t\t" << pose.quat.w() << " " << pose.quat.vec().transpose() << std::endl;
        //std::cout << conj_pose.pos.transpose() << "\t\t" << conj_pose.quat.w() << " " << conj_pose.quat.vec().transpose() << std::endl << std::endl;

        // timing
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

        // time taken to perform NLS solution
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1000000.0;
        //--------------------------------------------------------------------/

        //-- Performance Metrics & Storage -----------------------------------/

        // compute position and attitude scores
        double pos_score = Utilities::PositionScore(poseTrue.pos, pose.pos);
        double att_score = Utilities::AttitudeScore(poseTrue.quat, pose.quat);
        double conj_att_score = Utilities::AttitudeScore(poseTrue.quat, conj_pose.quat);

        // store info from ith run
        solved_poses.push_back(pose);
        solved_poses_conj.push_back(conj_pose);
        solution_times.push_back((double)duration);
        pos_scores.push_back(pos_score);
        //att_scores.push_back(att_score);
        att_scores.push_back(std::min(att_score, conj_att_score));
    }

    //-- Performance Metric Stats & Output -----------------------------------/

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