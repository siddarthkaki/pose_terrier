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

#include "cost_functor.h"
#include "Utilities.h"
#include "PoseSolver.h"

#include "third_party/json.hpp"

using Eigen::Vector3d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;
using Eigen::Quaterniond;
using nlohmann::json;

/**
 * @function main
 * @brief main function
 */
int main(int argc, char** argv)
{
    //google::InitGoogleLogging(argv[0]);

    //-- Read-in problem geometry and params ---------------------------------/

    // read params from JSON file
    std::ifstream input_stream("../params.json");
    json json_params;
    input_stream >> json_params;

    // specify rigid position vector of camera wrt chaser in chaser frame
    Vector3d rCamVec;
    for (unsigned int idx = 0; idx < 2; idx++)
    { rCamVec(idx) = json_params["rCamVec"].at(idx); }

    // specify camera focal length
    double focal_length = json_params["focal_length"];//5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    double meas_std = double(json_params["meas_std_deg"])*Utilities::DEG2RAD;

    // specify rigid position vector of feature points wrt target in target frame
    unsigned int num_features = json_params["rFeaMat"].size();
    MatrixXd rFeaMat(num_features,3);
    for (unsigned int idx = 0; idx < num_features; idx++)
    {   for (unsigned int jdx = 0; jdx < 3; jdx++)
        { rFeaMat(idx,jdx) = json_params["rFeaMat"][idx]["fea" + std::to_string(idx+1)][jdx]; } }
    
    //------------------------------------------------------------------------/

    // initial state guess
    double posArr0[3] = { 0.0, 0.0, 25.0 };
    double quatArr0[4] = {1.0, 0.0, 0.0, 0.0}; // w,x,y,z

    // convert initial state information from double arrays to Eigen
    Pose pose0;
    pose0.pos(0) = posArr0[0];
    pose0.pos(1) = posArr0[1];
    pose0.pos(2) = posArr0[2];
    pose0.quat.w() = quatArr0[0];
    pose0.quat.x() = quatArr0[1];
    pose0.quat.y() = quatArr0[2];
    pose0.quat.z() = quatArr0[3]; 
    pose0.quat.normalize();

    //-- Simulate Measurements -----------------------------------------------/

    // true state information
    double posArr [3] = {0.5377, 1.8339, 18.2235};
    //double quatArr [4] = {0.6937, -0.6773, 0.0642, 0.2365};
    double quatArr [4] = {1.2340,   -1.5971,    0.7174,    -0.2721};
    
    // convert true state information from double arrays to Eigen
    Pose poseTrue;
    poseTrue.pos(0) = posArr[0];
    poseTrue.pos(1) = posArr[1];
    poseTrue.pos(2) = posArr[2];
    poseTrue.quat.w() = quatArr[0];
    poseTrue.quat.x() = quatArr[1];
    poseTrue.quat.y() = quatArr[2];
    poseTrue.quat.z() = quatArr[3]; 
    poseTrue.quat.normalize();

    // express feature points in chaser frame at the specified pose
    MatrixXd rMat = Utilities::FeaPointsTargetToChaser(poseTrue, rCamVec, rFeaMat);

    // generate simulated measurements
    VectorXd yVec = Utilities::SimulateMeasurements(rMat, focal_length);

    // add Gaussian noise to simulated measurements
    VectorXd yVecNoise = Utilities::AddGaussianNoiseToVector(yVec, meas_std);

    //------------------------------------------------------------------------/

    //-- Solve for pose ------------------------------------------------------/

    // timing
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // solve for pose with ceres (via wrapper)
    PoseSolution poseSol = PoseSolver::SolvePoseReinit(pose0, yVecNoise, rCamVec, rFeaMat);

    Pose conj_pose_temp = Utilities::ConjugatePose(poseSol.pose);
    Pose conj_pose = PoseSolver::SolvePose(conj_pose_temp, yVecNoise, rCamVec, rFeaMat).pose;

    // timing
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    
    // time taken to perform NLS solution
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    //------------------------------------------------------------------------/

    //-- Performance Metrics & Outputs ---------------------------------------/

    // compute position and attitude scores
    double      pos_score = Utilities::PositionScore(poseTrue.pos , poseSol.pose.pos );
    double      att_score = Utilities::AttitudeScore(poseTrue.quat, poseSol.pose.quat);
    double conj_att_score = Utilities::AttitudeScore(poseTrue.quat,    conj_pose.quat);

    // print to command line
    std::cout << poseSol.summary.BriefReport() << "\n";
    //std::cout << poseSol.summary.FullReport() << "\n";
    
    
    std::cout << "posVec :\t" << pose0.pos.transpose();
    std::cout << "\t->\t";
    std::cout << poseSol.pose.pos.transpose() << " [m]" << std::endl;

    std::cout << "attVec :\t" << pose0.quat.w() << " " << pose0.quat.vec().transpose();
    std::cout << "\t\t->\t";
    std::cout << poseSol.pose.quat.w() << " " << poseSol.pose.quat.vec().transpose() << std::endl;
    std::cout << "\t\t\t\t\t" << conj_pose.quat.w() << " " << conj_pose.quat.vec().transpose() << std::endl;
    

    std::cout << "pos_score :\t\t" << pos_score << " [m]" << std::endl;
    std::cout << "att_score :\t\t" << att_score*Utilities::RAD2DEG << " [deg]"<< std::endl;
    std::cout << "conj_att_score :\t" << conj_att_score*Utilities::RAD2DEG << " [deg]"<< std::endl;

    std::cout << "Time taken by program is : "  << std::setprecision(9)
                                                << (double)duration
                                                << " [ms]"
                                                << std::endl; 

    //------------------------------------------------------------------------/

    return 0;
}