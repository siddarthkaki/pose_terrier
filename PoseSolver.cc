/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include "PoseSolver.h"

/**
 * @function SolvePose
 * @brief Non-linear Least-Squares Levenberg–Marquardt Solver for Pose based on
 *        Relative Bearing Measurements to Specified Feature Points with A-Priori 
 *        Known 2D-3D Correspondences 
 * @return VectorXd of estimate state (pose)
 */
PoseSolution PoseSolver::SolvePose(const Pose& pose0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat, const double bearing_meas_std)
{
    PoseSolution poseSol;

    VectorXd quatHatVec(4); 
    quatHatVec << pose0.quat.w(), pose0.quat.x(), pose0.quat.y(), pose0.quat.z();
    Vector3d  posHatVec = pose0.pos;
    
    // The variables to solve for with initial values.
    // The variables will be mutated in place by the solver.
    double* quatHatArr = quatHatVec.data();
    double*  posHatArr =  posHatVec.data();
    
    // number of feature points
    int numPts = rFeaMat.rows();

    // Build the problem.
    ceres::Problem problem;

    // Specify parameterisation for quaternion block.
    ceres::LocalParameterization* quaternion_parameterization = new ceres::QuaternionParameterization;

    // Specify loss function.
    ceres::LossFunction* loss_function = NULL; // new ceres::HuberLoss(0.1);

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).    
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<MeasResidCostFunctorQuat, ceres::DYNAMIC, 4, 3>(
            new MeasResidCostFunctorQuat(yVec, rFeaMat, rCamVec, bearing_meas_std), numPts*2);

    problem.AddResidualBlock(cost_function, loss_function, quatHatArr, posHatArr);

    problem.SetParameterization(quatHatArr, quaternion_parameterization);

    // Run the solver
    ceres::Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.linear_solver_type = ceres::DENSE_QR;
    //options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
    options.use_nonmonotonic_steps = true;
    options.num_threads = 1;
    options.use_inner_iterations = false;
    options.minimizer_progress_to_stdout = false;
    ceres::Solve(options, &problem, &poseSol.summary);

    // Retrieve parameter blocks
    // std::cout << problem.NumParameterBlocks() << std::endl;
    // std::vector<double*> parameter_blocks;
    // problem.GetParameterBlocks(&parameter_blocks);
    // double* posParamBlock = parameter_blocks.at(1);
    // double* quatParamBlock = parameter_blocks.at(0);

    // Covariance computation
    ceres::Covariance::Options cov_options;
    cov_options.algorithm_type = ceres::DENSE_SVD;
    //cov_options.null_space_rank = 1;
    ceres::Covariance covariance(cov_options);

    std::vector<std::pair<const double*, const double*> > covariance_blocks;    
    covariance_blocks.push_back(std::make_pair(quatHatArr, quatHatArr));
    covariance_blocks.push_back(std::make_pair(posHatArr, posHatArr));
    covariance_blocks.push_back(std::make_pair(quatHatArr, posHatArr));

    CHECK(covariance.Compute(covariance_blocks, &problem));
    Matrix3d_rm cov_xx = Matrix3d_rm::Zero();
    Matrix3d_rm cov_yy = Matrix3d_rm::Zero();
    Matrix3d_rm cov_xy = Matrix3d_rm::Zero();
    covariance.GetCovarianceBlockInTangentSpace(quatHatArr, quatHatArr, cov_xx.data());
    covariance.GetCovarianceBlock(posHatArr, posHatArr, cov_yy.data());
    covariance.GetCovarianceBlockInTangentSpace(quatHatArr, posHatArr, cov_xy.data());
    
    // convert estimated state information from double arrays to Eigen
    poseSol.pose.quat.w() = quatHatVec(0);
    poseSol.pose.quat.x() = quatHatVec(1);
    poseSol.pose.quat.y() = quatHatVec(2);
    poseSol.pose.quat.z() = quatHatVec(3);
    poseSol.pose.pos  = posHatVec;
    
    // store pose covariance estimate
    poseSol.cov_pose.topLeftCorner(3,3) = cov_xx;
    poseSol.cov_pose.topRightCorner(3,3) = cov_xy;
    poseSol.cov_pose.bottomLeftCorner(3,3) = cov_xy.transpose();
    poseSol.cov_pose.bottomRightCorner(3,3) = cov_yy;

    return poseSol;
}

/**
 * @function SolvePoseReinit
 * @brief Non-linear Least-Squares Levenberg–Marquardt Solver with Multiple
 *        Random Reinitialisations for Pose based on Relative Bearing
 *        Measurements to Specified Feature Points with A-Priori Known 2D-3D
 *        Correspondences 
 * @return VectorXd of estimate state (pose)
 */
PoseSolution PoseSolver::SolvePoseReinit(const Pose& pose0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat, const double bearing_meas_std)
{
    unsigned int num_init = 5;

    PoseSolution posSolOptimal;
    double min_cost = 100;

    for (unsigned int init_idx = 0; init_idx < num_init; init_idx++)
    {
        Pose pose0i = pose0;
        if (init_idx > 0)
        { pose0i.quat = Quaterniond::UnitRandom(); }
        PoseSolution posSoli = SolvePose(pose0i, yVec, rCamVec, rFeaMat, bearing_meas_std);   
    
        double curr_cost = posSoli.summary.final_cost;
        if (curr_cost < min_cost)
        {
            min_cost = curr_cost;
            posSolOptimal = posSoli;
        }
    }

    return posSolOptimal;
}

/**
 * @function TwoPointDiffTwistEstimator
 * @brief Estimates Twist (velocity + angular velocity) from Pose at current
 *        time j and Pose at time i = j - T, where T is the sampling period 
 * @return Twist struct
 */
Twist PoseSolver::TwoPointDiffTwistEstimator(const Pose& posei, const Pose& posej, const double& T)
{
    Twist twist_j;

    twist_j.vel =  (posej.pos - posei.pos) / T;

    Quaterniond dqdt = (posei.quat.conjugate()*posej.quat);
    dqdt.w() /= T;
    dqdt.vec() /= T;

    twist_j.ang_vel = 2*( dqdt*posej.quat.conjugate() ).vec();

    return twist_j;
}