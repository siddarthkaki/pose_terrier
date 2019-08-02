/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include "PoseSolver.h"

/**
 * @function SolvePose
 * @brief Non-linear Least-Squares Levenberg–Marquardt Solver for Pose based on
 *        Relative Bearing Measurements to Specified Feature Points with Known
 *        2D-3D Correspondences A-Priori
 * @return VectorXd of estimate state (pose)
 */
PoseSolution PoseSolver::SolvePose(const Pose& pose0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat)
{
    PoseSolution poseSol;

    Vector3d  posHatVec = pose0.pos;
    VectorXd quatHatVec(4); 
    quatHatVec << pose0.quat.w(), pose0.quat.x(), pose0.quat.y(), pose0.quat.z();

    // The variables to solve for with initial values.
    // The variables will be mutated in place by the solver.
    double*  posHatArr =  posHatVec.data();
    double* quatHatArr = quatHatVec.data();
    
    // number of feature points
    int numPts = rFeaMat.rows();

    // Build the problem.
    ceres::Problem problem;

    // Specify parameterisation for quaternion block.
    ceres::LocalParameterization *quaternion_parameterization = new ceres::QuaternionParameterization;

    // Specify loss function.
    ceres::LossFunction* loss_function = NULL; //new ceres::HuberLoss(0.1);

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).    
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<MeasResidCostFunctorQuat, ceres::DYNAMIC, 3, 4>(
            new MeasResidCostFunctorQuat(yVec, rFeaMat, rCamVec), numPts*2);

    problem.AddResidualBlock(cost_function, loss_function, posHatArr, quatHatArr);

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

    // convert estimated state information from double arrays to Eigen
    poseSol.pose.pos  = posHatVec;
    poseSol.pose.quat.w() = quatHatVec(0);
    poseSol.pose.quat.x() = quatHatVec(1);
    poseSol.pose.quat.y() = quatHatVec(2);
    poseSol.pose.quat.z() = quatHatVec(3);

    return poseSol;
}

/**
 * @function SolvePoseReinit
 * @brief Non-linear Least-Squares Levenberg–Marquardt Solver with Multiple
 *        Random Reinitialisationsfor Pose based on Relative Bearing
 *        Measurements to Specified Feature Points with A-Priori Known 2D-3D
 *        Correspondences 
 * @return VectorXd of estimate state (pose)
 */
PoseSolution PoseSolver::SolvePoseReinit(const Pose& pose0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat)
{
    unsigned int num_init = 5;

    PoseSolution posSolOptimal;
    double min_cost = 100;

    for (unsigned int init_idx = 0; init_idx < num_init; init_idx++)
    {
        Pose pose0i = pose0;
        if (init_idx > 0)
        { pose0i.quat = Quaterniond::UnitRandom(); }
        PoseSolution posSoli = SolvePose(pose0i, yVec, rCamVec, rFeaMat);   
    
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
}
