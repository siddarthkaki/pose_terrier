#include "PoseSolver.h"

/**
 * @function SolvePose
 * @brief Non-linear Least-Squares Levenberg–Marquardt Solver for Pose based on
 *        Relative Bearing Measurements to Specified Feature Points with Known
 *        2D-3D Correspondences A-Priori
 * @return VectorXd of estimate state (pose)
 */
PoseSolution PoseSolver::SolvePose(VectorXd stateVec0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat)
{
    PoseSolution poseSol;

    // The variables to solve for with initial values.
    // The variables will be mutated in place by the solver.
    double* posHatArr = stateVec0.head(3).data();
    double* eulHatArr = stateVec0.tail(3).data();
    
    //memcpy(posHatArr, posArr0, sizeof(posArr0));
    //memcpy(eulHatArr, eulArr0, sizeof(eulArr0));

    // number of feature points
    int numPts = rFeaMat.rows();

    // Build the problem.
    ceres::Problem problem;

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).
    //ceres::CostFunction* cost_function =
    //    MeasResidCostFunctor::Create(yVec, rFeaMat, rCamVec);
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<MeasResidCostFunctor, ceres::DYNAMIC, 3, 3>(
            new MeasResidCostFunctor(yVec, rFeaMat, rCamVec), numPts*2);
    
    problem.AddResidualBlock(cost_function, NULL, posHatArr, eulHatArr);

    // Run the solver
    ceres::Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.linear_solver_type = ceres::DENSE_QR;
    options.use_nonmonotonic_steps = true;
    options.num_threads = 1;
    options.use_inner_iterations = false;
    options.minimizer_progress_to_stdout = false;
    ceres::Solve(options, &problem, &poseSol.summary);

    // convert estimated state information from double arrays to Eigen
    poseSol.stateHatVec = VectorXd::Zero(6);
    poseSol.stateHatVec.head(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(posHatArr);
    poseSol.stateHatVec.tail(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(eulHatArr);

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
PoseSolution PoseSolver::SolvePoseReinit(const VectorXd& stateVec0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat)
{
    unsigned int num_init = 5;
    double reinit_att_noise_std = 2;

    PoseSolution posSolOptimal;
    double min_cost = 100;

    for (unsigned int init_idx = 0; init_idx < num_init; init_idx++)
    {
        VectorXd stateVec0i = stateVec0;
        stateVec0i.tail(3) = Utilities::AddGaussianNoiseToVector(stateVec0i.tail(3), reinit_att_noise_std);
        PoseSolution posSoli = SolvePose(stateVec0i, yVec, rCamVec, rFeaMat);   
    
        double curr_cost = posSoli.summary.final_cost;
        if (curr_cost < min_cost)
        {
            min_cost = curr_cost;
            posSolOptimal = posSoli;
        }
    }

    return posSolOptimal;
}