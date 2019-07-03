#include "PoseSolver.h"

/**
 * @function SolvePose
 * @brief Non-linear Least-Squares Levenberg–Marquardt Solver for Pose based on
 *        Relative Bearing Measurements to Specified Feature Points with Known
 *        2D-3D Correspondences A-Priori
 * @return VectorXd of estimate state (pose)
 */
VectorXd PoseSolver::SolvePose(VectorXd yVec, VectorXd stateVec0, Vector3d rCamVec, MatrixXd rFeaMat)
{
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
    options.num_threads = 4;
    options.use_inner_iterations = false;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // convert estimated state information from double arrays to Eigen
    VectorXd stateHatVec(6);
    stateHatVec.head(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(posHatArr);
    stateHatVec.tail(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(eulHatArr);

    std::cout << summary.BriefReport() << "\n";
    
    return stateHatVec;
}