#include <Eigen/Core>
#include <math.h>

#include "ceres/ceres.h"
#include "glog/logging.h"

#include "cost_functor.h"

using Eigen::Vector3d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;


int main(int argc, char** argv)
{
    //google::InitGoogleLogging(argv[0]);

    // The variables to solve for with its initial value.
    // The veriables will be mutated in place by the solver.
    double posArr [3] = { 5.0, -3.0, 27.0 };
    double eulArr [3] = { 8.0, 3.0, -2.0 };

    double posArr0 [3];
    memcpy(posArr, posArr0, sizeof(posArr));

    double eulArr0 [3];
    memcpy(eulArr, eulArr0, sizeof(eulArr));

    // specify measurements
    VectorXd yVec(8);
    yVec <<    -0.0500,   -0.0045,   -0.0552,    0.0606,   -0.0425,   -0.0975,   -0.0983,   -0.0343;
    //yVec <<    -0.0594,   -0.0365,   -0.0158,    0.0455,   -0.0481,   -0.0746,   -0.0907,   -0.0403; // true meas


    // specify rigid position vector of feature points wrt target in target frame
    MatrixXd rFeaMat(4,3);
    rFeaMat <<  0.0,    0.0,    0.5,
                0.0,    0.0,   -1.5,
                0.0,    1.0,    1.0,
                0.0,   -1.0,    1.0;

    // specify rigid position vector of camera wrt chaser in chaser frame
    Vector3d rCamVec;
    rCamVec << 0.0, 0.0, 0.0;

    int numPts = rFeaMat.rows();

    // Build the problem.
    ceres::Problem problem;

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).
    //ceres::CostFunction* cost_function =
    //    MeasResidCostFunctor::Create(yVec, rFeaMat, rCamVec);
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<MeasResidCostFunctor, ceres::DYNAMIC, 3, 3>(
            new MeasResidCostFunctor(yVec, rFeaMat, rCamVec), numPts*2);
    
    problem.AddResidualBlock(cost_function, NULL, posArr, eulArr);

    // Run the solver
    ceres::Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    
    // True state: -1.0443, -0.3456, 21.4858, -61.7028, 83.3595, -133.3508;

    std::cout << "posVec :\t"; // << posVec0 << " -> " << posVec << "\n";
    for (const auto& e : posArr0) { std::cout << e << ", "; }
    std::cout << "\t->\t";
    for (const auto& e : posArr)  { std::cout << e << ", "; }
    std::cout << "[m]" << std::endl;

    std::cout << "eulVec :\t";
    for (const auto& e : eulArr0) { std::cout << e*180.0/M_PI << ", "; }
    std::cout << "\t->\t";
    for (const auto& e : eulArr)  { std::cout << e*180.0/M_PI << ", "; }
    std::cout << "[deg]" << std::endl;

    return 0;
}