#include <Eigen/Core>
#include <math.h>

#include "ceres/ceres.h"
//#include "glog/logging.h"

#include "cost_functor.h"
#include "Utilities.h"

using Eigen::Vector3d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;

// Functions

/**
 * @function main
 * @brief main function
 */
int main(int argc, char** argv)
{
    //google::InitGoogleLogging(argv[0]);

    const double DEG2RAD = M_PI/180.0;
    const double RAD2DEG = 180.0/M_PI;

    // true state information
    //double posArr [3] = { -1.0443, -0.3456, 21.4858 };
    //double eulArr [3] = { -61.7028*DEG2RAD, 83.3595*DEG2RAD, -133.3508*DEG2RAD };

    double posArr [3] = { 0.5377, 1.8339, 18.2235 };
    double eulArr [3] = { 1.3543, 0.5007, -2.0541 };

    // initial state guess
    //double posArr0 [3] = { 6.2, -4.6, 30.0 };
    //double eulArr0 [3] = { -60.0*DEG2RAD, 80.0*DEG2RAD, -135.0*DEG2RAD };

    double posArr0 [3] = { -1.0, 3.0, 25.0 };
    double eulArr0 [3] = { 1.5, 0.4, -2.2 };

    // The variables to solve for with its initial value.
    // The veriables will be mutated in place by the solver.
    double posArrHat [3];
    double eulArrHat [3];
    
    memcpy(posArrHat, posArr0, sizeof(posArr));
    memcpy(eulArrHat, eulArr0, sizeof(eulArr));

    // TODO NOT WORKING!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // convert true state information from double arrays to Eigen
    VectorXd stateVec(6);
    stateVec.head(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(posArr);
    stateVec.tail(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(eulArr);

    // specify rigid position vector of camera wrt chaser in chaser frame
    Vector3d rCamVec;
    rCamVec << 0.0, 0.0, 0.0;

    // specify camera focal length
    double focal_length = 5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    double meas_std = 1.0*DEG2RAD;

    // specify rigid position vector of feature points wrt target in target frame
    MatrixXd rFeaMat(4,3);
    rFeaMat <<  0.0,    0.0,    0.5,
                0.0,    0.0,   -1.5,
                0.0,    1.0,    1.0,
                0.0,   -1.0,    1.0;

    // express feature points in chaser frame at the specified pose
    MatrixXd rMat = Utilities::FeaPointsTargetToChaser(stateVec, rCamVec, rFeaMat);

    // generate simulated measurements
    VectorXd yVec = Utilities::SimulateMeasurements(rMat, focal_length);

    // add Gaussian noise to simulated measurements
    VectorXd yVecNoise = Utilities::AddNoiseToMeasurements(yVec, meas_std);

    // specify true measurements
    //VectorXd yVec(8);
    //yVec <<    -0.0500,   -0.0045,   -0.0552,    0.0606,   -0.0425,   -0.0975,   -0.0983,   -0.0343;
    //yVec <<    -0.0594,   -0.0365,   -0.0158,    0.0455,   -0.0481,   -0.0746,   -0.0907,   -0.0403; // true meas
    

    int numPts = rFeaMat.rows();

    // Build the problem.
    ceres::Problem problem;

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).
    //ceres::CostFunction* cost_function =
    //    MeasResidCostFunctor::Create(yVec, rFeaMat, rCamVec);
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<MeasResidCostFunctor, ceres::DYNAMIC, 3, 3>(
            new MeasResidCostFunctor(yVec, rFeaMat, rCamVec), numPts*2);
    
    problem.AddResidualBlock(cost_function, NULL, posArrHat, eulArrHat);

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
    for (const auto& e : posArrHat)  { std::cout << e << ", "; }
    std::cout << "[m]" << std::endl;

    std::cout << "eulVec :\t";
    for (const auto& e : eulArr0) { std::cout << e*180.0/M_PI << ", "; }
    std::cout << "\t->\t";
    for (const auto& e : eulArrHat)  { std::cout << e*180.0/M_PI << ", "; }
    std::cout << "[deg]" << std::endl;

    return 0;
}