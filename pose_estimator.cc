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
    double posArr [3] = { 0.5377, 1.8339, 18.2235 };
    double eulArr [3] = { 1.3543, 0.5007, -2.0541 };

    // initial state guess
    double posArr0 [3] = { -1.0, 3.0, 25.0 };
    double eulArr0 [3] = { 1.5, 0.4, -2.2 };

    // The variables to solve for with its initial value.
    // The variables will be mutated in place by the solver.
    double posHatArr [3];
    double eulHatArr [3];
    
    memcpy(posHatArr, posArr0, sizeof(posArr));
    memcpy(eulHatArr, eulArr0, sizeof(eulArr));

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
    int numPts = rFeaMat.rows();

    // express feature points in chaser frame at the specified pose
    MatrixXd rMat = Utilities::FeaPointsTargetToChaser(stateVec, rCamVec, rFeaMat);

    // generate simulated measurements
    VectorXd yVec = Utilities::SimulateMeasurements(rMat, focal_length);

    // add Gaussian noise to simulated measurements
    VectorXd yVecNoise = Utilities::AddNoiseToMeasurements(yVec, meas_std);

    // Build the problem.
    ceres::Problem problem;

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).
    //ceres::CostFunction* cost_function =
    //    MeasResidCostFunctor::Create(yVecNoise, rFeaMat, rCamVec);
    ceres::CostFunction* cost_function = new ceres::AutoDiffCostFunction<MeasResidCostFunctor, ceres::DYNAMIC, 3, 3>(
            new MeasResidCostFunctor(yVecNoise, rFeaMat, rCamVec), numPts*2);
    
    problem.AddResidualBlock(cost_function, NULL, posHatArr, eulHatArr);

    // Run the solver
    ceres::Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.minimizer_progress_to_stdout = true;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    // convert estimated state information from double arrays to Eigen
    VectorXd stateHatVec(6);
    stateHatVec.head(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(posHatArr);
    stateHatVec.tail(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(eulHatArr);

    // compute position and attitude scores
    double pos_score = Utilities::PositionScore(stateVec, stateHatVec);
    double att_score = Utilities::AttitudeScore(stateVec, stateHatVec);

    // print to command line
    std::cout << summary.BriefReport() << "\n";
    
    std::cout << "posVec :\t"; // << posVec0 << " -> " << posVec << "\n";
    for (const auto& e : posArr0) { std::cout << e << ", "; }
    std::cout << "\t->\t";
    for (const auto& e : posHatArr)  { std::cout << e << ", "; }
    std::cout << "[m]" << std::endl;

    std::cout << "eulVec :\t";
    for (const auto& e : eulArr0) { std::cout << e*180.0/M_PI << ", "; }
    std::cout << "\t->\t";
    for (const auto& e : eulHatArr)  { std::cout << e*180.0/M_PI << ", "; }
    std::cout << "[deg]" << std::endl;

    std::cout << "pos_score :\t" << pos_score << " [m]" << std::endl;
    std::cout << "att_score :\t" << att_score*RAD2DEG << " [deg]"<< std::endl;

    return 0;
}