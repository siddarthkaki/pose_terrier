#include <Eigen/Core>
#include <chrono>
#include <iomanip>
#include <math.h>

#include "ceres/ceres.h"
//#include "glog/logging.h"

#include "cost_functor.h"
#include "Utilities.h"
#include "PoseSolver.h"

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

    //-- Set-up problem geometry and params ----------------------------------/

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

    //------------------------------------------------------------------------/

    // initial state guess
    double posArr0[3] = { -1.0, 3.0, 25.0 };
    double eulArr0[3] = {0.0, 0.0, 0.0};//{ 1.5, 0.4, -2.2 };

    // convert initial state information from double arrays to Eigen
    VectorXd stateVec0(6);
    stateVec0.head(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(posArr0);
    stateVec0.tail(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(eulArr0);

    //-- Simulate Measurements -----------------------------------------------/

    // true state information
    double posArr [3] = { 0.5377, 1.8339, 18.2235 };
    double eulArr [3] = {2.0, -1.75, -1.5};//{ 1.3543, 0.5007, -2.0541 };

    // convert true state information from double arrays to Eigen
    VectorXd stateVec(6);
    stateVec.head(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(posArr);
    stateVec.tail(3) = Eigen::Map<Eigen::Matrix<double,3,1>>(eulArr);

    // express feature points in chaser frame at the specified pose
    MatrixXd rMat = Utilities::FeaPointsTargetToChaser(stateVec, rCamVec, rFeaMat);

    // generate simulated measurements
    VectorXd yVec = Utilities::SimulateMeasurements(rMat, focal_length);

    // add Gaussian noise to simulated measurements
    VectorXd yVecNoise = Utilities::AddGaussianNoiseToVector(yVec, meas_std);

    //------------------------------------------------------------------------/

    //-- Solve for pose ------------------------------------------------------/

    // timing
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // solve for pose with ceres (via wrapper)
    PoseSolution poseSol = PoseSolver::SolvePoseReinit(yVecNoise, stateVec0, rCamVec, rFeaMat);

    // timing
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    
    // time taken to perform NLS solution
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    //------------------------------------------------------------------------/

    //-- Performance Metrics & Outputs ---------------------------------------/

    // compute position and attitude scores
    double pos_score = Utilities::PositionScore(stateVec, poseSol.stateHatVec);
    double att_score = Utilities::AttitudeScore(stateVec, poseSol.stateHatVec);

    // print to command line
    std::cout << poseSol.summary.BriefReport() << "\n";
    //std::cout << poseSol.summary.FullReport() << "\n";
    
    /*
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
    */

    std::cout << "pos_score :\t" << pos_score << " [m]" << std::endl;
    std::cout << "att_score :\t" << att_score*RAD2DEG << " [deg]"<< std::endl;

    std::cout << "Time taken by program is : "  << std::setprecision(9)
                                                << (double)duration
                                                << " [ms]"
                                                << std::endl; 

    //------------------------------------------------------------------------/

    return 0;
}