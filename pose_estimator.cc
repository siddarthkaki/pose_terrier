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
    double meas_std = 2.0*DEG2RAD;

    // specify rigid position vector of feature points wrt target in target frame
    MatrixXd rFeaMat(4,3);
    rFeaMat <<  0.0,    0.0,    0.5,
                0.0,    0.0,   -1.5,
                0.0,    1.0,    1.0,
                0.0,   -1.0,    1.0;
                //1.0,    0.0,    0.5;

    //------------------------------------------------------------------------/

    // initial state guess
    double posArr0[3] = { 0.0, 0.0, 25.0 };
    double quatArr0[4] = {1.0, 0.0, 0.0, 0.0}; // w,x,y,z

    // convert initial state information from double arrays to Eigen
    Pose state0;
    state0.pos(0) = posArr0[0];
    state0.pos(1) = posArr0[1];
    state0.pos(2) = posArr0[2];
    state0.quat.w() = quatArr0[0];
    state0.quat.x() = quatArr0[1];
    state0.quat.y() = quatArr0[2];
    state0.quat.z() = quatArr0[3]; 
    state0.quat.normalize();

    //-- Simulate Measurements -----------------------------------------------/

    // true state information
    double posArr [3] = {0.5377, 1.8339, 18.2235};
    //double quatArr [4] = {0.6937, -0.6773, 0.0642, 0.2365};
    double quatArr [4] = {0.6443,   -0.6230,    0.4369,    0.0767};
    
    // convert true state information from double arrays to Eigen
    Pose stateTrue;
    stateTrue.pos(0) = posArr[0];
    stateTrue.pos(1) = posArr[1];
    stateTrue.pos(2) = posArr[2];
    stateTrue.quat.w() = quatArr[0];
    stateTrue.quat.x() = quatArr[1];
    stateTrue.quat.y() = quatArr[2];
    stateTrue.quat.z() = quatArr[3]; 
    stateTrue.quat.normalize();

    // express feature points in chaser frame at the specified pose
    MatrixXd rMat = Utilities::FeaPointsTargetToChaser(stateTrue, rCamVec, rFeaMat);

    // generate simulated measurements
    VectorXd yVec = Utilities::SimulateMeasurements(rMat, focal_length);

    // add Gaussian noise to simulated measurements
    VectorXd yVecNoise = Utilities::AddGaussianNoiseToVector(yVec, meas_std);

    //------------------------------------------------------------------------/

    //-- Solve for pose ------------------------------------------------------/

    // timing
    std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

    // solve for pose with ceres (via wrapper)
    PoseSolution poseSol = PoseSolver::SolvePose(state0, yVecNoise, rCamVec, rFeaMat);

    Pose conj_state = Utilities::ConjugatePose(poseSol.state);

    // timing
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    
    // time taken to perform NLS solution
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>( t2 - t1 ).count();
    //------------------------------------------------------------------------/

    //-- Performance Metrics & Outputs ---------------------------------------/

    // compute position and attitude scores
    double      pos_score = Utilities::PositionScore(stateTrue.pos , poseSol.state.pos );
    double      att_score = Utilities::AttitudeScore(stateTrue.quat, poseSol.state.quat);
    double conj_att_score = Utilities::AttitudeScore(stateTrue.quat,    conj_state.quat);

    // print to command line
    std::cout << poseSol.summary.BriefReport() << "\n";
    //std::cout << poseSol.summary.FullReport() << "\n";
    
    
    std::cout << "posVec :\t" << state0.pos.transpose();
    std::cout << "\t->\t";
    std::cout << poseSol.state.pos.transpose() << " [m]" << std::endl;

    std::cout << "attVec :\t" << state0.quat.coeffs().transpose();
    std::cout << "\t\t->\t";
    std::cout << poseSol.state.quat.coeffs().transpose() << std::endl;
    std::cout << "\t\t\t\t\t" << conj_state.quat.coeffs().transpose() << std::endl;
    

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