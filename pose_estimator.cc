#include <Eigen/Core>
#include <math.h>

#include "ceres/ceres.h"
#include "glog/logging.h"

using ceres::AutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;

using Eigen::Vector3d;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;

// A templated cost functor for automatic differentiation
// Measurements: yVec
// Feature Points: rFeaMat
// Camera Location: rCamVec
// States: posVec, eulVec
class MeasResidCostFunctor
{
    public: MeasResidCostFunctor(VectorXd yVec, MatrixXd rFeaMat, Vector3d rCamVec): yVec_(yVec), rFeaMat_(rFeaMat), rCamVec_(rCamVec_) {}

    template <typename T>
    bool operator()(const T* const posArr,
                    const T* const eulArr,
                    T* residuals) const
    {

        Eigen::Matrix<T, 3, 1> posVec;
        posVec(0) = posArr[0];
        posVec(1) = posArr[1];
        posVec(2) = posArr[2];

        int numPts = 4;//rFeaMat_.rows();
        int numMeas = yVec_.rows();

        // x[0] = static_cast<T>(0) <- way to keep in Jet type

        for (int idx = 0; idx < numPts; idx++)
        {
            Vector3d rFeaVeci = rFeaMat_.block(idx,0,1,3).transpose();
            
            //residuals[2*idx+0] = measResidVec(2*idx+0);
            //residuals[2*idx+1] = measResidVec(2*idx+1);

            residuals[2*idx+0] = yVec_(2*idx+0) - static_cast<T>(atan2(posVec(0) - rCamVec_(0) + rFeaVeci(0)*(cos(eulArr[2])*cos(eulArr[1]) - sin(eulArr[2])*sin(eulArr[1])*sin(eulArr[0])) + rFeaVeci(1)*(cos(eulArr[1])*sin(eulArr[2]) + cos(eulArr[2])*sin(eulArr[1])*sin(eulArr[0])) - rFeaVeci(2)*sin(eulArr[1])*cos(eulArr[0]), posVec(2) - rCamVec_(2) + rFeaVeci(0)*(cos(eulArr[2])*sin(eulArr[1]) + cos(eulArr[1])*sin(eulArr[2])*sin(eulArr[0])) + rFeaVeci(1)*(sin(eulArr[2])*sin(eulArr[1]) - cos(eulArr[2])*cos(eulArr[1])*sin(eulArr[0])) + rFeaVeci(2)*cos(eulArr[1])*cos(eulArr[0])));
            residuals[2*idx+1] = yVec_(2*idx+1) - static_cast<T>(atan2(posVec(1) - rCamVec_(1) + rFeaVeci(2)*sin(eulArr[0]) + rFeaVeci(1)*cos(eulArr[2])*cos(eulArr[0]) - rFeaVeci(0)*sin(eulArr[2])*cos(eulArr[0]), posVec(2) - rCamVec_(2) + rFeaVeci(0)*(cos(eulArr[2])*sin(eulArr[1]) + cos(eulArr[1])*sin(eulArr[2])*sin(eulArr[0])) + rFeaVeci(1)*(sin(eulArr[2])*sin(eulArr[1]) - cos(eulArr[2])*cos(eulArr[1])*sin(eulArr[0])) + rFeaVeci(2)*cos(eulArr[1])*cos(eulArr[0])));
            
//            yHatVec(2*idx+0) = atan2(posVec(0) - rCamVec_(0) + rFeaVeci(0)*(cos(psi)*cos(theta) - sin(psi)*sin(theta)*sin(phi)) + rFeaVeci(1)*(cos(theta)*sin(psi) + cos(psi)*sin(theta)*sin(phi)) - rFeaVeci(2)*sin(theta)*cos(phi), posVec(2) - rCamVec_(2) + rFeaVeci(0)*(cos(psi)*sin(theta) + cos(theta)*sin(psi)*sin(phi)) + rFeaVeci(1)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(2)*cos(theta)*cos(phi));
//            yHatVec(2*idx+1) = atan2(posVec(1) - rCamVec_(1) + rFeaVeci(2)*sin(phi) + rFeaVeci(1)*cos(psi)*cos(phi) - rFeaVeci(0)*sin(psi)*cos(phi), posVec(2) - rCamVec_(2) + rFeaVeci(0)*(cos(psi)*sin(theta) + cos(theta)*sin(psi)*sin(phi)) + rFeaVeci(1)*(sin(psi)*sin(theta) - cos(psi)*cos(theta)*sin(phi)) + rFeaVeci(2)*cos(theta)*cos(phi));
        }

        return true;
    }

    private:
        VectorXd yVec_;
        MatrixXd rFeaMat_;
        Vector3d rCamVec_;

    public:
    // Factory to hide the construction of the CostFunction object from
    // the client code.
    static ceres::CostFunction* Create( const VectorXd yVec,
                                        const MatrixXd rFeaMat,
                                        const Vector3d rCamVec_ )
    {
        return (new ceres::AutoDiffCostFunction<MeasResidCostFunctor, 8, 3, 3>(
            new MeasResidCostFunctor(yVec, rFeaMat, rCamVec_)));

    }

};



int main(int argc, char** argv)
{
    google::InitGoogleLogging(argv[0]);

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

    int numMeas = rFeaMat.rows();

    // Build the problem.
    Problem problem;

    // Set up the only cost function (also known as residual). This uses
    // auto-differentiation to obtain the derivative (jacobian).
    //ceres::CostFunction* cost_function =
    //    MeasResidCostFunctor::Create(yVec, rFeaMat, rCamVec);
    CostFunction* cost_function = new AutoDiffCostFunction<MeasResidCostFunctor, 8, 3, 3>(
            new MeasResidCostFunctor(yVec, rFeaMat, rCamVec));
    
    problem.AddResidualBlock(cost_function, NULL, posArr, eulArr);

    // Run the solver!
    Solver::Options options;
    options.minimizer_type = ceres::TRUST_REGION;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
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