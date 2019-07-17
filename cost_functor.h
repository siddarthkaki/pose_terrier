#ifndef COST_FUNCTOR_H_
#define COST_FUNCTOR_H_

#include <Eigen/Core>
#include <math.h>

#include "ceres/ceres.h"

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

        int numPts = rFeaMat_.rows();
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
};


// A templated cost functor for automatic differentiation
// Measurements: yVec
// Feature Points: rFeaMat
// Camera Location: rCamVec
// States: posVec, quatVec
class MeasResidCostFunctorQuat
{
    public: MeasResidCostFunctorQuat(VectorXd yVec, MatrixXd rFeaMat, Vector3d rCamVec): yVec_(yVec), rFeaMat_(rFeaMat), rCamVec_(rCamVec_) {}

    template <typename T>
    bool operator()(const T* const posArr,
                    const T* const quatArr,
                    T* residuals) const
    {

        /*Eigen::Matrix<T, 3, 1> posVec;
        posVec(0) = posArr[0];
        posVec(1) = posArr[1];
        posVec(2) = posArr[2];*/

        // Map T* array to Eigen Vector3 with correct Scalar type
        Eigen::Matrix<T,3,1> posVec = Eigen::Map<const Eigen::Matrix<T,3,1>>(posArr);

        // Map T* array to Eigen Quaternion with correct Scalar type
        Eigen::Quaternion<T> quat = Eigen::Map<const Eigen::Quaternion<T>>(quatArr);
        quat.normalize();

        int numPts = rFeaMat_.rows();
        int numMeas = yVec_.rows();

        // x[0] = static_cast<T>(0) <- way to keep in Jet type

        for (int idx = 0; idx < numPts; idx++)
        {
            Vector3d rFeaVeci = rFeaMat_.block(idx,0,1,3).transpose();
            Eigen::Matrix<T,3,1> rFeaVeciJet;
            rFeaVeciJet << T(rFeaVeci(0)), T(rFeaVeci(1)), T(rFeaVeci(2));

            Eigen::Matrix<T,3,1> rFeaVeciRotated = quat*rFeaVeciJet;
            
            // position vector of feature point 1 wrt chaser in chaser frame
            Eigen::Matrix<T,3,1> rFeaVeciRotatedTranslated = posVec - rCamVec_ + rFeaVeciRotated;
            
            residuals[2*idx+0] = yVec_(2*idx+0) - static_cast<T>(atan2( rFeaVeciRotatedTranslated(0), rFeaVeciRotatedTranslated(2) ));
            residuals[2*idx+1] = yVec_(2*idx+1) - static_cast<T>(atan2( rFeaVeciRotatedTranslated(1), rFeaVeciRotatedTranslated(2) ));
            
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
    /*
    static ceres::CostFunction* Create( const VectorXd yVec,
                                        const MatrixXd rFeaMat,
                                        const Vector3d rCamVec )
    {
        return (new ceres::AutoDiffCostFunction<MeasResidCostFunctor, 8, 3, 3>(
            new MeasResidCostFunctor(yVec, rFeaMat, rCamVec)));

    }
    */

};

#endif // COST_FUNCTOR_H_