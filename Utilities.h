#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math.h>

#include "third_party/eigenmvn/eigenmvn.h"

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;

class Utilities
{
    public:
        static Matrix3d Euler2DCM_312(const Vector3d& eulVec);
        static MatrixXd FeaPointsTargetToChaser(const VectorXd& stateVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat);
        static Vector2d CameraProjection(const Vector3d& point3DVec, const double& f);
        static VectorXd SimulateMeasurements(const MatrixXd& rMat, const double& focal_length);
        static VectorXd AddGaussianNoiseToVector(const VectorXd& vec, const double& std);
        static double   PositionScore(const VectorXd& stateVec, const VectorXd& stateHatVec);
        static double   AttitudeScore(const VectorXd& stateVec, const VectorXd& stateHatVec);
};

#endif // UTILITIES_H_