#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math.h>

#include "ceres/ceres.h"
#include "cost_functor.h"

#include "third_party/eigenmvn/eigenmvn.h"

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;

class Utilities
{
    public:
        static Matrix3d Euler2DCM_312(Vector3d eulVec);
        static MatrixXd FeaPointsTargetToChaser(VectorXd stateVec, Vector3d rCamVec, MatrixXd rFeaMat);
        static Vector2d CameraProjection(Vector3d point3DVec, double f);
        static VectorXd SimulateMeasurements(MatrixXd rMat, double focal_length);
        static VectorXd AddNoiseToMeasurements(VectorXd yVec, double std);
        static double   PositionScore(VectorXd stateVec, VectorXd stateHatVec);
        static double   AttitudeScore(VectorXd stateVec, VectorXd stateHatVec);
        static VectorXd SolvePoseReinit(VectorXd yVec, VectorXd xHatVec0, Vector3d rCamVec, MatrixXd rFeaMat);
};

#endif // UTILITIES_H_