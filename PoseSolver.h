#ifndef POSESOLVER_H_
#define POSESOLVER_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math.h>

#include "ceres/ceres.h"
#include "cost_functor.h"
#include "Utilities.h"

#include "third_party/eigenmvn/eigenmvn.h"

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Quaterniond;

class PoseSolver
{
    public:
        static PoseSolution SolvePose(const Pose& state0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat);
        static PoseSolution SolvePoseReinit(const Pose& state0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat);
};

#endif // POSESOLVER_H_