#ifndef POSESOLVER_H_
#define POSESOLVER_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math.h>

#include "ceres/ceres.h"
#include "cost_functor.h"

#include "third_party/eigenmvn/eigenmvn.h"

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class PoseSolver
{
    public:
        static VectorXd SolvePose(VectorXd yVec, VectorXd stateVec0, Vector3d rCamVec, MatrixXd rFeaMat);
};

#endif // POSESOLVER_H_