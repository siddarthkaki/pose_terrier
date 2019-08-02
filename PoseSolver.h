/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

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
        static PoseSolution SolvePose(const Pose& pose0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat);
        static PoseSolution SolvePoseReinit(const Pose& pose0, const VectorXd& yVec, const Vector3d& rCamVec, const MatrixXd& rFeaMat);
        static Twist TwoPointDiffTwistEstimator(const Pose& posei, const Pose& posej, const double& T);
};

#endif // POSESOLVER_H_