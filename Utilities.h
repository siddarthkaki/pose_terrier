/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <chrono>
#include <random>
#include <algorithm>
#include <math.h>

#include "ceres/ceres.h"

#include "third_party/csvfile.h"
#include "third_party/CppRot/cpprot.h"
#include "third_party/eigenmvn/eigenmvn.h"

using Eigen::Vector2d;
using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
typedef Eigen::Matrix<double, 4, 1> Vector4d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 8, 1> Vector8d;
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> Matrix3d_rm;
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 6, 6, Eigen::RowMajor> Matrix6d_rm;
typedef Eigen::Matrix<double, 8, 8> Matrix8d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Matrix<double, 4, Eigen::Dynamic> MatrixQuat;
typedef Eigen::Matrix<double, Eigen::Dynamic, 9> JESTORBatch;

struct Pose
{
    Vector3d pos;
    Quaterniond quat;
};

struct PoseSolution
{
    Pose pose;
    Matrix6d cov_pose;
    ceres::Solver::Summary summary;
};

struct Twist
{
    Vector3d vel;
    Vector3d ang_vel;
};

class Utilities
{
    public:
        static Matrix3d Euler2DCM_312(const Vector3d& eulVec);
        static Vector3d DCM2Euler_312(const MatrixXd& DCM);
        static double UnwrapAngles(const double &old_angle, const double &new_angle);
        static MatrixXd FeaPointsTargetToChaser(const Pose& state, const Vector3d& rCamVec, const MatrixXd& rFeaMat);
        static Vector2d CameraProjection(const Vector3d& point3DVec, const double& f);
        static VectorXd SimulateMeasurements(const MatrixXd& rMat, const double& focal_length);
        static VectorXd AddGaussianNoiseToVector(const VectorXd& vec, const double& std);
        static Pose ConjugatePose(const Pose& state);
        static MatrixXd ConvertToEigenMatrix(double **data, unsigned int rows, unsigned int cols);
        static Vector4d QuatToVec4(const Quaterniond& quat);
        static Quaterniond Vec4ToQuat(const Vector4d& vec4);
        static double PositionScore(const Vector3d& pos, const Vector3d& posHat);
        static double AttitudeScore(const Quaterniond& quat, const Quaterniond& quatHat);
        static double StdVectorMean(const std::vector<double>& vec);
        static double StdVectorVar(const std::vector<double>& vec);
        static std::string WrapVarToPath(std::string varname);
        static void WritePosesToCSV(const std::vector<Pose>& vec, const std::string& filename, const bool& append_mode);
        static void WriteQuatsToCSV(const std::vector<Pose>& vec, const std::string& filename, const bool& append_mode);
        static void WriteKFStatesToCSV(const std::vector<VectorXd>& states, const std::string& filename, const bool& append_mode);
        static void WriteKFCovarsToCSV(const std::vector<MatrixXd>& covars, const std::string& filename, const bool& append_mode);
        static void WriteTimestampsToFile(const std::vector<double>& timestamps, const std::string& filename, const bool& append_mode);

        static constexpr double DEG2RAD = M_PI/180.0;
        static constexpr double RAD2DEG = 180.0/M_PI;
};

#endif // UTILITIES_H_