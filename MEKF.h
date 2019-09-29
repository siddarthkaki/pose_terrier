/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef MEKF_H_
#define MEKF_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <iostream>
#include <math.h>

#include "Utilities.h"

#include "third_party/CppRot/cpprot.h"

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::AngleAxisd;

namespace MEKF {

class MEKF
{
    private:
        bool processed_measurement_;
        Vector3d delta_gibbs_est_;

    public:
        unsigned int num_states_;
        unsigned int num_measurements_;
        double dt_;
        double tau_;

        MatrixXd Q_; // process_noise_covariance
        MatrixXd R_; // measurement_noise_covariance
        MatrixXd F_; // covariance_propagation_dynamics_model
        MatrixXd A_; // quaternion_propagation_dynamics_model
        MatrixXd G_; // input_model
        MatrixXd H_; // measurement_model

        Quaterniond quat_est_;
        Vector3d omega_est_;
        Vector3d alpha_est_;

        MatrixXd covar_est_;

        //MatrixXd covarkk_;
        //MatrixXd covark1k_;
        //MatrixXd covark1k1_;

        //std::vector<VectorXd> states;
        //std::vector<MatrixXd> covars;

        MEKF(const double &dt);
        MEKF(const unsigned int &num_states, const unsigned int &num_measurements, const double &dt);
        void Init(const double &process_noise_std, const double &measurement_noise_std, const double &dt);
        void SetInitialStateAndCovar(const Quaterniond &quat0, const Vector3d &omega0, const Vector3d &alpha0, const MatrixXd &covar0);
        void Predict();
        void Update(const VectorXd &measurement);
        void Reset();
        void StoreAndClean();

        //static VectorXd MeasResidFunction(const VectorXd &measurement, const VectorXd &statek1k, const double &dt);

        void PrintModelMatrices();
        //static double StdVectorVar(const std::vector<double>& vec);
};

} // end namespace

#endif // MEKF_H_