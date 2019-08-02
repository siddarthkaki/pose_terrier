/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef KALMANFILTER_H_
#define KALMANFILTER_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;

namespace KF {

class KalmanFilter
{
    public:
        unsigned int num_states_;
        unsigned int num_measurements_;
        unsigned int num_inputs_;
        double dt_;

        MatrixXd Q_; // process_noise_covariance
        MatrixXd R_; // measurement_noise_covariance
        MatrixXd F_; // dynamics_model
        MatrixXd G_; // input_model
        MatrixXd H_; // measurement_model

        VectorXd statekk_;
        VectorXd statek1k_;
        VectorXd statek1k1_;

        MatrixXd covarkk_;
        MatrixXd covark1k_;
        MatrixXd covark1k1_;

        std::vector<VectorXd> states;
        std::vector<MatrixXd> covars;

        KalmanFilter();
        KalmanFilter(const unsigned int &num_states, const unsigned int &num_measurements, const unsigned int &num_inputs, const double &dt);
        void InitLinearPoseTracking(const double &process_noise_std, const double &measurement_noise_std, const double &dt);
        void SetInitialStateAndCovar(const VectorXd &state0, MatrixXd &covar0);
        void Predict(const VectorXd &input);
        void Update(const VectorXd &measurement);
        void KFStep(const VectorXd &measurement);
        void KFStep(const VectorXd &measurement, const VectorXd &input);
        void StoreAndClean();

        void PrintModelMatrices();
        //static double StdVectorVar(const std::vector<double>& vec);
};

} // end namespace

#endif // KALMANFILTER_H_