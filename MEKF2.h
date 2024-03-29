/* Copyright (c) 2021 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef MEKF2_H_
#define MEKF2_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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

//namespace MEKF2 {

class MEKF2
{
    private:
        bool processed_measurement_;
        Vector3d delta_gibbs_est_;

        Matrix3d I33 = Matrix3d::Identity();
        MatrixXd I44 = MatrixXd::Identity(4, 4);

    public:
        
        double dt_;
        double tau_;
        double delta_min_;
        unsigned int num_att_states_;
        unsigned int num_pos_states_;
        unsigned int num_states_;

        unsigned int num_att_measurements_;
        unsigned int num_pos_measurements_;
        unsigned int num_measurements_;
        double max_flip_thresh_deg_;
        MatrixXd Q_;     // process_noise_covariance
        MatrixXd R_;     // measurement_noise_covariance
        MatrixXd F_pos_; // position_covariance_propagation_dynamics_model
        MatrixXd F_;     // covariance_propagation_dynamics_model
        MatrixXd A_;     // quaternion_propagation_dynamics_model
        MatrixXd H_;     // measurement_model

        Matrix3d J_;     // moment-of-inertia matrix

        //Vector3d alpha_est_;
        //VectorXd x_est_;

        VectorXd state_est_;

        Vector3d pos_est_;
        Quaterniond quat_est_;
        Vector3d omega_est_;
        MatrixXd covar_est_;
        Matrix3d omega_covar_est_;

        //MatrixXd covarkk_;
        //MatrixXd covark1k_;
        //MatrixXd covark1k1_;

        //std::vector<VectorXd> states;
        //std::vector<MatrixXd> covars;

        MEKF2(const double &dt);
        MEKF2(const unsigned int &num_states, const unsigned int &num_measurements, const double &dt);
        void Init(const double &process_noise_std, const double &measurement_noise_std, const double &dt);
        void Init(
            const double &process_noise_std, 
            const double &measurement_noise_std, 
            const double &dt, 
            const double &tau,  
            const double &qpsd, 
            const double &max_flip_thresh_deg
        );
        void SetInitialStateAndCovar(const Quaterniond &quat0, const Vector3d &omega0, const Vector3d &alpha0, const VectorXd &x0, const MatrixXd &covar0);
        void Predict();
        void PredictEuler();
        void Update(const VectorXd &measurement);
        void Reset();
        void StoreAndClean();
        void AngVelUpdate(const Vector3d &measurement, const Matrix3d &covar);

        Matrix3d EulerDynamicsJacobian(const Vector3d &omega, const Matrix3d &J);

        void PrintModelMatrices();
};

//} // end namespace

#endif // MEKF2_H_