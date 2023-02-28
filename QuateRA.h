/* Copyright (c) 2022 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef QUATERA_H_
#define QUATERA_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
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

//namespace QuateRA {

class QuateRA
{
    private:
        bool processed_measurement_;
        Vector3d delta_gibbs_est_;
        VectorXd state_est_;

        Matrix3d I33 = Matrix3d::Identity();
        MatrixXd I44 = MatrixXd::Identity(4, 4);

        MatrixXd Hk_;

        //Vector3d e1uv;

        Vector4d ProjectQuatToPlane(const Vector4d &quat, const Vector4d &u1, const Vector4d &u2);

    public:
    
        double dt_; // measurement time-step
        double tk_; // current time
        bool adapt_window_size_; // toggle window size adaptation
        unsigned int L_; // sliding window size 
        unsigned int Lmin_; // minimum feasible window size
        unsigned int Lmax_; // maximum feasible window size

        // double euler_noise_std_; // orientation measurement std. dev.

        double angle_noise_std_; // std. dev. of attitude estimation error
        double eps_mean_; // statistically determined mean of window size thresholding value
        double eps_std_; // statistically determined std. dev. of window size thresholding value
        double threshold_n_;

        std::vector<double> L_history;
        std::vector<double> eps_history;
        
        Matrix3d R_; // measurement noise covariance
        Matrix6d F_; // covariance propagation dynamics model
        Matrix6d A_; // orientation propagation dynamics model
        MatrixXd H_; // measurement model

        MatrixQuat Qmat_; // quaternion measurement storage
        //MatrixXd RUVmat_; // rotated unit vector storage 

        VectorXd tvec_; // measurement time storage

        Matrix6d FIM_; // Fisher Information Matrix

        Vector3d ang_vel_est_;
        Matrix6d covar_est_;

        QuateRA(const double &dt);
        void Init(
            const bool &adapt_window_size,
            const unsigned int &L, 
            const unsigned int &Lmin, 
            const unsigned int &Lmax,  
            const double &angle_noise_std_,
            const double &threshold_n
        );
        void InitMeasurement(const Vector4d &measurement, const Matrix3d &covar);
        void Update(const Vector4d &measurement, const Matrix3d &covar, const double &tk);
        
        void PrintModelMatrices();
};

//} // end namespace

#endif // QUATERA_H_