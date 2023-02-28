/* Copyright (c) 2023 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef JESTOR_H_
#define JESTOR_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/SVD>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <iostream>
#include <math.h>

#include "Utilities.h"
#include "OsqpEigen/OsqpEigen.h"

#include "third_party/CppRot/cpprot.h"
#include "third_party/osqp_cpp/osqp++.h"

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::AngleAxisd;

//namespace JESTOR {

class JESTOR
{
    private:
        Matrix3d I33 = Matrix3d::Identity();
        MatrixXd I44 = MatrixXd::Identity(4, 4);

        Vector8d lb_, ub_;

        osqp::OsqpInstance instance;

        JESTORBatch A_; // matrix to store batch measurements

        const double kInfinity = std::numeric_limits<double>::infinity();

        Vector8d ConstrainedQP(const Matrix8d &B_reduced_, const Vector8d &b1_vec_, const Vector8d &lb_, const Vector8d &ub_);
        Vector8d ConstrainedQP2(const Matrix8d &B_reduced_, const Vector8d &b1_vec_, const Vector8d &lb_, const Vector8d &ub_);        
        MatrixXd f_Omega_BF(const Vector3d &omega);

    public:
    
        double dt_; // measurement time-step
        double tk_; // current time
        bool adapt_window_size_; // toggle window size adaptation
        unsigned int L_; // sliding window size 
        unsigned int Lmin_; // minimum feasible window size
        unsigned int Lmax_; // maximum feasible window size

        Matrix3d J_est_; // reduced MOI estimate
 
        JESTOR();
        void Init(
            const bool &adapt_window_size,
            const unsigned int &L, 
            const unsigned int &Lmin, 
            const unsigned int &Lmax,  
            const double &angle_noise_std_,
            const double &threshold_n
        );
        void Update(const Quaterniond &m_quat, const Vector3d &m_omega);
};

//} // end namespace

#endif // JESTOR_H_