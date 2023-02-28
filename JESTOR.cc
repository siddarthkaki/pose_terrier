/* Copyright (c) 2023 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
#include "JESTOR.h"

//namespace JESTOR {

    JESTOR::JESTOR()
    {
        L_ = 10;
        lb_ << -kInfinity, -kInfinity,        0.0, -kInfinity,        0.0, -kInfinity, -kInfinity, -kInfinity;
        ub_ <<  kInfinity,  kInfinity,  kInfinity,  kInfinity,  kInfinity,  kInfinity,  kInfinity,  kInfinity; 
        //ub <<  kInfinity,  kInfinity,        1.0,  kInfinity,        1.0,  kInfinity,  kInfinity,  kInfinity; 
            
    }

    void JESTOR::Init(
            const bool &adapt_window_size,
            const unsigned int &L, 
            const unsigned int &Lmin, 
            const unsigned int &Lmax,  
            const double &angle_noise_std,
            const double &threshold_n
    )
    {
        adapt_window_size_ = adapt_window_size;
        L_ = L;
        Lmin_ = Lmin;
        Lmax_ = Lmax;
    }

    // Update step
    void JESTOR::Update(const Quaterniond &m_quat, const Vector3d &m_omega)
    {
        std::cout << "Quat : " << m_quat << std::endl;
        std::cout << "Omega: " << m_omega.transpose() << std::endl << std::endl;

        unsigned int al = A_.rows(); 
        A_.conservativeResize(al+3, 9);
        A_.bottomRows(3).leftCols(6) = f_Omega_BF(m_omega);
        A_.bottomRows(3).rightCols(3) = -CppRot::Quat2Tmat(m_quat);

        al = A_.rows();

        if (al > 3*8)
        {
            if (al > L_)
            {
                //std::cout << Qmat_.rightCols(1).transpose() << std::endl << std::endl;
                //std::cout << tvec_(ql-1) << " - " << tvec_(0) << " = " << tvec_(ql-1) - tvec_(0) << " " << std::endl << std::endl;

                // MatrixQuat Qmat_new = Qmat_.rightCols(L_);
                // //MatrixXd RUVmat_new = RUVmat_.rightCols(L_);
                // VectorXd tvec_new = tvec_.tail(L_);
                // Qmat_ = Qmat_new;
                // //RUVmat_ = RUVmat_new;
                // tvec_ = tvec_new;
                // ql = L_;

                //std::cout << Qmat_.rightCols(1).transpose() << std::endl << std::endl;
                //std::cout << tvec_(ql-1) << " - " << tvec_(0) << " = " << tvec_(ql-1) - tvec_(0) << " ";
            }

            Matrix9d B_ = A_.transpose()*A_;
            Eigen::JacobiSVD<MatrixXd> svd;
            svd.compute(B_, Eigen::ComputeFullU | Eigen::ComputeFullV);

            Vector9d sigmas = svd.singularValues();
            sigmas(8) = 0;

            Matrix9d S_tilde = MatrixXd::Zero(9, 9);
            S_tilde.diagonal() << sigmas;

            Matrix9d B_tilde = svd.matrixU()*S_tilde*svd.matrixV().transpose();

            Matrix8d B_reduced = B_tilde.bottomRightCorner(8,8);
            B_reduced = (B_reduced + B_reduced.transpose())/2.0;

            Vector8d b1_vec = B_tilde.col(0).tail(8);

            Vector8d x_reduced;

            // unconstrained solution
            //x_reduced = -B_reduced.inverse()*b1_vec;
    
            // constrained solution
            x_reduced = ConstrainedQP(B_reduced, b1_vec, lb_, ub_);

            Vector9d x_;
            x_ << 1.0, x_reduced;

            // updated estimate of reduced MOI matrix
            J_est_.row(0) << x_(0), x_(1), x_(2);
            J_est_.row(1) << x_(1), x_(3), x_(4);
            J_est_.row(2) << x_(2), x_(4), x_(5);

            std::cout << "J_est_diag: " << J_est_.diagonal().transpose() << std::endl << std::endl;
        }
    }

    // Constrained quadratic programming solver wrapper for osqp-cpp
    Vector8d JESTOR::ConstrainedQP(const Matrix8d &B_reduced_, const Vector8d &b1_vec_, const Vector8d &lb_, const Vector8d &ub_)
    {        
        Eigen::SparseMatrix<double>  objective_matrix =               B_reduced_.sparseView();
        Eigen::SparseMatrix<double>  objective_vector =                  b1_vec_.sparseView();
        Eigen::SparseMatrix<double> constraint_matrix = MatrixXd::Identity(8, 8).sparseView();

        //objective_matrix = Matrix8d::Identity().sparseView();
        //objective_vector = Vector8d::Zero().sparseView();

        instance.objective_matrix = objective_matrix;
        instance.objective_vector = objective_vector;
        instance.constraint_matrix = constraint_matrix;
        instance.lower_bounds = lb_;
        instance.upper_bounds = ub_;

        osqp::OsqpSolver solver;
        osqp::OsqpSettings settings;
        // Edit settings if appropriate.
        auto status = solver.Init(instance, settings);
        // Assuming status.ok().
        osqp::OsqpExitCode exit_code = solver.Solve();
        // Assuming exit_code == OsqpExitCode::kOptimal.
        double optimal_objective = solver.objective_value();
        Vector8d optimal_solution = solver.primal_solution();

        return optimal_solution;
    }

    // TODO
    MatrixXd JESTOR::f_Omega_BF(const Vector3d &omega)
    {
        MatrixXd Omega_BF = MatrixXd::Zero(3, 6);
        Omega_BF.row(0) << omega(0), omega(1), omega(2),      0.0,      0.0,      0.0;
        Omega_BF.row(1) <<      0.0, omega(0),      0.0, omega(1), omega(2),      0.0;
        Omega_BF.row(2) <<      0.0,      0.0, omega(0),      0.0, omega(1), omega(2);

        return Omega_BF;
    }

//} // end namespace
