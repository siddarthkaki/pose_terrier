/* Copyright (c) 2022 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
#include "QuateRA.h"

namespace QuateRA {

    QuateRA::QuateRA(const double &dt)
    {
        dt_ = dt;
        processed_measurement_ = false;
    }

    void QuateRA::Init(
            const unsigned int &L, 
            const unsigned int &Lmin, 
            const unsigned int &Lmax,  
            const double &eps_mean, 
            const double &eps_std,  
            const double &threshold_n
    )
    {
        L_ = L;
        Lmin_ = Lmin;
        Lmax_ = Lmax;
        eps_mean_ = eps_mean;
        eps_std_ = eps_std;
        threshold_n_ = threshold_n;

        Qmat_ = Vector4d::Zero(); // TODO double check matrix ops on matrices with uninitialised values

        tvec_ = VectorXd::Zero(1,1);

        H_ = MatrixXd::Zero(3,6);
        H_.block(0, 0, 3, 3) = I33;

        eps_cr_hi_ = eps_mean_ + threshold_n_*eps_std_;
        eps_cr_lw_ = eps_mean_ - threshold_n_*eps_std_;

        ang_vel_est_ = Vector3d::Zero();
        covar_est_ = 10000.0*Matrix6d::Identity();

        processed_measurement_ = false;
    }

    // Init Measurement
    void QuateRA::InitMeasurement(const Vector4d &measurement, const Matrix3d &covar)
    {
        Qmat_ = measurement;
        tvec_(0) = 0.0;

        R_ = covar;
        FIM_ = Matrix6d::Zero();
        FIM_.topLeftCorner(3,3) = R_.inverse();
    }

    // Update step
    void QuateRA::Update(const Vector4d &measurement, const Matrix3d &covar, const double &tk)
    {
        unsigned int ql = Qmat_.cols(); 
        Qmat_.conservativeResize(4,ql+1);
        Qmat_.rightCols(1) = measurement;

        tvec_.conservativeResize(ql+1);
        tvec_(ql) = tk;

        ql++;

        if (ql > 3)
        {

            if (ql > L_)
            {
                //std::cout << Qmat_.rightCols(1).transpose() << std::endl << std::endl;
                //std::cout << tvec_(ql-1) << " - " << tvec_(0) << " = " << tvec_(ql-1) - tvec_(0) << " " << std::endl << std::endl;

                MatrixQuat Qmat_new = Qmat_.rightCols(L_);
                VectorXd tvec_new = tvec_.tail(L_);
                Qmat_ = Qmat_new;
                tvec_ = tvec_new;
                ql = L_;

                //std::cout << Qmat_.rightCols(1).transpose() << std::endl << std::endl;
                //std::cout << tvec_(ql-1) << " - " << tvec_(0) << " = " << tvec_(ql-1) - tvec_(0) << " ";

            }

            //std::cout << Qmat_.rightCols(1).transpose() << std::endl << std::endl;
            //std::cout << tvec_.transpose() << std::endl << std::endl;


            Matrix4d Zmat_ = Qmat_*Qmat_.transpose();
            Eigen::JacobiSVD<MatrixXd> svd( Zmat_, Eigen::ComputeFullU | Eigen::ComputeFullV);

            Vector4d u1 = svd.matrixU().col(0);
            Vector4d u2 = svd.matrixU().col(1);

            Vector4d sigmas = svd.singularValues();

            // spin-axis direction
            Quaterniond u1_quat = Utilities::Vec4ToQuat(u1);
            Quaterniond u2_quat = Utilities::Vec4ToQuat(u2);
            Quaterniond SADk1_quat = CppRot::QuatMult_S( u2_quat, u1_quat.inverse() );

            Vector3d SADk1 = SADk1_quat.vec();

            // measurement re-projection onto spin-axis plane
            // MatrixQuat Qproj = MatrixQuat::Zero(4,ql);
            VectorXd Phivec = VectorXd::Zero(ql);
            MatrixXd Hmat = MatrixXd::Zero(ql,2);

            Vector4d q0_proj = ProjectQuatToPlane(Qmat_.col(0), u1, u2);

            Matrix4d SADk1_cross = Matrix4d::Zero();
            SADk1_cross.block(0,1,1,3) = -SADk1.transpose();
            SADk1_cross.block(1,0,3,1) = SADk1;
            SADk1_cross.bottomRightCorner(3,3) = -CppRot::CrossProductEquivalent(SADk1);
            Vector4d q0_proj_perp = SADk1_cross*q0_proj;

            for (unsigned int qdx = 0; qdx < ql; qdx++)
            {
                Vector4d qtemp = Qmat_.col(qdx);
                
                Vector4d qtemp_proj = ProjectQuatToPlane(qtemp, u1, u2);

                // Qproj.col(qdx) = qtemp_proj;
                
                double a1 = (qtemp_proj.transpose()*q0_proj).value();
                double a2 = (qtemp_proj.transpose()*q0_proj_perp).value();

                double Phi_temp = 2*atan2(a2, a1);
                //double Phi_temp = 2*atan2(qu2,qu1);
                Phi_temp = atan2(sin(Phi_temp), cos(Phi_temp));
                if (qdx == 0)
                {
                    Phivec(qdx) = Phi_temp;
                }
                else
                {
                    Phivec(qdx) = Utilities::UnwrapAngles(Phivec(qdx-1), Phi_temp);
                }
                Hmat.row(qdx) << 1, tvec_(qdx) - tvec_(0);

                //if (qdx == ql - 1) { std::cout << tvec_(qdx) << " - " << tvec_(0) << " = " << tvec_(qdx) - tvec_(0) << " "; }
            }
            //std::cout << std::endl << std::endl;

            // least-squares solution to spin rate
            VectorXd Xhat = ((Hmat.transpose()*Hmat).inverse())*Hmat.transpose()*Phivec;
            double OmegaEstNormk1 = Xhat(1);

            ang_vel_est_ = OmegaEstNormk1*SADk1;

            // covariance propagation
            R_ = covar;

            A_ = Matrix6d::Zero();
            A_.topLeftCorner(3,3) = -CppRot::CrossProductEquivalent(ang_vel_est_);
            A_.topRightCorner(3,3) = I33;

            double dt = tvec_(ql-1) - tvec_(ql-2);
            F_ = ( -A_*dt ).exp();
            Matrix6d Finv_ = F_.inverse();

            FIM_ = Finv_.transpose()*FIM_*Finv_ + H_.transpose()*R_.inverse()*H_;
            covar_est_ = FIM_.inverse();

            //std::cout << covar_est_ << std::endl << std::endl;

            // TODO adaptive sliding window length
        }

        
    }

    Vector4d QuateRA::ProjectQuatToPlane(const Vector4d &quat, const Vector4d &u1, const Vector4d &u2)
    {
        double qu1 = (quat.transpose()*u1).value();
        double qu2 = (quat.transpose()*u2).value();

        return (qu1*u1 + qu2*u2) / sqrt(qu1*qu1 + qu2*qu2);
    }


    void QuateRA::PrintModelMatrices()
    {
        std::cout << "R:\t" << std::endl << R_ << std::endl << std::endl;
        std::cout << "F:\t" << std::endl << F_ << std::endl << std::endl;
        std::cout << "A:\t" << std::endl << A_ << std::endl << std::endl;
        std::cout << "H:\t" << std::endl << H_ << std::endl << std::endl;
    }

} // end namespace
