/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
#include "MEKF.h"

namespace MEKF {

    MEKF::MEKF(const double &dt)
    {
        dt_ = dt;
        processed_measurement_ = false;
    }

    void MEKF::Init(const double &process_noise_std, const double &measurement_noise_std, const double &dt)
    {
        num_states_ = 9; // delta_gibbs(3), omega(3), alpha(3)
        num_measurements_ = 3; // delta_gibbs (3)
        dt_ = dt;
        tau_ = 1.0; // TODO PASS AS ARGUMENT

        Q_ = MatrixXd::Identity(num_states_, num_states_)*pow(process_noise_std,2); // process_noise_covariance
        Q_(6,6) *= 1e3;
        Q_(7,7) *= 1e3;
        Q_(8,8) *= 1e3;
        
        R_ = MatrixXd::Identity(num_measurements_, num_measurements_)*pow(measurement_noise_std,2); // measurement_noise_covariance

        A_ = MatrixXd::Identity(4, 4); // quaternion_propagation

        F_ = MatrixXd::Identity(num_states_, num_states_); // covariance_propagation

        H_ = MatrixXd::Zero(num_measurements_, num_states_); // measurement_model
        H_.block(0,0,num_measurements_,num_measurements_) = Matrix3d::Identity();

         quat_est_ = Quaterniond::Identity();
        omega_est_ = Vector3d::Zero();
        alpha_est_ = Vector3d::Zero();
        delta_gibbs_est_ = Vector3d::Zero();

        covar_est_ = MatrixXd::Zero(num_states_, num_states_);

        processed_measurement_ = false;
    }

    void MEKF::SetInitialStateAndCovar(const Quaterniond &quat0, const Vector3d &omega0, const Vector3d &alpha0, const MatrixXd &covar0)
    {
         quat_est_ =  quat0;
        omega_est_ = omega0;
        alpha_est_ = alpha0;
        covar_est_ = covar0;
    }

    // Prediction step
    void MEKF::Predict()
    {
        double omega_norm = omega_est_.norm();
        Vector3d omega_hat = omega_est_ / omega_norm;

        Matrix3d I33 = Matrix3d::Identity();
        MatrixXd I44 = MatrixXd::Identity(4, 4);
        MatrixXd omega_hat_44_equivalent = MatrixXd::Zero(4, 4);
        omega_hat_44_equivalent.block(0, 1, 1, 3) = -omega_hat.transpose();
        omega_hat_44_equivalent.block(1, 0, 3, 1) =  omega_hat;
        omega_hat_44_equivalent.block(1, 1, 3, 3) = -CppRot::CrossProductEquivalent(omega_hat);

        double phi = 0.5 * omega_norm * dt_;

        A_ = cos(phi) * I44 + sin(phi) * omega_hat_44_equivalent;
        F_.block(0, 0, 3, 3) = (-CppRot::CrossProductEquivalent(omega_est_) * dt_).exp();
        F_.block(3, 3, 3, 3) = I33;
        F_.block(3, 6, 3, 3) = I33 * dt_;
        F_.block(6, 6, 3, 3) = I33 * exp(-dt_ / tau_);

        // propagate quaternion
        quat_est_ = Utilities::Vec4ToQuat( A_ * Utilities::QuatToVec4(quat_est_) );

        // propagate omega vector
        omega_est_ = F_.block(3, 3, 3, 3)*omega_est_ + F_.block(3, 6, 3, 3)*alpha_est_;

        // propagate alpha vector
        alpha_est_ = F_.block(6, 6, 3, 3)*alpha_est_;

        // propagate covariance
        covar_est_ = F_ * covar_est_ * F_.transpose() + Q_;
    }

    // Update step
    void MEKF::Update(const VectorXd &measurement)
    {
        // Kalman gain
        MatrixXd K = covar_est_ * H_.transpose() * ((H_ * covar_est_ * H_.transpose() + R_).inverse());

        // delta Gibbs update
        Quaterniond quat_meas = Utilities::Vec4ToQuat(measurement.head(4));
        Quaterniond delta_quat = CppRot::QuatMult_S(quat_meas, quat_est_.inverse());
        Vector3d meas_innovation = 2.0 * delta_quat.vec() / delta_quat.w();
        delta_gibbs_est_ = K.block(0, 0, 3, 3)*meas_innovation;
        
        // Joseph update (general)
        MatrixXd I = MatrixXd::Identity(num_states_, num_states_);
        covar_est_ = (I - K * H_) * covar_est_ * ((I - K * H_).transpose()) + K * R_ * (K.transpose());
    }

    // Reset step
    void MEKF::Reset()
    {
        Quaterniond delta_quat = Quaterniond::Identity();
        delta_quat.w() = 1.0;
        delta_quat.vec() = 0.5*delta_gibbs_est_;
        Quaterniond quat_star = CppRot::QuatMult_S(delta_quat, quat_est_).normalized();
        
        
        // NOTE: heuristic method to ignore 180 deg pose ambiguities
        Quaterniond dq = CppRot::QuatMult_S(quat_est_, quat_star.inverse());
        double dangle = 2.0*acos( abs( dq.w() ) );
        if (dangle < 20.0*Utilities::DEG2RAD)
        {
            quat_est_ = quat_star;
        }
        
       
        processed_measurement_ = true;
    }

    void MEKF::StoreAndClean()
    {
        /*
        if(processed_measurement_)
        {
            last_state_estimate = statek1k1_;
            last_covar_estimate = covark1k1_;
            //states.push_back(err_statek1k1_);
            //covars.push_back(covark1k1_);

            statekk_ = statek1k1_;
            covarkk_ = covark1k1_;
        }
        else
        {
            last_state_estimate = statek1k_;
            last_covar_estimate = covark1k_;
            //states.push_back(err_statek1k_);
            //covars.push_back(covark1k_);

            statekk_ = statek1k_;
            covarkk_ = covark1k_;
        }
        */
        
        processed_measurement_ = false;

        /*
        statek1k_  = VectorXd::Zero(num_states_);
        statek1k1_ = VectorXd::Zero(num_states_);

        covark1k_  = MatrixXd::Zero(num_states_covar_, num_states_covar_);
        covark1k1_ = MatrixXd::Zero(num_states_covar_, num_states_covar_);
        */
    }

    void MEKF::PrintModelMatrices()
    {
        std::cout << "A:\t" << std::endl << A_ << std::endl << std::endl;
        std::cout << "F:\t" << std::endl << F_ << std::endl << std::endl;
        std::cout << "Q:\t" << std::endl << Q_ << std::endl << std::endl;
        std::cout << "H:\t" << std::endl << H_ << std::endl << std::endl;
        std::cout << "R:\t" << std::endl << R_ << std::endl << std::endl;
    }

} // end namespace