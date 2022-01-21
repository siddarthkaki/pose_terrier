/* Copyright (c) 2021 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
#include "MEKF2.h"

namespace MEKF2 {

    MEKF2::MEKF2(const double &dt)
    {
        dt_ = dt;
        processed_measurement_ = false;
    }

    void MEKF2::Init(
        const double &process_noise_std, 
        const double &measurement_noise_std, 
        const double &dt, 
        const double &tau,
        const double &qpsd, 
        const double &max_flip_thresh_deg,  
        const double &pos_uw_threshold,
        const double &pos_uw_pct
    )
    {
        num_att_states_ = 9; // delta_gibbs(3), omega(3), alpha(3)
        num_pos_states_ = 9; // pos(3), posDot(3), posDotDot(3)
        num_states_ = num_att_states_ + num_pos_states_;
        num_att_measurements_ = 3; // delta_gibbs(3)
        num_pos_measurements_ = 3; // pos(3)
        num_measurements_ = num_pos_measurements_ + num_att_measurements_;
        dt_ = dt;
        tau_ = tau;
        max_flip_thresh_deg_ = max_flip_thresh_deg;
        pos_uw_threshold_ = pos_uw_threshold;
        pos_uw_pct_ = pos_uw_pct;
        
        // TODO: spilt process & measurement noise std for pos and att

        Q_ = MatrixXd::Identity(num_states_, num_states_)*pow(process_noise_std,2); // process_noise_covariance
        double qpsd_ = qpsd;
        double pss = qpsd_*tau_/2.0;
        Q_(6,6) *= ( 1.0 - exp(-2.0*dt_/tau_) )*pss;
        Q_(7,7) *= ( 1.0 - exp(-2.0*dt_/tau_) )*pss;
        Q_(8,8) *= ( 1.0 - exp(-2.0*dt_/tau_) )*pss;
        MatrixXd Q_pos_ = MatrixXd::Identity(num_pos_states_, num_pos_states_); // position_process_noise_covariance
        Q_pos_(0,0) = 0.25*pow(dt,4);
        Q_pos_(1,1) = 0.25*pow(dt,4);
        Q_pos_(2,2) = 0.25*pow(dt,4);
        Q_pos_(3,3) = pow(dt,2);
        Q_pos_(4,4) = pow(dt,2);
        Q_pos_(5,5) = pow(dt,2);
        Q_pos_(6,6) = 1.0;
        Q_pos_(7,7) = 1.0;
        Q_pos_(8,8) = 1.0;
        Q_pos_(0,3) = 0.5*pow(dt,3);
        Q_pos_(1,4) = 0.5*pow(dt,3);
        Q_pos_(2,5) = 0.5*pow(dt,3);
        Q_pos_(3,6) = dt;
        Q_pos_(4,7) = dt;
        Q_pos_(5,8) = dt;
        Q_pos_(3,0) = 0.5*pow(dt,3);
        Q_pos_(4,1) = 0.5*pow(dt,3);
        Q_pos_(5,2) = 0.5*pow(dt,3);
        Q_pos_(6,3) = dt;
        Q_pos_(7,4) = dt;
        Q_pos_(8,5) = dt;
        Q_pos_(0,6) =0.5*pow(dt,2);
        Q_pos_(1,7) =0.5*pow(dt,2);
        Q_pos_(2,8) =0.5*pow(dt,2);
        Q_pos_(6,0) =0.5*pow(dt,2);
        Q_pos_(7,1) =0.5*pow(dt,2);
        Q_pos_(8,2) =0.5*pow(dt,2);
        Q_pos_ = Q_pos_*pow(process_noise_std,2);
        Q_.bottomRightCorner(num_pos_states_, num_pos_states_) = Q_pos_;
        
        R_ = MatrixXd::Identity(num_measurements_, num_measurements_)*pow(measurement_noise_std,2); // measurement_noise_covariance

        A_ = MatrixXd::Identity(4, 4); // quaternion_propagation

        F_ = MatrixXd::Identity(num_states_, num_states_); // convariance_dynamics_propagation
        F_.block(3, 3, 3, 3) = I33;
        F_.block(3, 6, 3, 3) = I33 * dt_;
        F_.block(6, 6, 3, 3) = I33 * exp(-dt_ / tau_);
        F_pos_ = MatrixXd::Identity(num_pos_states_, num_pos_states_); // position_dynamics_propagation
        F_pos_(0,3) = dt_;
        F_pos_(1,4) = dt_;
        F_pos_(2,5) = dt_;
        F_pos_(3,6) = dt_;
        F_pos_(4,7) = dt_;
        F_pos_(5,8) = dt_;
        F_pos_(0,6) = 0.5*pow(dt_,2);
        F_pos_(1,7) = 0.5*pow(dt_,2);
        F_pos_(2,8) = 0.5*pow(dt_,2);
        F_.bottomRightCorner(num_pos_states_, num_pos_states_) = F_pos_;

        H_ = MatrixXd::Zero(num_measurements_, num_states_); // measurement_model
        H_.block(0, 0, num_att_measurements_, num_att_measurements_) = Matrix3d::Identity();
        H_.block(num_att_measurements_, num_att_states_, num_pos_measurements_, num_pos_measurements_) = Matrix3d::Identity();

         pos_est_ = Vector3d::Zero();
        quat_est_ = Quaterniond::Identity();
        delta_gibbs_est_ = Vector3d::Zero();
        state_est_ = VectorXd::Zero(num_states_);

        covar_est_ = MatrixXd::Zero(num_states_, num_states_);

        processed_measurement_ = false;
    }

    void MEKF2::Init(const double &process_noise_std, const double &measurement_noise_std, const double &dt)
    {
        num_att_states_ = 9; // delta_gibbs(3), omega(3), alpha(3)
        num_pos_states_ = 9; // pos(3), posDot(3), posDotDot(3)
        num_states_ = num_att_states_ + num_pos_states_;
        num_att_measurements_ = 3; // delta_gibbs(3)
        num_pos_measurements_ = 3; // pos(3)
        num_measurements_ = num_pos_measurements_ + num_att_measurements_;
        dt_ = dt;
        tau_ = 1.0;
        max_flip_thresh_deg_ = 30;
        pos_uw_threshold_ = 1.5;
        pos_uw_pct_ = 0.2;
        // TODO: spilt process & measurement noise std for pos and att

        Q_ = MatrixXd::Identity(num_states_, num_states_)*pow(process_noise_std,2); // process_noise_covariance
        double qpsd = 1e5;
        double pss = qpsd*tau_/2.0;
        Q_(6,6) *= ( 1.0 - exp(-2.0*dt_/tau_) )*pss;
        Q_(7,7) *= ( 1.0 - exp(-2.0*dt_/tau_) )*pss;
        Q_(8,8) *= ( 1.0 - exp(-2.0*dt_/tau_) )*pss;
        MatrixXd Q_pos_ = MatrixXd::Identity(num_pos_states_, num_pos_states_); // position_process_noise_covariance
        Q_pos_(0,0) = 0.25*pow(dt,4);
        Q_pos_(1,1) = 0.25*pow(dt,4);
        Q_pos_(2,2) = 0.25*pow(dt,4);
        Q_pos_(3,3) = pow(dt,2);
        Q_pos_(4,4) = pow(dt,2);
        Q_pos_(5,5) = pow(dt,2);
        Q_pos_(6,6) = 1.0;
        Q_pos_(7,7) = 1.0;
        Q_pos_(8,8) = 1.0;
        Q_pos_(0,3) = 0.5*pow(dt,3);
        Q_pos_(1,4) = 0.5*pow(dt,3);
        Q_pos_(2,5) = 0.5*pow(dt,3);
        Q_pos_(3,6) = dt;
        Q_pos_(4,7) = dt;
        Q_pos_(5,8) = dt;
        Q_pos_(3,0) = 0.5*pow(dt,3);
        Q_pos_(4,1) = 0.5*pow(dt,3);
        Q_pos_(5,2) = 0.5*pow(dt,3);
        Q_pos_(6,3) = dt;
        Q_pos_(7,4) = dt;
        Q_pos_(8,5) = dt;
        Q_pos_(0,6) =0.5*pow(dt,2);
        Q_pos_(1,7) =0.5*pow(dt,2);
        Q_pos_(2,8) =0.5*pow(dt,2);
        Q_pos_(6,0) =0.5*pow(dt,2);
        Q_pos_(7,1) =0.5*pow(dt,2);
        Q_pos_(8,2) =0.5*pow(dt,2);
        Q_pos_ = Q_pos_*pow(process_noise_std,2);
        Q_.bottomRightCorner(num_pos_states_, num_pos_states_) = Q_pos_;
        
        R_ = MatrixXd::Identity(num_measurements_, num_measurements_)*pow(measurement_noise_std,2); // measurement_noise_covariance

        A_ = MatrixXd::Identity(4, 4); // quaternion_propagation

        F_ = MatrixXd::Identity(num_states_, num_states_); // convariance_dynamics_propagation
        F_.block(3, 3, 3, 3) = I33;
        F_.block(3, 6, 3, 3) = I33 * dt_;
        F_.block(6, 6, 3, 3) = I33 * exp(-dt_ / tau_);
        F_pos_ = MatrixXd::Identity(num_pos_states_, num_pos_states_); // position_dynamics_propagation
        F_pos_(0,3) = dt_;
        F_pos_(1,4) = dt_;
        F_pos_(2,5) = dt_;
        F_pos_(3,6) = dt_;
        F_pos_(4,7) = dt_;
        F_pos_(5,8) = dt_;
        F_pos_(0,6) = 0.5*pow(dt_,2);
        F_pos_(1,7) = 0.5*pow(dt_,2);
        F_pos_(2,8) = 0.5*pow(dt_,2);
        F_.bottomRightCorner(num_pos_states_, num_pos_states_) = F_pos_;

        H_ = MatrixXd::Zero(num_measurements_, num_states_); // measurement_model
        H_.block(0, 0, num_att_measurements_, num_att_measurements_) = Matrix3d::Identity();
        H_.block(num_att_measurements_, num_att_states_, num_pos_measurements_, num_pos_measurements_) = Matrix3d::Identity();

        pos_est_ = Vector3d::Zero();
        quat_est_ = Quaterniond::Identity();
        delta_gibbs_est_ = Vector3d::Zero();
        state_est_ = VectorXd::Zero(num_states_);

        covar_est_ = MatrixXd::Zero(num_states_, num_states_);

        processed_measurement_ = false;
    }
    
    void MEKF2::SetInitialStateAndCovar(const Quaterniond &quat0, const Vector3d &omega0, const Vector3d &alpha0, const VectorXd &x0, const MatrixXd &covar0)
    {
         quat_est_ =  quat0;
        state_est_.segment(3, 3) = omega0;
        state_est_.segment(6, 3) = alpha0;
        state_est_.segment(num_att_states_, num_pos_states_) = x0;
        covar_est_ = covar0;
    }

    // Prediction step
    void MEKF2::Predict()
    {
        Vector3d omega_est_ = state_est_.segment(3, 3);

        double omega_norm = omega_est_.norm();
        Vector3d omega_hat = omega_est_ / omega_norm;

        MatrixXd omega_hat_44_equivalent = MatrixXd::Zero(4, 4);
        omega_hat_44_equivalent.block(0, 1, 1, 3) = -omega_hat.transpose();
        omega_hat_44_equivalent.block(1, 0, 3, 1) =  omega_hat;
        omega_hat_44_equivalent.block(1, 1, 3, 3) = -CppRot::CrossProductEquivalent(omega_hat);

        double phi = 0.5 * omega_norm * dt_;

        A_ = cos(phi) * I44 + sin(phi) * omega_hat_44_equivalent;
        F_.block(0, 0, 3, 3) = (-CppRot::CrossProductEquivalent(omega_est_) * dt_).exp();

        // propagate quaternion
        quat_est_ = Utilities::Vec4ToQuat( A_ * Utilities::QuatToVec4(quat_est_) );

        // propagate rest of the state
        state_est_.tail(num_states_-3) = F_.bottomRightCorner(num_states_-3, num_states_-3)*state_est_.tail(num_states_-3);
        pos_est_ = state_est_.segment(num_att_states_, 3);

        //std::cout << pos_est_.transpose() << std::endl << std::endl;
        
        // propagate covariance
        covar_est_ = F_ * covar_est_ * F_.transpose() + Q_;

    }

    // Update step
    void MEKF2::Update(const VectorXd &measurement)
    {
        // Kalman gain
        MatrixXd K = covar_est_ * H_.transpose() * ((H_ * covar_est_ * H_.transpose() + R_).inverse());
        
        Matrix3d pos_covar_est = covar_est_.block(num_att_states_, num_att_states_, 3, 3);
        Matrix3d R_pos = R_.bottomRightCorner(3, 3);
        
        // calculate underweighting factor (Lear's method)
        float alpha_u = sqrt(pos_covar_est.trace());
        float beta_u = 0.0;
        if(alpha_u > pos_uw_threshold_) // TODO: need to tune alpha and beta
        { 
            beta_u = pos_uw_pct_;
            std::cout << "UW ALPHA: " << alpha_u << std::endl;
        }

        Matrix3d K_pos = pos_covar_est * I33.transpose() * (( (1 + beta_u) * (I33 * pos_covar_est * I33.transpose()) + R_pos).inverse());


        // delta Gibbs update
        Quaterniond quat_meas = Utilities::Vec4ToQuat(measurement.head(4));
        Quaterniond delta_quat = CppRot::QuatMult_S(quat_meas, quat_est_.inverse());
        Vector3d meas_innovation = 2.0 * delta_quat.vec() / delta_quat.w();
        delta_gibbs_est_ = K.block(0, 0, 3, 3)*meas_innovation;

        // position update
        pos_est_ = pos_est_ + K_pos*( measurement.tail(3) - I33*pos_est_ );

        //std::cout << pos_est_.transpose() << std::endl << std::endl;

        state_est_.segment(num_att_states_, 3) = pos_est_;
        
        // Joseph update (general)
        MatrixXd I = MatrixXd::Identity(num_states_, num_states_);
        covar_est_ = (I - K * H_) * covar_est_ * ((I - K * H_).transpose()) + K * R_ * (K.transpose());
    }

    // Reset step
    void MEKF2::Reset()
    {
        Quaterniond delta_quat = Quaterniond::Identity();
        delta_quat.w() = 1.0;
        delta_quat.vec() = 0.5*delta_gibbs_est_;
        Quaterniond quat_star = CppRot::QuatMult_S(delta_quat, quat_est_).normalized();
        
        
        // NOTE: heuristic method to ignore 180 deg pose ambiguities
        Quaterniond dq = CppRot::QuatMult_S(quat_est_, quat_star.inverse());
        double dangle = 2.0*acos( abs( dq.w() ) );
        /*MatrixXd att_covar_est = covar_est_.block(0, 0, num_att_states_, num_att_states_);
        double att_covar_rm = sqrt(att_covar_est.trace() / num_att_states_);
        
        std::cout << "ATT COV: "<< att_covar_rm * Utilities::RAD2DEG << std::endl;
        if (att_covar_rm < 10 * Utilities::DEG2RAD)*/
        if (dangle < max_flip_thresh_deg_ * Utilities::DEG2RAD)
        {
            quat_est_ = quat_star;
            //std::cout << "flip heuristic triggered: "<< dangle*Utilities::RAD2DEG << std::endl;
        }
        else{
            std::cout << "flip heuristic triggered: "<< dangle* Utilities::RAD2DEG << std::endl;
        }
        //quat_est_ = quat_star;

        delta_gibbs_est_ = Vector3d::Zero();        
       
        processed_measurement_ = true;
    }

    void MEKF2::StoreAndClean()
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

    void MEKF2::PrintModelMatrices()
    {
        std::cout << "A:\t" << std::endl << A_ << std::endl << std::endl;
        std::cout << "F:\t" << std::endl << F_ << std::endl << std::endl;
        std::cout << "Q:\t" << std::endl << Q_ << std::endl << std::endl;
        std::cout << "H:\t" << std::endl << H_ << std::endl << std::endl;
        std::cout << "R:\t" << std::endl << R_ << std::endl << std::endl;
    }

} // end namespace
