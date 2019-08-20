/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */
 
#include "KalmanFilter.h"

namespace KF {

    KalmanFilter::KalmanFilter()
    {
        processed_measurement_ = false;
    }

    KalmanFilter::KalmanFilter(const unsigned int &num_states, const unsigned int &num_measurements, const unsigned int &num_inputs, const double &dt)
    {
        num_states_ = num_states;
        num_measurements_ = num_measurements;
        num_inputs_ = num_inputs;
        dt_ = dt;

        Q_ = MatrixXd::Identity(num_states_, num_states_); // process_noise_covariance
        R_ = MatrixXd::Identity(num_measurements_, num_measurements_); // measurement_noise_covariance
        F_ = MatrixXd::Identity(num_states_, num_states_); // dynamics_model
        G_ = MatrixXd::Identity(num_states_, num_inputs_);; // input_model
        H_ = MatrixXd::Identity(num_measurements_, num_states_); // measurement_model
        
        statekk_   = VectorXd::Zero(num_states_);
        statek1k_  = VectorXd::Zero(num_states_);
        statek1k1_ = VectorXd::Zero(num_states_);

        covarkk_   = MatrixXd::Zero(num_states_, num_states_);
        covark1k_  = MatrixXd::Zero(num_states_, num_states_);
        covark1k1_ = MatrixXd::Zero(num_states_, num_states_);

        processed_measurement_ = false;
    }

    // Assumes a discrete Wiener process acceleration model
    void KalmanFilter::InitLinearPoseTracking(const double &process_noise_std, const double &measurement_noise_std, const double &dt)
    {
        num_states_ = 18; // pos (x, y, z), posDot, posDotDot, eul (phi, theta, psi), eulDot, eulDotDot
        num_measurements_ = 6; // pose (x, y, z, phi, theta, psi)
        num_inputs_ = 0;
        dt_ = dt;

        Q_ = MatrixXd::Identity(num_states_, num_states_);//*pow(process_noise_std,2); // process_noise_covariance
        Q_(0,0) = 0.25*pow(dt,4);
        Q_(1,1) = 0.25*pow(dt,4);
        Q_(2,2) = 0.25*pow(dt,4);
        Q_(3,3) = pow(dt,2);
        Q_(4,4) = pow(dt,2);
        Q_(5,5) = pow(dt,2);
        Q_(6,6) = 1.0;
        Q_(7,7) = 1.0;
        Q_(8,8) = 1.0;
        Q_(0,3) = 0.5*pow(dt,3);
        Q_(1,4) = 0.5*pow(dt,3);
        Q_(2,5) = 0.5*pow(dt,3);
        Q_(3,6) = dt;
        Q_(4,7) = dt;
        Q_(5,8) = dt;
        Q_(3,0) = 0.5*pow(dt,3);
        Q_(4,1) = 0.5*pow(dt,3);
        Q_(5,2) = 0.5*pow(dt,3);
        Q_(6,3) = dt;
        Q_(7,4) = dt;
        Q_(8,5) = dt;
        Q_(0,6) =0.5*pow(dt,2);
        Q_(1,7) =0.5*pow(dt,2);
        Q_(2,8) =0.5*pow(dt,2);
        Q_(6,0) =0.5*pow(dt,2);
        Q_(7,1) =0.5*pow(dt,2);
        Q_(8,2) =0.5*pow(dt,2);
        Q_.block(9,9,9,9) = Q_.block(0,0,9,9);
        Q_ = Q_*pow(process_noise_std,2);
        
        R_ = MatrixXd::Identity(num_measurements_, num_measurements_)*pow(measurement_noise_std,2); // measurement_noise_covariance

        F_ = MatrixXd::Identity(num_states_, num_states_); // dynamics_model
        // position
        F_(0,3) = dt_;
        F_(1,4) = dt_;
        F_(2,5) = dt_;
        F_(3,6) = dt_;
        F_(4,7) = dt_;
        F_(5,8) = dt_;
        F_(0,6) = 0.5*pow(dt_,2);
        F_(1,7) = 0.5*pow(dt_,2);
        F_(2,8) = 0.5*pow(dt_,2);
        // orientation
        F_( 9,12) = dt_;
        F_(10,13) = dt_;
        F_(11,14) = dt_;
        F_(12,15) = dt_;
        F_(13,16) = dt_;
        F_(14,17) = dt_;
        F_( 9,15) = 0.5*pow(dt_,2);
        F_(10,16) = 0.5*pow(dt_,2);
        F_(11,17) = 0.5*pow(dt_,2);

        G_ = MatrixXd::Zero(num_states_, num_inputs_); // input_model

        H_ = MatrixXd::Zero(num_measurements_, num_states_); // measurement_model
        H_(0,0) = 1;  // x
        H_(1,1) = 1;  // y
        H_(2,2) = 1;  // z
        H_(3,9) = 1;  // roll (phi)
        H_(4,10) = 1; // pitch (theta)
        H_(5,11) = 1; // yaw (psi)

        statekk_   = VectorXd::Zero(num_states_);
        statek1k_  = VectorXd::Zero(num_states_);
        statek1k1_ = VectorXd::Zero(num_states_);

        covarkk_   = MatrixXd::Zero(num_states_, num_states_);
        covark1k_  = MatrixXd::Zero(num_states_, num_states_);
        covark1k1_ = MatrixXd::Zero(num_states_, num_states_);

        processed_measurement_ = false;
    }

    void KalmanFilter::InitNonLinearPoseTracking(const double &process_noise_std, const double &measurement_noise_std, const double &dt)
    {
        num_states_ = 19; // pos (x,y,z), posDot, posDotDot, quat(4), omega(3), omegaDot(3) 
        num_measurements_ = 7; // Pose (x, y, z, quat)
        num_inputs_ = 0;
        dt_ = dt;

        Q_ = MatrixXd::Identity(num_states_, num_states_)*pow(process_noise_std,2); // process_noise_covariance
        
        R_ = MatrixXd::Identity(num_measurements_, num_measurements_)*pow(measurement_noise_std,2); // measurement_noise_covariance

        F_ = MatrixXd::Identity(num_states_, num_states_); // dynamics_model init; MUST EXTERNALLY UPDATE

        G_ = MatrixXd::Zero(num_states_, num_inputs_); // input_model

        H_ = MatrixXd::Zero(num_measurements_, num_states_); // measurement_model init; MUST EXTERNALLY UPDATE

        statekk_   = VectorXd::Zero(num_states_);
        statek1k_  = VectorXd::Zero(num_states_);
        statek1k1_ = VectorXd::Zero(num_states_);

        covarkk_   = MatrixXd::Zero(num_states_, num_states_);
        covark1k_  = MatrixXd::Zero(num_states_, num_states_);
        covark1k1_ = MatrixXd::Zero(num_states_, num_states_);

        processed_measurement_ = false;
    }

    void KalmanFilter::SetInitialStateAndCovar(const VectorXd &state0, MatrixXd &covar0)
    {
        statekk_ = state0;
        covarkk_ = covar0;

        last_state_estimate = statekk_;
        last_covar_estimate = covarkk_;
        //states.push_back(statekk_);
        //covars.push_back(covarkk_);
    }

    // Linear prediction step
    void KalmanFilter::Predict(const VectorXd &input)
    {
        statek1k_ = F_*statekk_ + G_*input;
        covark1k_ = F_*covarkk_*F_.transpose() + Q_;

        statekk_ = statek1k_;
        covarkk_ = covark1k_;
    }

    // Non-linear prediction step (takes in pointer to function f that performs nonlinear state propagation)
    void KalmanFilter::Predict(const VectorXd &input, VectorXd (*f)(VectorXd, double))
    {
        // TODO automatically determine F_ for time-step k from (*f)
        statek1k_ = (*f)(statekk_, dt_) + G_*input;
        covark1k_ = F_*covarkk_*F_.transpose() + Q_;

        statekk_ = statek1k_;
        covarkk_ = covark1k_;
    }

    // Linear update step
    void KalmanFilter::Update(const VectorXd &measurement)
    {
        MatrixXd K = covark1k_*H_.transpose()*( ( H_*covark1k_*H_.transpose() + R_ ).inverse() );
        
        statek1k1_ = statek1k_ + K*( measurement - H_*statek1k_ );
        //statek1k1_ = statek1k_ + covark1k_*H_.transpose()*( ( H_*covark1k_*H_.transpose() + R_ ).colPivHouseholderQr().solve(measurement - H_*statek1k_) );

        // Update for linear Gaussian systems
        covark1k1_ = covark1k_ - K*H_*covark1k_;

        processed_measurement_ = true;
    }

    // Non-linear update step (takes in pointer to function h that performs nonlinear measurement update)
    void KalmanFilter::Update(const VectorXd &measurement, VectorXd (*h)(VectorXd, double))
    {
        // TODO automatically determine H_ for time-step k from (*h)
        MatrixXd K = covark1k_*H_.transpose()*( ( H_*covark1k_*H_.transpose() + R_ ).inverse() );
        
        statek1k1_ = statek1k_ + K*( measurement - (*h)(statek1k_, dt_) );
                
        // Joseph update (general)
        MatrixXd I = MatrixXd::Identity(num_states_, num_states_);
        covark1k1_ = (I - K*H_)*covark1k_*((I - K*H_).transpose()) + K*R_*(K.transpose());

        processed_measurement_ = true;
    }

    void KalmanFilter::KFStep(const VectorXd &measurement)
    {
        KFStep(measurement, VectorXd::Zero(num_inputs_));
    }

    void KalmanFilter::KFStep(const VectorXd &measurement, const VectorXd &input)
    {
        Predict(input);
        Update(measurement);

        StoreAndClean();
    }

    void KalmanFilter::StoreAndClean()
    {
        if(processed_measurement_)
        {
            last_state_estimate = statek1k1_;
            last_covar_estimate = covark1k1_;
            //states.push_back(statek1k1_);
            //covars.push_back(covark1k1_);

            statekk_ = statek1k1_;
            covarkk_ = covark1k1_;
        }
        else
        {
            last_state_estimate = statek1k_;
            last_covar_estimate = covark1k_;
            //states.push_back(statek1k_);
            //covars.push_back(covark1k_);

            statekk_ = statek1k_;
            covarkk_ = covark1k_;
        }
        
        processed_measurement_ = false;

        statek1k_  = VectorXd::Zero(num_states_);
        statek1k1_ = VectorXd::Zero(num_states_);

        covark1k_  = MatrixXd::Zero(num_states_, num_states_);
        covark1k1_ = MatrixXd::Zero(num_states_, num_states_);
    }

    void KalmanFilter::PrintModelMatrices()
    {
        std::cout << "F:\t" << std::endl << F_ << std::endl << std::endl;
        std::cout << "Q:\t" << std::endl << Q_ << std::endl << std::endl;
        std::cout << "H:\t" << std::endl << H_ << std::endl << std::endl;
        std::cout << "R:\t" << std::endl << R_ << std::endl << std::endl;
    }

} // end namespace