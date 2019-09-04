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

/*
    MEKF::MEKF(const unsigned int &num_states, const unsigned int &num_measurements, const unsigned int &num_inputs, const double &dt)
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
        
        err_statekk_   = VectorXd::Zero(num_states_);
        err_statek1k_  = VectorXd::Zero(num_states_);
        err_statek1k1_ = VectorXd::Zero(num_states_);

        covarkk_   = MatrixXd::Zero(num_states_, num_states_);
        covark1k_  = MatrixXd::Zero(num_states_, num_states_);
        covark1k1_ = MatrixXd::Zero(num_states_, num_states_);

        processed_measurement_ = false;
    }
*/

    void MEKF::Init(const double &process_noise_std, const double &measurement_noise_std, const double &dt)
    {
        num_states_ = 10; // quat(4), omega(3), alpha(3)
        num_states_covar_ = num_states_ - 1; // eul(3), omega(3), alpha(3)
        num_measurements_ = 3; // attitude error (delta_eul)
        num_inputs_ = 0;
        dt_ = dt;

        Q_ = MatrixXd::Identity(num_states_covar_, num_states_covar_)*pow(process_noise_std,2); // process_noise_covariance
        
        R_ = MatrixXd::Identity(num_measurements_, num_measurements_)*pow(measurement_noise_std,2); // measurement_noise_covariance

        F_ = MatrixXd::Identity(num_states_covar_, num_states_covar_); // dynamics_model init; MUST EXTERNALLY UPDATE

        G_ = MatrixXd::Zero(num_states_, num_inputs_); // input_model

        H_ = MatrixXd::Zero(num_measurements_, num_states_covar_); // measurement_model init
        H_.block(0,0,num_measurements_,num_measurements_) = Matrix3d::Identity();

        statekk_   = VectorXd::Zero(num_states_);
        statek1k_  = VectorXd::Zero(num_states_);
        statek1k1_ = VectorXd::Zero(num_states_);

        covarkk_   = MatrixXd::Zero(num_states_covar_, num_states_covar_);
        covark1k_  = MatrixXd::Zero(num_states_covar_, num_states_covar_);
        covark1k1_ = MatrixXd::Zero(num_states_covar_, num_states_covar_);

        processed_measurement_ = false;
    }

    void MEKF::SetInitialStateAndCovar(const VectorXd &state0, const MatrixXd &covar0)
    {
        statekk_ = state0;
        covarkk_ = covar0;

        last_state_estimate = statekk_;
        last_covar_estimate = covarkk_;
        //states.push_back(err_statekk_);
        //covars.push_back(covarkk_);
    }

    // Non-linear prediction step (takes in pointer to function f that performs nonlinear state propagation)
    void MEKF::Predict(const VectorXd &input, VectorXd (*f)(VectorXd, double))
    {
        // TODO automatically determine F_ for time-step k from (*f)
        statek1k_ = (*f)(statekk_, dt_) + G_*input;
        covark1k_ = F_*covarkk_*F_.transpose() + Q_;

        statekk_ = statek1k_;
        covarkk_ = covark1k_;
    }

    // Non-linear update step (takes in pointer to function meas_resid that performs nonlinear measurement update residual computation)
    void MEKF::Update(const VectorXd &measurement, VectorXd (*meas_resid)(const VectorXd&, const VectorXd&, const double&))
    {
        // TODO automatically determine H_ for time-step k from (*h)
        MatrixXd K = covark1k_*H_.transpose()*( ( H_*covark1k_*H_.transpose() + R_ ).inverse() );
        
        VectorXd state_correction = K*( (*meas_resid)(measurement, statek1k_, dt_) );
        
        // implicit reset as part of measurement update for non-attitude states (omega, alpha)
        statek1k1_.tail(6) = statek1k_.tail(6) + state_correction.tail(6);

        delta_eul_angles = state_correction.head(3);
                
        // Joseph update (general)
        MatrixXd I = MatrixXd::Identity(num_states_covar_, num_states_covar_);
        covark1k1_ = (I - K*H_)*covark1k_*((I - K*H_).transpose()) + K*R_*(K.transpose());
    }

    // Reset step
    void MEKF::Reset()
    {
        Quaterniond delta_quat = AngleAxisd(delta_eul_angles(0), Vector3d::UnitX()) *
                                 AngleAxisd(delta_eul_angles(1), Vector3d::UnitY()) *
                                 AngleAxisd(delta_eul_angles(2), Vector3d::UnitZ());

        Quaterniond quatk1k;
        quatk1k.w() = statek1k_(0);
        quatk1k.x() = statek1k_(1);
        quatk1k.y() = statek1k_(2);
        quatk1k.z() = statek1k_(3);

        Quaterniond quatk1k1 = quatk1k * delta_quat;

        statek1k1_(0) = quatk1k1.w();
        statek1k1_(1) = quatk1k1.x();
        statek1k1_(2) = quatk1k1.y();
        statek1k1_(3) = quatk1k1.z();

        delta_eul_angles = Vector3d::Zero();

        processed_measurement_ = true;
    }

    void MEKF::StoreAndClean()
    {
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
        
        processed_measurement_ = false;

        statek1k_  = VectorXd::Zero(num_states_);
        statek1k1_ = VectorXd::Zero(num_states_);

        covark1k_  = MatrixXd::Zero(num_states_covar_, num_states_covar_);
        covark1k1_ = MatrixXd::Zero(num_states_covar_, num_states_covar_);
    }

    VectorXd MEKF::MeasResidFunction(const VectorXd &measurement, const VectorXd &statek1k, const double &dt)
    {
        Quaterniond quat_meas;
        quat_meas.w() = measurement(0);
        quat_meas.x() = measurement(1);
        quat_meas.y() = measurement(2);
        quat_meas.z() = measurement(3);

        Quaterniond quatk1k;
        quatk1k.w() = statek1k(0);
        quatk1k.x() = statek1k(1);
        quatk1k.y() = statek1k(2);
        quatk1k.z() = statek1k(3);

        Quaterniond quat_resid = quat_meas * (quatk1k.inverse());

        Vector3d eul_resid = quat_resid.toRotationMatrix().eulerAngles(0, 1, 2);
        VectorXd eul_resid_return = eul_resid;

        return eul_resid_return;
    }

    void MEKF::PrintModelMatrices()
    {
        std::cout << "F:\t" << std::endl << F_ << std::endl << std::endl;
        std::cout << "Q:\t" << std::endl << Q_ << std::endl << std::endl;
        std::cout << "H:\t" << std::endl << H_ << std::endl << std::endl;
        std::cout << "R:\t" << std::endl << R_ << std::endl << std::endl;
    }

} // end namespace