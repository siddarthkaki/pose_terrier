/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef MEKF_F_FUNCTOR_H_
#define MEKF_F_FUNCTOR_H_

#include <Eigen/Core>
#include <math.h>

#include "ceres/ceres.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;

// A templated cost functor for automatic differentiation
class MEKF_f_Functor
{
    public: MEKF_f_Functor(double dt): dt_(dt) {}

    template <typename T>
    bool operator()(const T* const stateArr,
                    T* residuals) const
    {
        //unsigned int num_states = 9; // eul (3), omega (3), alpha (3)
        unsigned int num_states = 6; // eul (3), omega (3)
        Eigen::Matrix<T, Eigen::Dynamic, 1> statekk(num_states);
        for (unsigned int idx = 0; idx < num_states; idx++)
        { statekk(idx) = stateArr[idx]; }

        Eigen::Matrix<T, Eigen::Dynamic, 1> statek1k = this->MEKF_f(statekk, dt_);

        for (unsigned int idx = 0; idx < num_states; idx++)
        { residuals[idx] = statek1k(idx); }

        return true;
    }


    // function performing non-linear dynamics propagation of state; also `residual` function
    template <typename T> static Eigen::Matrix<T, Eigen::Dynamic, 1> MEKF_f(Eigen::Matrix<T, Eigen::Dynamic, 1> statekk_, double dt_)
    {
        typedef Eigen::Matrix<T, 3, 1> Vector3T;
        typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXT;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
        typedef Eigen::Quaternion<T> QuaternionT;
        typedef Eigen::AngleAxis<T> AngleAxisT;

        VectorXT statek1k_(statekk_.size());

        QuaternionT quatkk_ =   AngleAxisT(statekk_(0), Vector3T::UnitX()) *
                                AngleAxisT(statekk_(1), Vector3T::UnitY()) *
                                AngleAxisT(statekk_(2), Vector3T::UnitZ());

        //TODO: 1st order gauss markov model for omega dynamics /////////////////////////////////////////////////////////////////////////////
        // dynamics_model for omega and alpha
        double tau = 0.075;
        MatrixXT F_ = MatrixXT::Identity(3, 3) * static_cast<T>(exp( - dt_ / tau ));
        //MatrixXT F_ = MatrixXT::Identity(6, 6);
        //F_(0,3) = static_cast<T>(dt_);
        //F_(1,4) = static_cast<T>(dt_);
        //F_(2,5) = static_cast<T>(dt_);

        // propagation
        T delta_eul1 = statekk_(3)*static_cast<T>(dt_);// + statekk_(6)*static_cast<T>( pow(dt_,2)/2.0 );
        T delta_eul2 = statekk_(4)*static_cast<T>(dt_);// + statekk_(7)*static_cast<T>( pow(dt_,2)/2.0 );
        T delta_eul3 = statekk_(5)*static_cast<T>(dt_);// + statekk_(8)*static_cast<T>( pow(dt_,2)/2.0 );

        QuaternionT quat_step = AngleAxisT(delta_eul1, Vector3T::UnitX()) *
                                AngleAxisT(delta_eul2, Vector3T::UnitY()) *
                                AngleAxisT(delta_eul3, Vector3T::UnitZ()); 	

        QuaternionT quatk1k_ = quatkk_ * quat_step;

        Vector3T eulk1k_ = quatk1k_.toRotationMatrix().eulerAngles(0, 1, 2);

        statek1k_(0) = eulk1k_(0);
        statek1k_(1) = eulk1k_(1);
        statek1k_(2) = eulk1k_(2);
        //statek1k_.tail(6) = F_*statekk_.tail(6);
        statek1k_.tail(3) = F_*statekk_.tail(3);

        return statek1k_;
    }

    // function performing non-linear dynamics propagation of state with attitude in quaternion form
    template <typename T> static Eigen::Matrix<T, Eigen::Dynamic, 1> MEKF_f_quat(Eigen::Matrix<T, Eigen::Dynamic, 1> statekk_, double dt_)
    {
        typedef Eigen::Matrix<T, 3, 1> Vector3T;
        typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXT;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
        typedef Eigen::Quaternion<T> QuaternionT;
        typedef Eigen::AngleAxis<T> AngleAxisT;

        // num_states = 10: quat (4), omega (3), alpha (3)
        VectorXT statek1k_(statekk_.size());

        QuaternionT quatkk_;
        quatkk_.w() = statekk_(0);
        quatkk_.x() = statekk_(1);
        quatkk_.y() = statekk_(2);
        quatkk_.z() = statekk_(3);


        // dynamics_model for omega and alpha
        double tau = 0.075;
        MatrixXT F_ = MatrixXT::Identity(3, 3) * static_cast<T>(exp( - dt_ / tau ));
        //MatrixXT F_ = MatrixXT::Identity(6, 6);
        //F_(0,3) = static_cast<T>(dt_);
        //F_(1,4) = static_cast<T>(dt_);
        //F_(2,5) = static_cast<T>(dt_);

        // propagation
        T delta_eul1 = statekk_(4)*static_cast<T>(dt_);// + statekk_(7)*static_cast<T>( pow(dt_,2)/2.0 );
        T delta_eul2 = statekk_(5)*static_cast<T>(dt_);// + statekk_(8)*static_cast<T>( pow(dt_,2)/2.0 );
        T delta_eul3 = statekk_(6)*static_cast<T>(dt_);// + statekk_(9)*static_cast<T>( pow(dt_,2)/2.0 );

        QuaternionT quat_step = AngleAxisT(delta_eul1, Vector3T::UnitX()) *
                                AngleAxisT(delta_eul2, Vector3T::UnitY()) *
                                AngleAxisT(delta_eul3, Vector3T::UnitZ()); 	
        
        QuaternionT quatk1k_ = quatkk_ * quat_step;

        statek1k_(0) = quatk1k_.w();
        statek1k_(1) = quatk1k_.x();
        statek1k_(2) = quatk1k_.y();
        statek1k_(3) = quatk1k_.z();
        //statek1k_.tail(6) = F_*statekk_.tail(6);
        statek1k_.tail(3) = F_*statekk_.tail(3);

        return statek1k_;
    }



    private:
        double dt_;
};

#endif // MEKF_F_FUNCTOR_H_