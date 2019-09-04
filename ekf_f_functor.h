/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#ifndef EKF_F_FUNCTOR_H_
#define EKF_F_FUNCTOR_H_

#include <Eigen/Core>
#include <math.h>

#include "ceres/ceres.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Vector3d;

// A templated cost functor for automatic differentiation
class EKF_f_Functor
{
    public: EKF_f_Functor(double dt): dt_(dt) {}

    template <typename T>
    bool operator()(const T* const stateArr,
                    T* residuals) const
    {

        Eigen::Matrix<T, Eigen::Dynamic, 1> statekk(19);
        for (unsigned int idx = 0; idx < 19; idx++)
        { statekk(idx) = stateArr[idx]; }

        Eigen::Matrix<T, Eigen::Dynamic, 1> statek1k = this->KF_NL_f(statekk, dt_);

        for (unsigned int idx = 0; idx < 19; idx++)
        { residuals[idx] = statek1k(idx); }

        return true;
    }


    // function performing non-linear dynamics propagation of state; also `residual` function
    template <typename T> static Eigen::Matrix<T, Eigen::Dynamic, 1> KF_NL_f(Eigen::Matrix<T, Eigen::Dynamic, 1> statekk_, double dt_)
    {
        typedef Eigen::Matrix<T, 3, 1> Vector3T;
        typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXT;
        typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
        typedef Eigen::Quaternion<T> QuaternionT;

        VectorXT statek1k_(statekk_.size());

        // position
        MatrixXT F_pos_ = MatrixXT::Identity(9, 9); // position dynamics_model
        F_pos_(0,3) = static_cast<T>(dt_);
        F_pos_(1,4) = static_cast<T>(dt_);
        F_pos_(2,5) = static_cast<T>(dt_);
        F_pos_(3,6) = static_cast<T>(dt_);
        F_pos_(4,7) = static_cast<T>(dt_);
        F_pos_(5,8) = static_cast<T>(dt_);
        F_pos_(0,6) = static_cast<T>(0.5*pow(dt_,2));
        F_pos_(1,7) = static_cast<T>(0.5*pow(dt_,2));
        F_pos_(2,8) = static_cast<T>(0.5*pow(dt_,2));
        statek1k_.head(9) = F_pos_*statekk_.head(9);

        // orientation
        QuaternionT quatk_;
        quatk_.w() = statekk_(9);
        quatk_.vec() = statekk_.segment(10,3);
        quatk_.normalize();
        QuaternionT dquatk_ =   Eigen::AngleAxis<T>(statekk_(13)*dt_, Vector3T::UnitX())*
                                Eigen::AngleAxis<T>(statekk_(14)*dt_, Vector3T::UnitY())*
                                Eigen::AngleAxis<T>(statekk_(15)*dt_, Vector3T::UnitZ());
        QuaternionT quatk1k_= quatk_*dquatk_;
        statek1k_(9) = quatk1k_.w();
        statek1k_.segment(10,3) = quatk1k_.vec();
        
        MatrixXT F_att_ = F_pos_.block(0,0, 6,6);
        statek1k_.tail(6) = F_att_*statekk_.tail(6);

        return statek1k_;
    }



    private:
        double dt_;
};

#endif // EKF_F_FUNCTOR_H_