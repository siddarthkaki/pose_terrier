#ifndef CPPROT_H_
#define CPPROT_H_

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <vector>
#include <iostream>
#include <math.h>

using Eigen::Vector3d;
using Eigen::VectorXd;
using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Quaterniond;

class CppRot
{
    public:
        /**
         * @function CrossProductEquivalent
         * 
         * @brief   Outputs the cross-product-equivalent matrix u_cross such that 
         *          for arbitrary 3-by-1 vectors u and v: cross(u,v) = u_cross*v.
         * 
         * @return  3-by-3 skew-symmetric cross-product equivalent matrix
         */
        static inline Matrix3d CrossProductEquivalent(const Vector3d &u)
        {
            Matrix3d u_cross;
            u_cross <<     0, -u(2),  u(1),
                        u(2),     0, -u(0),
                       -u(1),  u(0),     0;

            return u_cross;
        };

        /**
         * @function QuatMult_S
         * 
         * @brief   Computes q_return = q1 * q2, where * is the quaternion 
         *          product defined by Malcolm Schuster, which is the 
         *          composition formula for rotating-space successive 
         *          rotations. q2 represents the first rotation applied, 
         *          and q1 represents the second rotation applied.
         * 
         * @reference   Zanetti, Renato. "Rotations, Transformations, Left 
         *              Quaternions, Right Quaternions?." The Journal of the 
         *              Astronautical Sciences (2019): 1-21.
         * 
         * @return  Equivalent quaternion representation for the composite 
         *          rotating-space rotation: q_return = q1 * q2
         */
        static inline Quaterniond QuatMult_S(const Quaterniond &q1, const Quaterniond &q2)
        {
            double q1w = q1.w(); // scalar
            double q2w = q2.w(); // scalar
            
            Vector3d q1vec = q1.vec(); // vector
            Vector3d q2vec = q2.vec(); // vector
            
            Quaterniond q_return;

            double q3w = q1w*q2w - q1vec.dot(q2vec);
            Vector3d q3vec = q1w*q2vec + q2w*q1vec - q1vec.cross(q2vec);

            q_return.w() = q3w;
            q_return.vec() = q3vec;

            return q_return;
        };

        /**
         * @function QuatMult_H
         * 
         * @brief   Computes q_return = q2 * q1, where * is the quaternion 
         *          product defined by William Hamilton, which is the 
         *          composition formula for fixed-space successive 
         *          rotations. q1 represents the first rotation applied, 
         *          and q2 represents the second rotation applied.
         * 
         * @reference   Zanetti, Renato. "Rotations, Transformations, Left 
         *              Quaternions, Right Quaternions?." The Journal of the 
         *              Astronautical Sciences (2019): 1-21.
         * 
         * @return  Equivalent quaternion representation for the composite 
         *          fixed-space rotation: q_return = q2 * q1
         */
        static inline Quaterniond QuatMult_H(const Quaterniond &q2, const Quaterniond &q1)
        {            
            double q1w = q1.w(); // scalar
            double q2w = q2.w(); // scalar
            
            Vector3d q1vec = q1.vec(); // vector
            Vector3d q2vec = q2.vec(); // vector
            
            Quaterniond q_return;

            double q3w = q1w*q2w - q2vec.dot(q1vec);
            Vector3d q3vec = q1w*q2vec + q2w*q1vec + q2vec.cross(q1vec);

            q_return.w() = q3w;
            q_return.vec() = q3vec;

            return q_return;
        };

        /**
         * @function Quat2Tmat
         * 
         * @brief   Computes the (passive) transformation matrix represented by
         *          the input quaternion
         * 
         * @reference   Zanetti, Renato. "Rotations, Transformations, Left 
         *              Quaternions, Right Quaternions?." The Journal of the 
         *              Astronautical Sciences (2019): 1-21.
         * 
         * @return  3-by-3 (passive) transformation matrix
         */
        static inline Matrix3d Quat2Tmat(const Quaterniond &quat)
        {
            double qw = quat.w();
            Vector3d qvec = quat.vec();
            MatrixXd qcross = CrossProductEquivalent(qvec);
            Matrix3d Tmat = Matrix3d::Identity() - 2*qw*qcross + 2*qcross*qcross;

            return Tmat;
        };

        /**
         * @function Quat2Rmat
         * 
         * @brief   Computes the (active) rotation matrix represented by
         *          the input quaternion
         * 
         * @reference   Zanetti, Renato. "Rotations, Transformations, Left 
         *              Quaternions, Right Quaternions?." The Journal of the 
         *              Astronautical Sciences (2019): 1-21.
         * 
         * @return  3-by-3 (active) rotation matrix
         */
        static inline Matrix3d Quat2Rmat(const Quaterniond &quat)
        {
            double qw = quat.w();
            Vector3d qvec = quat.vec();
            MatrixXd qcross = CrossProductEquivalent(qvec);
            Matrix3d Rmat = Matrix3d::Identity() + 2*qw*qcross + 2*qcross*qcross;

            return Rmat;
        };
        
        /**
         * @function AngleAxis2Quat
         * 
         * @brief   Computes the quaternion associated with the (active)
         *          axis-angle rotation of angle `theta` (in rad) about unit
         *          vector `axis`
         * 
         * @reference   Zanetti, Renato. "Rotations, Transformations, Left 
         *              Quaternions, Right Quaternions?." The Journal of the 
         *              Astronautical Sciences (2019): 1-21.
         * 
         * @return  Quaternion representing input (active) axis-angle rotation
         */
        static inline Quaterniond AngleAxis2Quat(const double &theta, const Vector3d &axis)
        {
            Quaterniond quat;
            quat.w() = cos(theta/2);
            quat.vec() = sin(theta/2)*axis.normalized();
            
            return quat;
        };

        /**
         * @function AngleAxis2Tmat
         * 
         * @brief   Computes the (passive) transformation matrix associated
         *          with the (active) axis-angle rotation of angle `theta`
         *          (in rad) about unit vector `axis`
         * 
         * @reference   Zanetti, Renato. "Rotations, Transformations, Left 
         *              Quaternions, Right Quaternions?." The Journal of the 
         *              Astronautical Sciences (2019): 1-21.
         * 
         * @return  3-by-3 (passive) transformation matrix representing input
         *          (active) axis-angle rotation
         */
        static inline Matrix3d AngleAxis2Tmat(const double &theta, const Vector3d &axis)
        {
            Matrix3d Tmat = Quat2Tmat(AngleAxis2Quat(theta, axis));

            return Tmat;
        };

        /**
         * @function AngleAxis2Rmat
         * 
         * @brief   Computes the (active) rotation matrix associated
         *          with the (active) axis-angle rotation of angle `theta`
         *          (in rad) about unit vector `axis`
         * 
         * @reference   Zanetti, Renato. "Rotations, Transformations, Left 
         *              Quaternions, Right Quaternions?." The Journal of the 
         *              Astronautical Sciences (2019): 1-21.
         * 
         * @return  3-by-3 (active) rotation matrix representing input
         *          (active) axis-angle rotation
         */
        static inline Matrix3d AngleAxis2Rmat(const double &theta, const Vector3d &axis)
        {
            Matrix3d Rmat = Quat2Rmat(AngleAxis2Quat(theta, axis));

            return Rmat;
        };
};

#endif // CPPROT_H_