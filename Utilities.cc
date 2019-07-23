#include "Utilities.h"

using Eigen::AngleAxisd;
using Eigen::Quaterniond;

/**
 * @function Euler2DCM_312
 * @brief Converts euler angles phi, theta, psi (RPY) to DCM with 3-1-2 rotation
 * @return DCM corresponding to 3-1-2 euler angle rotation
 */
Matrix3d Utilities::Euler2DCM_312(const Vector3d& eulVec)
{
    //double phi   = eulVec(0);
    //double theta = eulVec(1);
    //double psi   = eulVec(2);

    Matrix3d DCM;
    DCM =   AngleAxisd( eulVec(1), Vector3d::UnitY() ).toRotationMatrix().transpose()*
            AngleAxisd( eulVec(0), Vector3d::UnitX() ).toRotationMatrix().transpose()*
            AngleAxisd( eulVec(2), Vector3d::UnitZ() ).toRotationMatrix().transpose();

    return DCM;
}

/**
 * @function DCM2Euler_312
 * @brief Converts DCM to euler angles phi, theta, psi (RPY) representing 3-1-2 rotation
 * @return Euler angles corresponding to 3-1-2 rotation
 */
Vector3d Utilities::DCM2Euler_312(const MatrixXd& DCM)
{
    double phi = asin( DCM(1,2) );

    double theta = atan2( -DCM(0,2), DCM(2,2) );

    double psi = atan2( -DCM(1,0), DCM(1,1) );

    Vector3d eulVec;
    eulVec << phi, theta, psi;

    return eulVec;
}

/**
 * @function FeaPointsTargetToChaser
 * @brief Transforms feature points in target frame wrt target to chaser frame wrt chaser
 * @return MatrixXd of feature points in chaser frame wrt chaser
 */
MatrixXd Utilities::FeaPointsTargetToChaser(const Pose& state, const Vector3d& rCamVec, const MatrixXd& rFeaMat)
{
    unsigned int numPts = rFeaMat.rows(); // number of feature points

    MatrixXd rMat_c = MatrixXd::Zero(numPts, 3);

    for (unsigned int idx = 0; idx < numPts; idx++)
    {
        Quaterniond rFeaVeciQuat;
        rFeaVeciQuat.w() = 0;
        rFeaVeciQuat.vec() = rFeaMat.row(idx);
        Vector3d rFeaVeci_c = ( (state.quat)*rFeaVeciQuat*(state.quat.conjugate()) ).vec();

        // position vector of feature point i wrt chaser in chaser frame
        Vector3d rVeci = state.pos - rCamVec + rFeaVeci_c;

        rMat_c.row(idx) = rVeci.transpose();
    }

    return rMat_c;
}

/**
 * @function CameraProjection
 * @brief Performs simple camera projection of 3D point to image plane
 * @return Vector2d of projected 2D point
 */
Vector2d Utilities::CameraProjection(const Vector3d& point3DVec, const double& f)
{
    MatrixXd PMat(3,4);
    PMat << f, 0, 0, 0,
            0, f, 0, 0,
            0, 0, 1, 0;

    VectorXd homoFeaPtVec(4);
    homoFeaPtVec << point3DVec,
                    1;

    Vector3d projVec = PMat*homoFeaPtVec;

    double x = projVec(0)/projVec(2);
    double y = projVec(1)/projVec(2);

    if (!std::isfinite(x)) { x = 0; }
    if (!std::isfinite(y)) { y = 0; }

    Vector2d point2DVec;
    point2DVec << x, y;

    return point2DVec;
}

/**
 * @function SimulateMeasurements
 * @brief Simulates measurements from given pose
 * @return VectorXd of measurements
 */
VectorXd Utilities::SimulateMeasurements(const MatrixXd& rMat, const double& focal_length)
{
    unsigned int numPts = rMat.rows();
   
    // project feature points to image plane
    MatrixXd imgPtMat = MatrixXd::Zero(numPts, 2);
    for (unsigned int idx = 0; idx < numPts; idx++)
    {
        Vector3d rVeci = rMat.row(idx);//.transpose;
        Vector2d imgPti = CameraProjection(rVeci,focal_length);
        imgPtMat.row(idx) = imgPti.transpose();
    }

    // azimuth & elevation measurements for feature points
    VectorXd azVec = VectorXd::Zero(numPts);
    VectorXd elVec = VectorXd::Zero(numPts);
    for (unsigned int idx = 0; idx < numPts; idx++)
    {
        Vector2d imgPti = imgPtMat.row(idx);
        
        azVec(idx) = atan2(imgPti(0),focal_length);
        elVec(idx) = atan2(imgPti(1),focal_length);
    }

    VectorXd yVec(numPts*2); // vector of measurements
    for (unsigned int idx = 0; idx < numPts; idx++)
    {
        yVec(idx*2+0) = azVec(idx);
        yVec(idx*2+1) = elVec(idx);
    }

    return yVec;
}

/**
 * @function AddGaussianNoiseToVector
 * @brief Adds zero-mean Gaussian noise with provided std to VectorXd
 * @return Input VectorXd values with additive Gaussian noise
 */
VectorXd Utilities::AddGaussianNoiseToVector(const VectorXd& vec, const double& std)
{
    const unsigned int numMeas = vec.size();

    MatrixXd covarMat = pow(std,2)*MatrixXd::Identity(numMeas, numMeas);

    Eigen::EigenMultivariateNormal<double> normX_solver(vec, covarMat);
    
    VectorXd vecNoise = normX_solver.samples(1);

    return vecNoise;
}

/**
 * @function UniformRandomAttitude
 * @brief Generates a quaternion representing a uniform random rotation in SO(3)
 * @return Quaterniond contain uniform random rotation
 */
Quaterniond Utilities::UniformRandomAttitude()
{
    Eigen::EigenMultivariateNormal<double> normX_solver(Vector3d::Zero(), MatrixXd::Identity(3,3));
    
    Vector3d rand_axis = normX_solver.samples(1);
    rand_axis.normalize();

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 2.0*M_PI);
    double rand_ang = distribution(generator);

    Quaterniond rand_att;
    rand_att = AngleAxisd(rand_ang, rand_axis);
    rand_att.normalize();

    return rand_att;
}

/**
 * @function ConjugatePose
 * @brief Finds the so-called "conjugate pose" for a given pose
 * @return Pose struct of conjugate pose
 */
Pose Utilities::ConjugatePose(const Pose& state)
{
    /*
    Vector3d p0 = state.quat*Vector3d::Zero();
    Vector3d p1 = state.quat*Vector3d::UnitZ();
    Vector3d p2 = state.quat*Vector3d::UnitY();

    Vector3d p01 = p1 - p0;
    Vector3d p02 = p2 - p0;
    Vector3d n = p01.cross(p02);
    
    if ( n(2) < 0 ) { n = -n; }
    */
   //Vector3d n = state.quat*Vector3d::UnitZ();

    Quaterniond UnitZQuat;
    UnitZQuat.w() = 0.0;
    UnitZQuat.vec() = Vector3d::UnitZ();

    Vector3d n = ( (state.quat)*UnitZQuat*(state.quat.conjugate()) ).vec();

    double phi  = acos( n.dot(Vector3d::UnitZ())/n.norm() );
    Vector3d aHat = n.cross(Vector3d::UnitZ()).normalized();

    Quaterniond conj_quat;
    conj_quat = AngleAxisd( 2*phi, aHat );
    conj_quat = conj_quat.conjugate()*state.quat;
    
    Pose conj_state;
    conj_state.pos = state.pos;
    conj_state.quat = conj_quat;

    return conj_state;
}

/**
 * @function PositionScore
 * @brief Computes the position score for a position state estimate
 * @return Position score value for provided state estimate
 */
double Utilities::PositionScore(const Vector3d& pos, const Vector3d& posHat)
{
    Vector3d posErrVec = pos - posHat;

    //double pos_score = posErrVec.squaredNorm()/pos.squaredNorm();
    double pos_score = posErrVec.norm();
    return pos_score;
}

/**
 * @function AttitudeScore
 * @brief computes the attitude score for an attitude state estimate
 * @return Attitude score value for provided state estimate
 */
double Utilities::AttitudeScore(const Quaterniond& quat, const Quaterniond& quatHat)
{
    Quaterniond dqVec = (quat.normalized())*(quatHat.normalized().conjugate());

    double att_score = 2*acos( std::abs( dqVec.w() ) );
    return att_score;
}