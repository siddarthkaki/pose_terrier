%% housekeeping
clear; close all; clc;

%% init
rot = rotation;

%% read in truth and "measurement" data
truePosesMat = f_read_poses("data/true_poses.csv");

solvedPosesMat = f_read_poses("data/solved_poses.csv");

%% convert euler angle rotation parameterisation to quaternion parameterision
[num_poses,~] = size(truePosesMat);

rMatTruth = truePosesMat(:,1:3);
rMatSolved = solvedPosesMat(:,1:3);

eulMatTruth = truePosesMat(:,4:6);
eulMatSolved = solvedPosesMat(:,4:6);
eulMatFiltered = zeros(num_poses+1,3);

quatMatTruth = zeros(num_poses,4);
quatMatSolved = zeros(num_poses,4);
quatMatFiltered = zeros(num_poses+1,4);

for idx = 1:num_poses,
    eulVecTruth = eulMatTruth(idx,:).';
    eulVecSolved = eulMatSolved(idx,:).';
    
    quatTruth = dcm2quat(rot.euler2dcm_312(eulVecTruth)).';
    quatSolved = dcm2quat(rot.euler2dcm_312(eulVecSolved)).';
    
    quatMatTruth(idx,:) = rot.normedQuat(quatTruth).';
    quatMatSolved(idx,:) = rot.normedQuat(quatSolved).';
end

%% init MEKF
mekf_processed_measurement = false;
mekf_process_noise_std = 0.01;
mekf_measurement_noise_std = 0.05;
mekf_dt = 0.01;
mekf_tau = 1; % 1st-order Gauss-Markov Process for angular velocity


mekf_num_states = 9; % delta_gibbs (3), omega (3), alpha (3)
mekf_num_measurements = 3; % delta_gibbs (3)

% system matrices
mekf_Q = []; % process_noise_covariance
mekf_R = []; % measurement_noise_covariance
mekf_A = []; % quaternion_propagation
mekf_F = []; % dynamics_model_covariance_propagation
mekf_H = []; % measurement_model

% states + covariance matrix
mekf_state_quat = [];
mekf_covar = [];
mekf_state_omega = [];
mekf_state_alpha = [];
mekf_state_delta_gibbs = [];

% process_noise_covariance
mekf_Q = eye(mekf_num_states, mekf_num_states)*(mekf_process_noise_std^2);
% mekf_Q(4,4) = mekf_Q(4,4)*1e2;
% mekf_Q(5,5) = mekf_Q(5,5)*1e2;
% mekf_Q(6,6) = mekf_Q(6,6)*1e2;
mekf_Q(7,7) = mekf_Q(7,7)*1e3;
mekf_Q(8,8) = mekf_Q(8,8)*1e3;
mekf_Q(9,9) = mekf_Q(9,9)*1e3;

mekf_R = eye(mekf_num_measurements, mekf_num_measurements)*(mekf_measurement_noise_std^2); % measurement_noise_covariance

mekf_A = eye(4, 4); % quaternion_propagation

mekf_F = eye(mekf_num_states, mekf_num_states); % dynamics_model_covariance_propagation

mekf_H = zeros(mekf_num_measurements, mekf_num_states); % measurement_model init
mekf_H(1:3, 1:3) = eye(3);

% init states + covariance matrix
mekf_state_quat = quatMatSolved(1,:).';%eul2quat(rand(1,3)).';
mekf_state_omega = 0.01*ones(3,1);
mekf_state_alpha = 0.1*ones(3,1);
mekf_state_delta_gibbs = zeros(3,1);

mekf_covar = eye(mekf_num_states, mekf_num_states);
mekf_covar(1,1) = 0.1;
mekf_covar(2,2) = 0.1;
mekf_covar(3,3) = 0.1;

quatMatFiltered(1,:) = mekf_state_quat.';
eulMatFiltered(1,:) = rot.dcm2euler_312(rot.quat2Tmat(mekf_state_quat)).';


%% loop

for idx = 1:num_poses,
    %% propagation (state propagation in terms of quaternions, covariance propagation in terms of gibbs vector)
    Omega = norm(mekf_state_omega);
    omega_hat = mekf_state_omega / Omega;
    I44 = eye(4,4);
    omega_hat_44_equivalent = zeros(4,4);
    omega_hat_44_equivalent(1, 2:4) = -omega_hat.';
    omega_hat_44_equivalent(2:4, 1) = omega_hat;
    omega_hat_44_equivalent(2:4, 2:4) = -rot.crossProductEquivalent(omega_hat);
    
    phi = 0.5*Omega*mekf_dt;
    
    mekf_A = cos(phi)*I44 + sin(phi)*omega_hat_44_equivalent;
    mekf_F(1:3, 1:3) = expm(-rot.crossProductEquivalent(mekf_state_omega)*mekf_dt);
    mekf_F(4:6, 4:6) = eye(3);
    mekf_F(4:6, 7:9) = eye(3)*mekf_dt;
    mekf_F(7:9, 7:9) = eye(3)*exp(-mekf_dt / mekf_tau); 
    
    % propagate quaternion
    mekf_state_quat = mekf_A*mekf_state_quat;
    
    % propagate omega vector
    mekf_state_omega = mekf_F(4:6, 4:9)*[mekf_state_omega; mekf_state_alpha];
    
    % propagate alpha vector
    mekf_state_alpha = mekf_F(7:9, 7:9)*mekf_state_alpha;
    
    % propagate covariance
    mekf_covar = mekf_F*mekf_covar*mekf_F.' + mekf_Q;
    
    %% measurement update (updates delta gibbs vector, and omega vector)
    mekf_K = (mekf_covar * mekf_H.') / (mekf_H * mekf_covar * mekf_H.' + mekf_R);
    
    % delta Gibbs update
    mekf_quat_meas = quatMatSolved(idx,:).';
    mekf_delta_quat = rot.quatmult_S(mekf_quat_meas, quatinv(mekf_state_quat.').');
    mekf_meas_innovation = 2 * mekf_delta_quat(2:4) / mekf_delta_quat(1);
    mekf_state_delta_gibbs = mekf_K(1:3,1:3)*mekf_meas_innovation;
    
    % Joseph covariance update (general)
    mekf_I = eye(mekf_num_states, mekf_num_states);
    mekf_covar = (mekf_I -mekf_K * mekf_H) * mekf_covar * ((mekf_I - mekf_K * mekf_H).') + mekf_K * mekf_R * (mekf_K.');
    
    %% reset
    mekf_quat_star_temp = rot.quatmult_S([1; 0.5*mekf_state_delta_gibbs], mekf_state_quat);
    mekf_state_quat_temp = mekf_quat_star_temp / norm(mekf_quat_star_temp);
    
    % NOTE: heuristic method to ignore 180 deg pose ambiguities
    dquat = rot.quatmult_S( rot.normedQuat(mekf_state_quat), quatconj(rot.normedQuat(mekf_state_quat_temp).').' );
    dangle = 2*acos( abs( dquat(1) ) );
    if dangle < deg2rad(20) && idx > 10,
        mekf_state_quat = mekf_state_quat_temp;
    end
    
    mekf_state_delta_gibbs = zeros(3,1);
    
    %% storage
    quatMatFiltered(idx+1,:) = mekf_state_quat.';
    eulMatFiltered(idx+1,:) = rot.dcm2euler_312(rot.quat2Tmat(mekf_state_quat)).';
end

%% post-analysis
% compute attitude score
attScoreVec         = 1000*ones(num_poses,1);
attScoreVecFiltered = 1000*ones(num_poses,1);

for idx = 1:num_poses,

    quat            = quatMatTruth(idx,:);
    quatHat         = quatMatSolved(idx,:);
    quatHatFiltered = quatMatFiltered(idx+1,:);
    
    dquat         = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHat)) );
    dquatFiltered = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHatFiltered)) );
    %dqVec = (quat.normalized())*(quatHat.normalized().conjugate());

    attScoreVec(idx)         = 2*acos( abs( dquat(1) ) ); % rad
    attScoreVecFiltered(idx) = 2*acos( abs( dquatFiltered(1) ) ); % rad
end

%% plotting attitude
tVec = mekf_dt:mekf_dt:mekf_dt*num_poses;
tVecF = 0:mekf_dt:mekf_dt*num_poses;

figure(2)
subplot(4,1,1)
plot(tVec, rad2deg(eulMatTruth(:,1)))
title('Attitude')
hold on
plot(tVec, rad2deg(eulMatSolved(:,1)))
plot(tVecF, rad2deg(eulMatFiltered(:,1)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','NLS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\phi [deg]')

subplot(4,1,2)
plot(tVec, rad2deg(eulMatTruth(:,2)))
hold on
plot(tVec, rad2deg(eulMatSolved(:,2)))
plot(tVecF, rad2deg(eulMatFiltered(:,2)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','NLS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\theta [deg]')

subplot(4,1,3)
plot(tVec, rad2deg(eulMatTruth(:,3)))
hold on
plot(tVec, rad2deg(eulMatSolved(:,3)))
plot(tVecF, rad2deg(eulMatFiltered(:,3)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','NLS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\psi [deg]')

subplot(4,1,4)
plot(tVec, rad2deg(attScoreVec),'color',[0.8500, 0.3250, 0.0980])
hold on
plot(tVec, rad2deg(attScoreVecFiltered),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('NLS','KF')
xlabel('time [s]')
ylabel('att\_score [deg]')