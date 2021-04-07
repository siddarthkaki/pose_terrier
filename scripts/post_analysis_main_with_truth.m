%% housekeeping
clear; close all; clc;

%% init
prefix = "../data/" + "1617822191" + "_";

tVec = f_read_timestamps(prefix + "timestamps.csv");

truePosesMat = f_read_poses("../data/true_poses.csv");

solvedPosesMat = f_read_poses(prefix + "solved_poses.csv");

filteredPosesMat = f_read_poses(prefix + "filtered_poses.csv");

%statesMat = f_read_covars(prefix + "kf_states.csv");

%covarsMat = f_read_covars(prefix + "kf_covars.csv");

[num_poses,~] = size(filteredPosesMat);

tVecMeas = f_read_timestamps("../data/meas_timestamps.csv");

%% truth interpolation
truePosesMatInterp = zeros(size(filteredPosesMat));
for idx = 1:6,
    v = truePosesMat(:,idx);
    vq1 = interp1(tVecMeas,v,tVec);
    truePosesMatInterp(:,idx) = vq1;
end

%% compute position score

posScoreVec         = NaN*ones(num_poses,1);
posScoreVecFiltered = NaN*ones(num_poses,1);

for idx = 1:num_poses,

    posErr         = truePosesMatInterp(idx,1:3) - solvedPosesMat(idx,1:3);
    posErrFiltered = truePosesMatInterp(idx,1:3) - filteredPosesMat(idx,1:3);

    posScoreVec(idx)         = norm(posErr);%/norm(truePosesMat(idx,1:3));
    posScoreVecFiltered(idx) = norm(posErrFiltered);%/norm(truePosesMat(idx,1:3));
end

%% compute attitude score

attScoreVec         = 1000*ones(num_poses,1);
attScoreVecFiltered = 1000*ones(num_poses,1);

for idx = 1:num_poses,

    quat            = angle2quat(truePosesMatInterp(idx,4), truePosesMatInterp(idx,5), truePosesMatInterp(idx,6) );
    quatHat         = angle2quat(    solvedPosesMat(idx,4),     solvedPosesMat(idx,5),     solvedPosesMat(idx,6) );
    quatHatFiltered = angle2quat(  filteredPosesMat(idx,4),   filteredPosesMat(idx,5),   filteredPosesMat(idx,6) );
    
    dquat         = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHat)) );
    dquatFiltered = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHatFiltered)) );
    %dqVec = (quat.normalized())*(quatHat.normalized().conjugate());

    attScoreVec(idx)         = 2*acos( abs( dquat(1) ) ); % rad
    attScoreVecFiltered(idx) = 2*acos( abs( dquatFiltered(1) ) ); % rad
end

%% plotting position
figure('DefaultAxesFontSize',15)
subplot(4,1,1)
plot(tVecMeas, truePosesMat(:,1))
title('Position')
hold on
plot(tVec, solvedPosesMat(:,1))
plot(tVec, filteredPosesMat(:,1),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','NLS','KF','Location','Southwest')
xlabel('time [s]')
ylabel('x [m]')

subplot(4,1,2)
plot(tVecMeas, truePosesMat(:,2))
hold on
plot(tVec, solvedPosesMat(:,2))
plot(tVec, filteredPosesMat(:,2),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','NLS','KF','Location','Northeast')
xlabel('time [s]')
ylabel('y [m]')

subplot(4,1,3)
plot(tVecMeas, truePosesMat(:,3))
hold on
plot(tVec, solvedPosesMat(:,3))
plot(tVec, filteredPosesMat(:,3),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','NLS','KF','Location','Southwest')
xlabel('time [s]')
ylabel('z [m]')
%ylim([0,50])

subplot(4,1,4)
plot(tVec, posScoreVec,'color',[0.8500, 0.3250, 0.0980])
hold on
plot(tVec, posScoreVecFiltered,'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('NLS','KF')
xlabel('time [s]')
ylabel('pos\_score [m]')

boldify;

%% plotting attitude
figure('DefaultAxesFontSize',15)
subplot(4,1,1)
plot(tVecMeas, rad2deg(truePosesMat(:,4)))
title('Attitude')
hold on
plot(tVec, rad2deg(solvedPosesMat(:,4)))
plot(tVec, rad2deg(filteredPosesMat(:,4)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','NLS','KF','Location','Southeast')
xlabel('time [s]')
ylabel('\phi [deg]')
%ylim( [min(rad2deg(truePosesMat(:,4)))-5, max(rad2deg(truePosesMat(:,4)))+5] )

subplot(4,1,2)
plot(tVecMeas, rad2deg(truePosesMat(:,5)))
hold on
plot(tVec, rad2deg(solvedPosesMat(:,5)))
plot(tVec, rad2deg(filteredPosesMat(:,5)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','NLS','KF','Location','Southeast')
xlabel('time [s]')
ylabel('\theta [deg]')

subplot(4,1,3)
plot(tVecMeas, rad2deg(truePosesMat(:,6)))
hold on
plot(tVec, rad2deg(solvedPosesMat(:,6)))
plot(tVec, rad2deg(filteredPosesMat(:,6)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','NLS','KF','Location','Southeast')
xlabel('time [s]')
ylabel('\psi [deg]')

subplot(4,1,4)
plot(tVec, rad2deg(attScoreVec),'color',[0.8500, 0.3250, 0.0980])
hold on
plot(tVec, rad2deg(attScoreVecFiltered),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
ylim([0 25])
legend('NLS','KF')
xlabel('time [s]')
ylabel('att\_score [deg]')

boldify;