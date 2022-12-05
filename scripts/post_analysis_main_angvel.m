%% housekeeping
clear; close all; clc;

%% init

pose_type = 'NLS'; % 'PnP'

% folder = "long_test/";
% folder = "short_test/";
% folder = "new_test/";
folder = "";

% prefix = "../data/long_test/" + "1661875312" + "_";
% prefix = "../data/short_test/" + "1661872620" + "_";
% prefix = "../data/new_test/" + "1661975852" + "_";
prefix = "../data/" + "1663343139" + "_";


tVec = f_read_timestamps(prefix + "timestamps.csv");

truePosesMat = f_read_poses("../data/" + folder + "true_poses.csv");

solvedPosesMat = f_read_poses(prefix + "solved_poses.csv");

filteredPosesMat = f_read_poses(prefix + "filtered_poses.csv");

trueAngVelMat = f_read_omegas("../data/" + folder + "true_omegas.csv");

solvedAngVelMat = f_read_omegas(prefix + "solved_omegas.csv");

filteredAngVelMat = f_read_omegas(prefix + "filtered_omegas.csv");

filteredCovarsDiagMat = f_read_covars(prefix + "filtered_covar_diag.csv");

%statesMat = f_read_covars(prefix + "kf_states.csv");

%covarsMat = f_read_covars(prefix + "kf_covars.csv");

[num_poses,~] = size(filteredPosesMat);

tVecMeas = f_read_timestamps("../data/" + folder + "meas_timestamps.csv");

windowSize = f_read_timestamps(prefix + folder + "window_size.csv");

epsHistory = f_read_timestamps(prefix + folder + "eps_history.csv");


%% temp fix
% [num_true_poses,~] = size(truePosesMat);
% truePosesMat = [truePosesMat; zeros(num_poses - num_true_poses,6)];
% tVecMeas = 0:0.04:num_poses-1;

truePosesMat(:,5) = wrapTo2Pi(truePosesMat(:,5));
solvedPosesMat(:,5) = wrapTo2Pi(solvedPosesMat(:,5));
filteredPosesMat(:,5) = wrapTo2Pi(filteredPosesMat(:,5));

%% covar -> std extraction
attStdMat = sqrt(filteredCovarsDiagMat(:,1:3));
posStdMat = sqrt(filteredCovarsDiagMat(:,10:12));

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

%%
epsmean = mean(epsHistory(100:end))
epsstd = std(epsHistory(100:end))

%% TEST plotting 3sigmas
f100 = figure(100);
subplot(2,1,1)
plot(tVec, 3*rad2deg(attStdMat));
grid on
ylabel('deg')
ylim([0 10])
legend('\phi', '\theta', '\psi')

subplot(2,1,2)
plot(tVec, 3*posStdMat);
grid on
ylabel('m')
legend('x', 'y', 'z')

boldify;
set(gcf, 'PaperPositionMode', 'auto');

%% plotting position
f1 = figure(1);
subplot(4,1,1)
plot(tVecMeas, truePosesMat(:,1))
title('Position')
hold on
plot(tVec, solvedPosesMat(:,1))
plot(tVec, filteredPosesMat(:,1),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true',pose_type,'MEKF','Location','Southwest')
xlabel('time [s]')
ylabel('x [m]')
xlim([tVec(1) tVec(end)])

subplot(4,1,2)
plot(tVecMeas, truePosesMat(:,2))
hold on
plot(tVec, solvedPosesMat(:,2))
plot(tVec, filteredPosesMat(:,2),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
%legend('true',pose_type,'MEKF','Location','Northeast')
xlabel('time [s]')
ylabel('y [m]')
xlim([tVec(1) tVec(end)])

subplot(4,1,3)
plot(tVecMeas, truePosesMat(:,3))
hold on
plot(tVec, solvedPosesMat(:,3))
plot(tVec, filteredPosesMat(:,3),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
%legend('true',pose_type,'MEKF','Location','Southwest')
xlabel('time [s]')
ylabel('z [m]')
xlim([tVec(1) tVec(end)])
%ylim([0,50])

subplot(4,1,4)
plot(tVec, posScoreVec,'color',[0.8500, 0.3250, 0.0980])
hold on
plot(tVec, posScoreVecFiltered,'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend(pose_type,'MEKF')
xlabel('time [s]')
ylabel('pos\_score [m]')
xlim([tVec(1) tVec(end)])

figpos = [0 0 1920/2 1080];
set(f1,'Position',figpos)
boldify;
set(gcf, 'PaperPositionMode', 'auto');

%% plotting attitude
f2 = figure(2);
subplot(4,1,1)
plot(tVecMeas, rad2deg(truePosesMat(:,4)))
title('Attitude')
hold on
plot(tVec, rad2deg(solvedPosesMat(:,4)))
plot(tVec, rad2deg(filteredPosesMat(:,4)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true',pose_type,'MEKF','Location','Southeast')
xlabel('time [s]')
ylabel('\phi [deg]')
xlim([tVec(1) tVec(end)])
%ylim( [min(rad2deg(truePosesMat(:,4)))-5, max(rad2deg(truePosesMat(:,4)))+5] )

subplot(4,1,2)
plot(tVecMeas, rad2deg(truePosesMat(:,5)))
hold on
plot(tVec, rad2deg(solvedPosesMat(:,5)))
plot(tVec, rad2deg(filteredPosesMat(:,5)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
%legend('true',pose_type,'MEKF','Location','Southeast')
xlim([tVec(1) tVec(end)])
xlabel('time [s]')
ylabel('\theta [deg]')

subplot(4,1,3)
plot(tVecMeas, rad2deg(wrapTo2Pi(truePosesMat(:,6))))
hold on
plot(tVec, rad2deg(wrapTo2Pi(solvedPosesMat(:,6))))
plot(tVec, rad2deg(wrapTo2Pi(filteredPosesMat(:,6))),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
%legend('true',pose_type,'MEKF','Location','Southeast')
xlim([tVec(1) tVec(end)])
xlabel('time [s]')
ylabel('\psi [deg]')

subplot(4,1,4)
plot(tVec, rad2deg(attScoreVec),'color',[0.8500, 0.3250, 0.0980])
hold on
plot(tVec, rad2deg(attScoreVecFiltered),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
xlim([tVec(1) tVec(end)])
ylim([0 15])
legend(pose_type,'MEKF')
xlabel('time [s]')
ylabel('att\_score [deg]')


figpos = [0 0 1920/2 1080];
set(f2,'Position',figpos)
boldify;
set(gcf, 'PaperPositionMode', 'auto');



%%

f3 = figure(3);
for idx = 1:3,
    subplot(3,1,idx)
    plot(tVecMeas, rad2deg(trueAngVelMat(:,idx)),'-');
    hold on
    plot(tVec, rad2deg(solvedAngVelMat(:,idx)),'-');
    plot(tVec, rad2deg(filteredAngVelMat(:,idx)),'--','color',[0.4940, 0.1840, 0.5560]);%,'LineWidth',2);
    hold off
    grid on
    xlabel('time [s]')
    ylabel('deg/s')
end

% ylim([0 10])
legend('true', 'quatera', 'filtered')

boldify
set(gcf, 'PaperPositionMode', 'auto');

%%
f4 = figure(4);
subplot(2,1,1)
plot(epsHistory./(windowSize - 2))
title('eps')
grid on
boldify

subplot(2,1,2)
plot(windowSize)
title('window size')
grid on
boldify

%%
% window_start = 40;
% 
% lambda3_mean = mean(epsHistory(window_start:end))
% lambda3_std = std(epsHistory(window_start:end))
% 
% varepsSq_est = mean(epsHistory(window_start:end)./(windowSize(window_start:end)-2))
% vareps_est = mean(sqrt(epsHistory(window_start:end)./(windowSize(window_start:end)-2)))
