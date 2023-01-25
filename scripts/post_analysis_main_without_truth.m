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
% prefix = "../data/data/cygnus_jp_smooth2_1664229680_";
prefix = "../data/data/cygnus_jp_smooth2_med_augment_1664229841_";
% prefix = "../data/data/cygnus_jp_smooth2_heavy_augment_1664230004_";


tVec = f_read_timestamps(prefix + "timestamps.csv");

solvedPosesMat = f_read_poses(prefix + "solved_poses.csv");

filteredPosesMat = f_read_poses(prefix + "filtered_poses.csv");

filteredCovarsDiagMat = f_read_covars(prefix + "filtered_covar_diag.csv");

[num_poses,~] = size(filteredPosesMat);


%% temp fix
solvedPosesMat(:,5) = wrapTo2Pi(solvedPosesMat(:,5));
filteredPosesMat(:,5) = wrapTo2Pi(filteredPosesMat(:,5));

%% covar -> std extraction
attStdMat = sqrt(filteredCovarsDiagMat(:,1:3));
angVelStdMat = sqrt(filteredCovarsDiagMat(:,4:6));
angAccStdMat = sqrt(filteredCovarsDiagMat(:,7:9));
posStdMat = sqrt(filteredCovarsDiagMat(:,10:12));

% %% truth interpolation
% truePosesMatInterp = zeros(size(filteredPosesMat));
% for idx = 1:6,
%     v = truePosesMat(:,idx);
%     vq1 = interp1(tVecMeas,v,tVec);
%     truePosesMatInterp(:,idx) = vq1;
% end

% %% compute position score
% 
% posScoreVec         = NaN*ones(num_poses,1);
% posScoreVecFiltered = NaN*ones(num_poses,1);
% 
% for idx = 1:num_poses,
% 
%     posErr         = truePosesMatInterp(idx,1:3) - solvedPosesMat(idx,1:3);
%     posErrFiltered = truePosesMatInterp(idx,1:3) - filteredPosesMat(idx,1:3);
% 
%     posScoreVec(idx)         = norm(posErr);%/norm(truePosesMat(idx,1:3));
%     posScoreVecFiltered(idx) = norm(posErrFiltered);%/norm(truePosesMat(idx,1:3));
% end

% %% compute attitude score
% 
% attScoreVec         = 1000*ones(num_poses,1);
% attScoreVecFiltered = 1000*ones(num_poses,1);
% 
% for idx = 1:num_poses,
% 
%     quat            = angle2quat(truePosesMatInterp(idx,4), truePosesMatInterp(idx,5), truePosesMatInterp(idx,6) );
%     quatHat         = angle2quat(    solvedPosesMat(idx,4),     solvedPosesMat(idx,5),     solvedPosesMat(idx,6) );
%     quatHatFiltered = angle2quat(  filteredPosesMat(idx,4),   filteredPosesMat(idx,5),   filteredPosesMat(idx,6) );
%     
%     dquat         = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHat)) );
%     dquatFiltered = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHatFiltered)) );
%     %dqVec = (quat.normalized())*(quatHat.normalized().conjugate());
% 
%     attScoreVec(idx)         = 2*acos( abs( dquat(1) ) ); % rad
%     attScoreVecFiltered(idx) = 2*acos( abs( dquatFiltered(1) ) ); % rad
% end

%% TEST plotting 3sigmas
f100 = figure(100);
subplot(2,1,1)
plot(tVec, 3*rad2deg(attStdMat));
grid on
ylabel('deg')
ylim([0 10])
title('3\sigma')
legend('\phi', '\theta', '\psi')

subplot(2,1,2)
plot(tVec, 3*posStdMat);
grid on
ylabel('m')
legend('x', 'y', 'z')

boldify;
set(gcf, 'PaperPositionMode', 'auto');


%% att 3sigmas
f101 = figure(101);
subplot(3,1,1)
plot(tVec, 3*rad2deg(attStdMat));
grid on
ylabel('deg')
ylim([0 15])
title('3\sigma')

subplot(3,1,2)
plot(tVec, 3*rad2deg(angVelStdMat));
grid on
ylabel('deg/s')
% ylim([0 10])

subplot(3,1,3)
plot(tVec, 3*rad2deg(angAccStdMat));
grid on
ylabel('deg/s^2')
% ylim([0 10])

legend('\phi', '\theta', '\psi')

boldify;
set(gcf, 'PaperPositionMode', 'auto');

%% plotting position
f1 = figure(1);
subplot(3,1,1)
title('Position')
hold on
plot(tVec, solvedPosesMat(:,1))
plot(tVec, filteredPosesMat(:,1),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend(pose_type,'MEKF','Location','Southwest')
xlabel('time [s]')
ylabel('x [m]')
xlim([tVec(1) tVec(end)])

subplot(3,1,2)
hold on
plot(tVec, solvedPosesMat(:,2))
plot(tVec, filteredPosesMat(:,2),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
%legend('true',pose_type,'MEKF','Location','Northeast')
xlabel('time [s]')
ylabel('y [m]')
xlim([tVec(1) tVec(end)])

subplot(3,1,3)
hold on
plot(tVec, solvedPosesMat(:,3))
plot(tVec, filteredPosesMat(:,3),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
%legend('true',pose_type,'MEKF','Location','Southwest')
xlabel('time [s]')
ylabel('z [m]')
xlim([tVec(1) tVec(end)])
%ylim([0,50])

boldify;
set(gcf, 'PaperPositionMode', 'auto');


%% plotting attitude
f2 = figure(2);
subplot(3,1,1)
title('Attitude')
hold on
plot(tVec, rad2deg(solvedPosesMat(:,4)))
plot(tVec, rad2deg(filteredPosesMat(:,4)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend(pose_type,'MEKF','Location','Southeast')
xlabel('time [s]')
ylabel('\phi [deg]')
xlim([tVec(1) tVec(end)])
%ylim( [min(rad2deg(truePosesMat(:,4)))-5, max(rad2deg(truePosesMat(:,4)))+5] )

subplot(3,1,2)
hold on
plot(tVec, rad2deg(solvedPosesMat(:,5)))
plot(tVec, rad2deg(filteredPosesMat(:,5)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
%legend('true',pose_type,'MEKF','Location','Southeast')
xlim([tVec(1) tVec(end)])
xlabel('time [s]')
ylabel('\theta [deg]')

subplot(3,1,3)
hold on
plot(tVec, rad2deg(wrapTo2Pi(solvedPosesMat(:,6))))
plot(tVec, rad2deg(wrapTo2Pi(filteredPosesMat(:,6))),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
%legend('true',pose_type,'MEKF','Location','Southeast')
xlim([tVec(1) tVec(end)])
xlabel('time [s]')
ylabel('\psi [deg]')

figpos = [0 0 1920/2 1080];
set(f2,'Position',figpos)
boldify;
set(gcf, 'PaperPositionMode', 'auto');