%% housekeeping
clear; close all; clc;

%% init

pose_type = 'NLS'; % 'PnP'

% folder = "long_test/";
% folder = "short_test/";
% folder = "new_test/";
folder = "../data/";%500s_noadapt_quatera/";

% prefix = "../data/long_test/" + "1661875312" + "_";
% prefix = "../data/short_test/" + "1661872620" + "_";
% prefix = "../data/new_test/" + "1661975852" + "_";
prefix = folder + "1672513095" + "_";

ptl_data = f_read_filter_logs(prefix);

tVec = f_read_timestamps(prefix + "timestamps.csv");

truePosesMat = f_read_poses(folder + "true_poses.csv");

solvedPosesMat = f_read_poses(prefix + "solved_poses.csv");

filteredPosesMat = f_read_poses(prefix + "filtered_poses.csv");

trueAngVelMat = f_read_omegas(folder + "true_omegas.csv");

solvedAngVelMat = f_read_omegas(prefix + "solved_omegas.csv");

filteredAngVelMat = f_read_omegas(prefix + "filtered_omegas.csv");

filteredAngAccMat = f_read_omegas(prefix + "filtered_alphas.csv");

filteredPosStates = f_read_pos_states(prefix + "filtered_pos_states.csv");
filteredVelMat = filteredPosStates(:,4:6);
filteredAccMat = filteredPosStates(:,7:9);

filteredCovarsDiagMat = f_read_covars(prefix + "filtered_covar_diag.csv");

% statesMat = f_read_covars(prefix + "kf_states.csv");

% covarsMat = f_read_covars(prefix + "kf_covars.csv");

[num_poses,~] = size(filteredPosesMat);

tVecMeas = f_read_timestamps(folder + "meas_timestamps.csv");

windowSize = f_read_timestamps(prefix + "window_size.csv");

epsHistory = f_read_timestamps(prefix + "eps_history.csv");


%% temp fix
% [num_true_poses,~] = size(truePosesMat);
% truePosesMat = [truePosesMat; zeros(num_poses - num_true_poses,6)];
% tVecMeas = 0:0.04:num_poses-1;

for ddx = 4:6,
    truePosesMat(:,ddx) = wrapToPi(truePosesMat(:,ddx));
    solvedPosesMat(:,ddx) = wrapToPi(solvedPosesMat(:,ddx));
    filteredPosesMat(:,ddx) = wrapToPi(filteredPosesMat(:,ddx));
end

%% covar -> std extraction
attStdMat = ptl_data.FilteredStd(:,1:3);
angVelStdMat = ptl_data.FilteredStd(:,4:6);
angAccStdMat = ptl_data.FilteredStd(:,7:9);
posStdMat = ptl_data.FilteredStd(:,10:12);

%% truth interpolation
truePosesMatInterp = zeros(size(filteredPosesMat));
for idx = 1:6,
    v = truePosesMat(:,idx);
    vq1 = interp1(tVecMeas,v, ptl_data.Time);
    truePosesMatInterp(:,idx) = vq1;
    truePosesMatInterp(:,4:6) = wrapToPi(truePosesMatInterp(:,4:6));
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

trueQuatMat = NaN(num_poses,4);
solvedQuatMat = NaN(num_poses,4);

c_err_est_mat = NaN(num_poses,4);

for idx = 1:num_poses,

    temp_Tmat     = rotation.angleaxis2Rmat(truePosesMatInterp(idx,4),[1 0 0]) ...
                  * rotation.angleaxis2Rmat(truePosesMatInterp(idx,5),[0 1 0]) ...
                  * rotation.angleaxis2Rmat(truePosesMatInterp(idx,6),[0 0 1]);
    quat = dcm2quat(temp_Tmat);
    trueQuatMat(idx,:) = quat;

    temp_Tmat     = rotation.angleaxis2Rmat(solvedPosesMat(idx,4),[1 0 0]) ...
                  * rotation.angleaxis2Rmat(solvedPosesMat(idx,5),[0 1 0]) ...
                  * rotation.angleaxis2Rmat(solvedPosesMat(idx,6),[0 0 1]);
    quatHat = dcm2quat(temp_Tmat);
    solvedQuatMat(idx,:) = quatHat;

    temp_Tmat     = rotation.angleaxis2Rmat(filteredPosesMat(idx,4),[1 0 0]) ...
                  * rotation.angleaxis2Rmat(filteredPosesMat(idx,5),[0 1 0]) ...
                  * rotation.angleaxis2Rmat(filteredPosesMat(idx,6),[0 0 1]);
    quatHatFiltered = dcm2quat(temp_Tmat);

    %quat            = angle2quat(truePosesMatInterp(idx,4), truePosesMatInterp(idx,5), truePosesMatInterp(idx,6) );
    %quatHat         = angle2quat(    solvedPosesMat(idx,4),     solvedPosesMat(idx,5),     solvedPosesMat(idx,6) );
    %quatHatFiltered = angle2quat(  filteredPosesMat(idx,4),   filteredPosesMat(idx,5),   filteredPosesMat(idx,6) );
    
    dquat         = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHat)) );
    dquatFiltered = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHatFiltered)) );
    %dqVec = (quat.normalized())*(quatHat.normalized().conjugate());

    attScoreVec(idx)         = 2*acos( abs( dquat(1) ) ); % rad
    attScoreVecFiltered(idx) = 2*acos( abs( dquatFiltered(1) ) ); % rad


    %%%
    a = quat.';
    c = quatHat;
    b = rotation.quatmult_S(c, quatinv(a.').');

    c_err_est(1,1) = -b(2:4).'*a(2:4);
    c_err_est(2,1) = a(1)*b(2) + a(3)*b(4) - a(4)*b(3);
    c_err_est(3,1) = a(1)*b(3) - a(2)*b(4) + a(4)*b(2);
    c_err_est(4,1) = a(1)*b(4) + a(2)*b(3) - a(3)*b(2);

    c_err_est_mat(idx,:) = c_err_est;
end

f1111 = figure(1111);
for idx = 1:4,
    c_err_est_mat_denan(:,idx) = denan(c_err_est_mat(:,idx));
    r = c_err_est_mat(:,idx);
    pd = fitdist(r,'Normal');
    subplot(4,1,idx)
    histfit(r,100)
    grid on
    ylabel("\Delta q_{\eta}(" + string(idx-1) + ")")
    set(gca,'FontWeight','bold')

    %if idx == 1, xlim(1.8*[-1 1]*10^-3); end
    %if idx ~= 1, xlim(0.08*[-1 1]); end
end
fontsize(f1111,scale=2)
set(gcf, 'PaperPositionMode', 'auto');

dqm = mean(c_err_est_mat_denan).'
dqs = std(c_err_est_mat_denan).'

%%
ufm = rad2deg(mean(denan(attScoreVec)))
flm = rad2deg(mean(denan(attScoreVecFiltered)))
ufs = rad2deg(std(denan(attScoreVec)))
fls = rad2deg(std(denan(attScoreVecFiltered)))
% epsmean = meanr(epsHistory(100:end))
% epsstd = std(epsHistory(100:end))

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


%% att 3sigmas
f101 = figure(101);
subplot(3,1,1)
plot(tVec, 3*rad2deg(attStdMat));
grid on
ylabel('deg')
ylim([0 10])
title('3\sigma')

subplot(3,1,2)
plot(tVec, 3*rad2deg(angVelStdMat));
grid on
ylabel('deg/s')
ylim([0 10])

subplot(3,1,3)
plot(tVec, 3*rad2deg(angAccStdMat));
grid on
ylabel('deg/s^2')
ylim([0 10])

legend('\phi', '\theta', '\psi')

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
legend('true',pose_type,'MEKF','Location','Northeast')
xlabel('time [s]')
ylabel('y [m]')
xlim([tVec(1) tVec(end)])

subplot(4,1,3)
plot(tVecMeas, truePosesMat(:,3))
hold on
plot(tVec, solvedPosesMat(:,3))
plot(tVec, filteredPosesMat(:,3),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true',pose_type,'MEKF','Location','Southwest')
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
ylim([0 5])
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
legend('true',pose_type,'MEKF','Location','Southeast')
xlim([tVec(1) tVec(end)])
xlabel('time [s]')
ylabel('\theta [deg]')

subplot(4,1,3)
plot(tVecMeas, rad2deg(wrapTo2Pi(truePosesMat(:,6))))
hold on
plot(tVec, rad2deg(wrapTo2Pi(solvedPosesMat(:,6))))
plot(tVec, rad2deg(wrapTo2Pi(filteredPosesMat(:,6))),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true',pose_type,'MEKF','Location','Southeast')
xlim([tVec(1) tVec(end)])
xlabel('time [s]')
ylabel('\psi [deg]')

subplot(4,1,4)
plot(tVec, rad2deg(attScoreVec),'color',[0.8500, 0.3250, 0.0980])
hold on
plot(tVec, rad2deg(attScoreVecFiltered),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
xlim([tVec(1) tVec(end)])
ylim([0 40])
legend(pose_type,'MEKF')
xlabel('time [s]')
ylabel('att\_score [deg]')


figpos = [0 0 1920/2 1080];
set(f2,'Position',figpos)
boldify;
set(gcf, 'PaperPositionMode', 'auto');



%% plotting angular velocity

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
    if idx == 1, title('Angular Velocity'); end
end

% ylim([0 10])
legend('true', 'quatera', 'filtered')

boldify
set(gcf, 'PaperPositionMode', 'auto');

%% plotting angular acceleration

f13 = figure(13);
for idx = 1:3,
    subplot(3,1,idx)
    plot(tVec, rad2deg(filteredAngAccMat(:,idx)),'--','color',[0.4940, 0.1840, 0.5560]);%,'LineWidth',2);
    hold off
    grid on
    xlabel('time [s]')
    ylabel('deg/s^2')
    if idx == 1, title('Angular Acceleration'); end
end

% ylim([0 10])
legend('filtered')

boldify
set(gcf, 'PaperPositionMode', 'auto');

%% plotting translational velocity

f4 = figure(4);
for idx = 1:3,
    subplot(3,1,idx)
    plot(tVec, filteredVelMat(:,idx),'--','color',[0.4940, 0.1840, 0.5560]);%,'LineWidth',2);
    hold off
    grid on
    xlabel('time [s]')
    ylabel('m')
    if idx == 1, title('Translational Velocity'); end
end

% ylim([0 10])
legend('filtered')

boldify
set(gcf, 'PaperPositionMode', 'auto');

%% plotting translational acceleration

f14 = figure(14);
for idx = 1:3,
    subplot(3,1,idx)
    plot(tVec, filteredAccMat(:,idx),'--','color',[0.4940, 0.1840, 0.5560]);%,'LineWidth',2);
    hold off
    grid on
    xlabel('time [s]')
    ylabel('m/s^2')
    if idx == 1, title('Translational Acceleration'); end
end

% ylim([0 10])
legend('filtered')

boldify
set(gcf, 'PaperPositionMode', 'auto');

%%
f5 = figure(5);
subplot(2,1,1)
plot(epsHistory./(windowSize - 2)*2)
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
