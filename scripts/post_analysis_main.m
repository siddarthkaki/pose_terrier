%% housekeeping
clear; close all; clc;

%% init
prefix = "../data/" + "1566312462" + "_";

tVec = f_read_timestamps(prefix + "timestamps.csv");

solvedPosesMat = f_read_poses(prefix + "solved_poses.csv");

filteredPosesMat = f_read_poses(prefix + "filtered_poses.csv");

statesMat = f_read_covars(prefix + "kf_states.csv");

covarsMat = f_read_covars(prefix + "kf_covars.csv");

[num_poses,~] = size(filteredPosesMat);

%% plotting position
figure(1)
subplot(3,1,1)
plot(tVec, solvedPosesMat(:,1))
title('Position')
hold on
plot(tVec, filteredPosesMat(:,1),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('x [m]')

subplot(3,1,2)
plot(tVec, solvedPosesMat(:,2))
hold on
plot(tVec, filteredPosesMat(:,2),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('LS','KF','Location','Southwest')
xlabel('time [s]')
ylabel('y [m]')

subplot(3,1,3)
plot(tVec, solvedPosesMat(:,3))
hold on
plot(tVec, filteredPosesMat(:,3),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('z [m]')
%ylim([0,50])


%% plotting attitude
figure(2)
subplot(3,1,1)
plot(tVec, rad2deg(solvedPosesMat(:,4)))
title('Attitude')
hold on
plot(tVec, rad2deg(filteredPosesMat(:,4)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\phi [deg]')

subplot(3,1,2)
plot(tVec, rad2deg(solvedPosesMat(:,5)))
hold on
plot(tVec, rad2deg(filteredPosesMat(:,5)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\theta [deg]')

subplot(3,1,3)
plot(tVec, rad2deg(solvedPosesMat(:,6)))
hold on
plot(tVec, rad2deg(filteredPosesMat(:,6)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\psi [deg]')