%% housekeeping
clear; close all; clc;

%% init

prefix = "../data/" + "1569774849" + "_";

tVec = f_read_timestamps(prefix + "timestamps.csv");

truePosesMat = f_read_poses("../data/true_poses.csv");

solvedPosesMat = f_read_poses(prefix + "solved_poses.csv");

filteredPosesMat = f_read_poses(prefix + "filtered_poses.csv");

%statesMat = f_read_covars(prefix + "kf_states.csv");

%covarsMat = f_read_covars(prefix + "kf_covars.csv");

[num_poses,~] = size(filteredPosesMat);

tVecMeas = f_read_timestamps("../data/meas_timestamps.csv");

%% plotting position
figure(1)
subplot(3,1,1)
plot(tVecMeas, truePosesMat(:,1))
title('Position')
hold on
plot(tVec, solvedPosesMat(:,1))
plot(tVec, filteredPosesMat(:,1),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('x [m]')

subplot(3,1,2)
plot(tVecMeas, truePosesMat(:,2))
hold on
plot(tVec, solvedPosesMat(:,2))
plot(tVec, filteredPosesMat(:,2),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','LS','KF','Location','Southwest')
xlabel('time [s]')
ylabel('y [m]')

subplot(3,1,3)
plot(tVecMeas, truePosesMat(:,3))
hold on
plot(tVec, solvedPosesMat(:,3))
plot(tVec, filteredPosesMat(:,3),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('z [m]')
%ylim([0,50])

%% plotting attitude
figure(2)
subplot(3,1,1)
plot(tVecMeas, rad2deg(truePosesMat(:,4)))
title('Attitude')
hold on
plot(tVec, rad2deg(solvedPosesMat(:,4)))
plot(tVec, rad2deg(filteredPosesMat(:,4)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\phi [deg]')

subplot(3,1,2)
plot(tVecMeas, rad2deg(truePosesMat(:,5)))
hold on
plot(tVec, rad2deg(solvedPosesMat(:,5)))
plot(tVec, rad2deg(filteredPosesMat(:,5)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\theta [deg]')

subplot(3,1,3)
plot(tVecMeas, rad2deg(truePosesMat(:,6)))
hold on
plot(tVec, rad2deg(solvedPosesMat(:,6)))
plot(tVec, rad2deg(filteredPosesMat(:,6)),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\psi [deg]')