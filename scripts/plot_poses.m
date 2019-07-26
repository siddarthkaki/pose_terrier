%% housekeeping
clear; close all; clc;

%% init

truePosesMat = f_read_poses("../data/true_poses.csv");

solvedPosesMat = f_read_poses("../data/solved_poses.csv");

filteredPosesMat = f_read_poses("../data/filtered_poses.csv");

[num_poses,~] = size(truePosesMat);
dt = 0.0005;
tVec = 0:dt:(num_poses-1)*dt;

%% plotting position
figure(1)
subplot(3,1,1)
plot(tVec, truePosesMat(:,1))
hold on
plot(tVec, solvedPosesMat(:,1))
plot(tVec, filteredPosesMat(:,1))
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('x [m]')

subplot(3,1,2)
plot(tVec, truePosesMat(:,2))
hold on
plot(tVec, solvedPosesMat(:,2))
plot(tVec, filteredPosesMat(:,2))
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('y [m]')

subplot(3,1,3)
plot(tVec, truePosesMat(:,3))
hold on
plot(tVec, solvedPosesMat(:,3))
plot(tVec, filteredPosesMat(:,3))
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('z [m]')

%% plotting attitude
figure(2)
subplot(3,1,1)
plot(tVec, truePosesMat(:,4))
hold on
plot(tVec, solvedPosesMat(:,4))
plot(tVec, filteredPosesMat(:,4))
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\phi [rad]')

subplot(3,1,2)
plot(tVec, truePosesMat(:,5))
hold on
plot(tVec, solvedPosesMat(:,5))
plot(tVec, filteredPosesMat(:,5))
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\theta [rad]')

subplot(3,1,3)
plot(tVec, truePosesMat(:,6))
hold on
plot(tVec, solvedPosesMat(:,6))
plot(tVec, filteredPosesMat(:,6))
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('\psi [rad]')