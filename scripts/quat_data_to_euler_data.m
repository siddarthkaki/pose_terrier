clear; close all; clc;

rot = rotation;

truePosesQuatMat = f_read_poses("../data/gt_poses.csv");
truePosesQuatMat = truePosesQuatMat(:,1:4);

[n_poses, ~] = size(truePosesQuatMat);

truePosesEulMat = zeros(n_poses,6);

for idx = 1:n_poses,
    quat = truePosesQuatMat(idx,:);
    TMat = rot.quat2Tmat(quat);
    eul = rot.dcm2euler_312(TMat);
    truePosesEulMat(idx,4:6) = eul;
end

writematrix(truePosesEulMat,'../data/true_poses.csv') 