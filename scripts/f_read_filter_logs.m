function [data] = f_read_filter_logs(prefix)

%% Import data from text file

tVec = f_read_timestamps(prefix + "timestamps.csv");

% truePosesMat = f_read_poses(folder + "true_poses.csv");

solvedPosesMat = f_read_poses(prefix + "solved_poses.csv");
solvedPosMat = solvedPosesMat(:,1:3);
solvedEulMat = solvedPosesMat(:,4:6);

filteredPosesMat = f_read_poses(prefix + "filtered_poses.csv");
filteredPosMat = filteredPosesMat(:,1:3);
filteredEulMat = filteredPosesMat(:,4:6);

% trueAngVelMat = f_read_omegas(folder + "true_omegas.csv");

filteredAngVelMat = f_read_omegas(prefix + "filtered_omegas.csv");

filteredAngAccMat = f_read_omegas(prefix + "filtered_alphas.csv");

filteredPosStates = f_read_pos_states(prefix + "filtered_pos_states.csv");
%filteredPosMat = filteredPosStates(:,1:3);
filteredVelMat = filteredPosStates(:,4:6);
filteredAccMat = filteredPosStates(:,7:9);

filteredCovarsDiagMat = f_read_covars(prefix + "filtered_covar_diag.csv");

n_elem = length(tVec);

%%
SolvedRotation = NaN(n_elem,3);
FilteredRotation = NaN(n_elem,3);
FilteredQuat = NaN(n_elem,4);

for idx = 1:n_elem,
    eulvec = solvedEulMat(idx,:);

    temp_Tmat     = rotation.angleaxis2Rmat(eulvec(1),[1 0 0]) ...
                  * rotation.angleaxis2Rmat(eulvec(2),[0 1 0]) ...
                  * rotation.angleaxis2Rmat(eulvec(3),[0 0 1]);

    SolvedRotation(idx,:) = rotation.dcm2euler_312(temp_Tmat);

    %%%

    eulvec = filteredEulMat(idx,:);

    temp_Tmat     = rotation.angleaxis2Rmat(eulvec(1),[1 0 0]) ...
                  * rotation.angleaxis2Rmat(eulvec(2),[0 1 0]) ...
                  * rotation.angleaxis2Rmat(eulvec(3),[0 0 1]);

    FilteredRotation(idx,:) = rotation.dcm2euler_312(temp_Tmat);

    FilteredQuat(idx,:) = dcm2quat(temp_Tmat);
end

%%
data.SolvedRotation = rad2deg(SolvedRotation);
data.FilteredRotation = rad2deg(FilteredRotation);
data.FilteredQuat = FilteredQuat;
data.FilteredAngVel = rad2deg(filteredAngVelMat);
data.FilteredAngAcc = rad2deg(filteredAngAccMat);

data.SolvedTranslation = solvedPosMat;
data.FilteredTranslation = filteredPosMat;
data.FilteredVelocity = filteredVelMat;
data.FilteredAcceleration = filteredAccMat;
data.FilteredStd = sqrt(filteredCovarsDiagMat);
% data.FilteredStd(:,1:9) = rad2deg(data.FilteredStd(:,1:9));
data.Time = tVec;
