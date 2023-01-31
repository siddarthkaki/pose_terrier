function [data] = f_read_truth_logs(prefix)

%% Import data from text file

tVec = f_read_timestamps(prefix + "meas_timestamps.csv");

TruePosesMat = f_read_poses(prefix + "true_poses.csv");
TruePosMat = TruePosesMat(:,1:3);
TrueEulMat = TruePosesMat(:,4:6);

TrueAngVelMat = f_read_omegas(prefix + "true_omegas.csv");


n_elem = length(tVec);

%%
TrueRotation = NaN(n_elem,3);
TrueQuat = NaN(n_elem,4);

for idx = 1:n_elem,
    eulvec = TrueEulMat(idx,:);

    temp_Tmat     = rotation.angleaxis2Rmat(eulvec(1),[1 0 0]) ...
                  * rotation.angleaxis2Rmat(eulvec(2),[0 1 0]) ...
                  * rotation.angleaxis2Rmat(eulvec(3),[0 0 1]);

    TrueRotation(idx,:) = rotation.dcm2euler_312(temp_Tmat);
    TrueQuat(idx,:) = dcm2quat(temp_Tmat);
end

%%
data.TrueRotation = rad2deg(TrueRotation);
data.TrueQuat = TrueQuat;
data.TrueAngVel = rad2deg(TrueAngVelMat);
data.TrueTranslation = TruePosMat;

data.Time = tVec;
