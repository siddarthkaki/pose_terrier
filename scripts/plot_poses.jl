using DelimitedFiles
using ReferenceFrameRotations
using Plots

## init
truePosesMat = readdlm("data/true_poses.csv", ',', Float64, '\n');

solvedPosesMat = readdlm("data/solved_poses.csv", ',', Float64, '\n');

filteredPosesMat = readdlm("data/filtered_poses.csv", ',', Float64, '\n');

statesMat = readdlm("data/kf_states.csv", ',', Float64, '\n');

covarsMat = readdlm("data/kf_covars.csv", ',', Float64, '\n');

num_poses = size(truePosesMat)[1];
dt = 0.5;
tVec = 0:dt:(num_poses - 1) * dt;


## compute attitude score
attScoreVec         = 1000 * ones(num_poses, 1);
attScoreVecFiltered = 1000 * ones(num_poses, 1);

for idx = 1:num_poses
    # TODO migrate to julia quat object
    #q = angle_to_quat(pi/2, pi/4, pi/3, :ZYX)
    quat            = angle_to_quat(truePosesMat[idx,4],     truePosesMat[idx,5],     truePosesMat[idx,6], :ZYX);
    quatHat         = angle_to_quat(solvedPosesMat[idx,4],   solvedPosesMat[idx,5],   solvedPosesMat[idx,6], :ZYX);
    quatHatFiltered = angle_to_quat(filteredPosesMat[idx,4], filteredPosesMat[idx,5], filteredPosesMat[idx,6], :ZYX);

    dquat         = quat * conj(quatHat);
    dquatFiltered = quat * conj(quatHatFiltered);
    #dquat         = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHat)) );
    #dquatFiltered = quatmultiply( quatnormalize(quat), quatconj(quatnormalize(quatHatFiltered)) );

    attScoreVec[idx]         = 2 * acos(abs(dquat.q0)); # rad
    attScoreVecFiltered[idx] = 2 * acos(abs(dquatFiltered.q0)); # rad
end

## plotting
plotlyjs(); # Choose a backend

# plotting position
labels = ["true" "NLS" "KF"];
datap1 = [truePosesMat[:,1] solvedPosesMat[:,1] filteredPosesMat[:,1]];
posp1 = plot(tVec, datap1,  xlims = (0,maximum(tVec)+dt),
                            ylabel = "x [m]",
                            label = labels);

datap2 = [truePosesMat[:,2] solvedPosesMat[:,2] filteredPosesMat[:,2]];
posp2 = plot(tVec, datap2,  xlims = (0,maximum(tVec)+dt),
                            ylabel = "y [m]",
                            label = "");

datap3 = [truePosesMat[:,3] solvedPosesMat[:,3] filteredPosesMat[:,3]];
posp3 = plot(tVec, datap3,  xlims = (0,maximum(tVec)+dt),
                            xlabel = "time [s]",
                            ylabel = "z [m]",
                            label = "");

plot(posp1, posp2, posp3, layout = (3,1))

#=
figure(1)
subplot(3,1,1)
plot(tVec, truePosesMat(:,1))
title('Position')
hold on
plot(tVec, solvedPosesMat(:,1))
plot(tVec, filteredPosesMat(:,1),'color',[0.4940, 0.1840, 0.5560],'LineWidth',2)
grid on
legend('true','LS','KF','Location','Northwest')
xlabel('time [s]')
ylabel('x [m]')
=#
