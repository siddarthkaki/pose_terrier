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

## compute position score
posScoreVec         = 1000 * ones(num_poses,1);
posScoreVecFiltered = 1000 * ones(num_poses,1);

for idx = 1:num_poses,
    posErr         = truePosesMat[idx,1:3] - solvedPosesMat[idx,1:3];
    posErrFiltered = truePosesMat[idx,1:3] - filteredPosesMat[idx,1:3];

    posScoreVec[idx]         = norm(posErr);#/norm(truePosesMat[idx,1:3]);
    posScoreVecFiltered[idx] = norm(posErrFiltered);#/norm(truePosesMat[idx,1:3]);
end

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
gr() # Choose a backend

# shared attributes
xlims_ = (0,maximum(tVec)+dt);
labels = ["true" "NLS" "KF"];

# plotting position
datap1 = [truePosesMat[:,1] solvedPosesMat[:,1] filteredPosesMat[:,1]];
posp1 = plot(tVec, datap1, label = labels, legend = :bottomleft);
xlims!(xlims_)
ylabel!("x [m]")

datap2 = [truePosesMat[:,2] solvedPosesMat[:,2] filteredPosesMat[:,2]];
posp2 = plot(tVec, datap2, label = labels, legend = :bottomleft);
xlims!(xlims_)
ylabel!("y [m]")

datap3 = [truePosesMat[:,3] solvedPosesMat[:,3] filteredPosesMat[:,3]];
posp3 = plot(tVec, datap3, label = labels, legend = :bottomleft);
xlims!(xlims_)
ylabel!("z [m]")

datap4 = [posScoreVec posScoreVecFiltered];
posp4 = plot(tVec, datap4, label = "")
xlims!(xlims_)
xlabel!("time [s]")
ylabel!("pos_score [ ]")

plot(posp1, posp2, posp3, posp4, layout = (4,1))


#=
# plotting attitude
labels = ["true" "NLS" "KF"];
datap5 = (180.0/π)*[truePosesMat[:,4] solvedPosesMat[:,4] filteredPosesMat[:,4]];
posp5 = plot(tVec, datap5,  xlims = xlims_,
                            ylabel = "ϕ [deg]",
                            label = labels);

datap6 = (180.0/π)*[truePosesMat[:,5] solvedPosesMat[:,5] filteredPosesMat[:,5]];
posp6 = plot(tVec, datap6,  xlims = xlims_,
                            ylabel = "θ [deg]",
                            label = "");

datap7 = (180.0/π)*[truePosesMat[:,6] solvedPosesMat[:,6] filteredPosesMat[:,6]];
posp7 = plot(tVec, datap7,  xlims = xlims_,
                            ylabel = "ψ [deg]",
                            label = "");

datap8 = (180.0/π)*[attScoreVec attScoreVecFiltered];
posp8 = plot(tVec, datap8,  xlims = xlims_,
                            xlabel = "time [s]",
                            ylabel = "att_score [ ]",
                            label = "");

plot(posp5, posp6, posp7, posp8, layout = (4,1))
=#
