/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <Eigen/Core>
#include <chrono>
#include <iomanip>
#include <fstream>
#include <fcntl.h>
#include <math.h>
#include <sys/stat.h>
#include "Utilities.h"
#include "measurement.pb.h"

#include "third_party/json.hpp"

using Eigen::AngleAxisd;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXd;
using nlohmann::json;

sig_atomic_t volatile finished = 0;

/**
 * @function exit_handler
 * @brief exit_handler function
 */
void exit_handler(int signal)
{
    printf("Caught signal %d\n", signal);
    finished = 1;
}

/**
 * @function main
 * @brief main function
 */
int main(int argc, char **argv)
{

    //google::InitGoogleLogging(argv[0]);

    //-- Read-in problem geometry and params ---------------------------------/

    // read params from JSON file
    std::ifstream input_stream("params.json");
    json json_params;
    try
    {
        input_stream >> json_params;
    }
    catch (json::parse_error &e)
    {
        // output exception information
        std::cout << "message: " << e.what() << '\n'
                  << "exception id: " << e.id << '\n'
                  << "byte position of error: " << e.byte << std::endl;
    }

    // specify rigid position vector of camera wrt chaser in chaser frame
    Vector3d rCamVec;
    for (unsigned int idx = 0; idx < 2; idx++)
    {
        rCamVec(idx) = json_params["rCamVec"].at(idx);
    }

    // specify camera focal length
    double focal_length = json_params["focal_length"]; //5.5*pow(10,-3);

    // specify measurement noise standard deviation (rad)
    double meas_std = double(json_params["meas_std_deg"]) * Utilities::DEG2RAD;

    // specify rigid position vector of feature points wrt target in target frame
    unsigned int num_features = json_params["rFeaMat"].size();
    MatrixXd rFeaMat(num_features, 3);
    for (unsigned int idx = 0; idx < num_features; idx++)
    {
        for (unsigned int jdx = 0; jdx < 3; jdx++)
        {
            rFeaMat(idx, jdx) = json_params["rFeaMat"][idx]["fea" + std::to_string(idx + 1)][jdx];
        }
    }

    const unsigned int num_poses_test = json_params["num_poses_test"];
    const unsigned int vector_reserve_size = json_params["vector_reserve_size"];
    const bool log_to_file = json_params["log_to_file"];

    // FIFO pipe
    int fd;
    std::string pipe_path_meas = json_params["pipe_path_input"];
    const char *myfifo = pipe_path_meas.c_str();

    //-- Init ----------------------------------------------------------------/

    // declare vectors for storage
    std::vector<Pose> true_poses;
    std::vector<double> timestamps;

    // pre-allocate memory
    true_poses.reserve(vector_reserve_size);
    timestamps.reserve(vector_reserve_size);
    
    // clock object
    auto init_t = std::chrono::high_resolution_clock::now();
    auto curr_t = std::chrono::high_resolution_clock::now();
    double curr_elapsed_t = 0.0;

    // true pose
    Pose pose_true;
    pose_true.pos << 0.5, -0.25, 30.0;
    //pose_true.quat = Quaterniond::UnitRandom();
    pose_true.quat.w() = 1.0;
    pose_true.quat.vec() = Vector3d::Zero();

    // if pipe does not exist, create it
    struct stat buf;
    if (stat(myfifo, &buf) != 0)
    {
        mkfifo(myfifo, 0666);
    }

    // set-up exit handler
    struct sigaction sigIntHandler;
    sigIntHandler.sa_handler = exit_handler;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;
    sigaction(SIGINT, &sigIntHandler, NULL);

    // log path name prefixing and postfixing
    //std::string init_time_str = std::to_string(std::time(nullptr));
    std::string prefix = "../data/";
    std::string postfix = ".csv";

    unsigned int meas_count = 1;

    //-- Loop ----------------------------------------------------------------/

    while (1) //unsigned int pose_idx = 0; pose_idx < num_poses_test; pose_idx++)
    {
        //-- Simulate Measurements -------------------------------------------/

        // generate true pose values for ith run
        pose_true.pos(0) += 0.001;
        pose_true.pos(1) -= 0.001;
        pose_true.pos(2) += 0.01;
        Quaterniond quat_step = AngleAxisd(0.001, Vector3d::UnitX()) *
                                AngleAxisd(-0.001, Vector3d::UnitY()) *
                                AngleAxisd(0.001, Vector3d::UnitZ());
        pose_true.quat = pose_true.quat * quat_step;

        // express feature points in chaser frame at the specified pose
        MatrixXd rMat = Utilities::FeaPointsTargetToChaser(pose_true, rCamVec, rFeaMat);

        // generate simulated measurements
        VectorXd yVec = Utilities::SimulateMeasurements(rMat, focal_length);

        // add Gaussian noise to simulated measurements
        VectorXd yVecNoise = Utilities::AddGaussianNoiseToVector(yVec, meas_std);

        //-- Package Measurements into ProtoBuf ------------------------------/

        // deserialise from buffer array
        ProtoMeas::Measurements measurements;

        measurements.set_num_feature_points(num_features);

        for (unsigned int idx = 0; idx < num_features; idx++)
        {
            ProtoMeas::Position *feature_point = measurements.add_feature_points();
            feature_point->set_x(rFeaMat(idx, 0));
            feature_point->set_y(rFeaMat(idx, 1));
            feature_point->set_z(rFeaMat(idx, 2));

            ProtoMeas::Bearing *bearing = measurements.add_bearings();
            bearing->set_az(yVecNoise(2 * idx + 0));
            bearing->set_el(yVecNoise(2 * idx + 1));
        }

        //-- Write to Pipe ---------------------------------------------------/

        // Open FIFO for write only, without blocking
        fd = open(myfifo, O_WRONLY | O_NONBLOCK);

        // store byte size of pose object
        size_t size = measurements.ByteSize();

        // allocate sufficient buffer space
        void *buffer = malloc(size);

        // serialise to the buffer array
        measurements.SerializeToArray(buffer, size);

        // write size of pose object to pipe
        write(fd, &size, sizeof(size));

        // write serialised pose object to pipe
        write(fd, buffer, size);

        // close FIFO
        close(fd);

        // free memory
        free(buffer);

        // timing
        curr_t = std::chrono::high_resolution_clock::now();
        curr_elapsed_t = (double)std::chrono::duration_cast<std::chrono::nanoseconds>(curr_t - init_t).count();
        curr_elapsed_t *= pow(10.0, -9.0);

        // store data
        true_poses.push_back(pose_true);
        timestamps.push_back(curr_elapsed_t);

        // write to console
        std::cout << "Sent measurement: " << meas_count << std::endl;
        meas_count++;

        //-- Handling for Program Exit ---------------------------------------/
        if (finished)
        {
            if (log_to_file)
            {
                bool append_mode = false;

                // write to csv files
                Utilities::WritePosesToCSV(true_poses, prefix + "true_poses" + postfix, append_mode);
                Utilities::WriteTimestampsToFile(timestamps, prefix + "meas_timestamps" + postfix, append_mode);

                printf("Logged data to file.\n");
            }

            // close pipe
            close(fd);

            printf("Exiting....\n");
            exit(1);
        }

        // sleep for 1.0 sec
        usleep(1000000);
    }

    //-- Close-out -----------------------------------------------------------/

    return 0;
}