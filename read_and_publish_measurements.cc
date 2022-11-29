/* Copyright (c) 2022 Siddarth Kaki
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

#include "glog/logging.h"
#include "Utilities.h"
#include "pose.pb.h"
#include "measurement.pb.h"

#include "third_party/CppRot/cpprot.h"
#include "third_party/json.hpp"

using Eigen::AngleAxisd;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXd;
using nlohmann::json;

typedef std::vector< double > state_type;

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

    std::srand((unsigned int) time(NULL));

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

    double dt = json_params["meas_sim_dt"]; // sec

    const unsigned int num_poses_test = json_params["num_poses_test"];
    const unsigned int vector_reserve_size = json_params["vector_reserve_size"];

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
    
    // FIFO pipe
    int fd;
    std::string pipe_path_meas = json_params["pipe_path_input"];
    const char *myfifo = pipe_path_meas.c_str();

    //-- Init ----------------------------------------------------------------/

    // declare vectors for storage
    std::vector<VectorXd> noisy_measurements;
    std::vector<double> timestamps;

    // pre-allocate memory
    noisy_measurements.reserve(vector_reserve_size);
    timestamps.reserve(vector_reserve_size);
    
    // clock object
    auto init_t = std::chrono::high_resolution_clock::now();
    auto curr_t = std::chrono::high_resolution_clock::now();
    double curr_elapsed_t = 0.0;

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

    // specify input file
    std::ifstream input_stream_csv(prefix + "noisy_measurements" + postfix);
    
    // Make sure the file is open
    if(!input_stream_csv.is_open()) throw std::runtime_error("Could not open file");

    // Helper vars
    std::string line, colname;
    double val;

    //-- Loop ----------------------------------------------------------------/

    // Read data, line by line
    while(std::getline(input_stream_csv, line))
    {
        // Create a stringstream of the current line
        std::stringstream ss(line);
       
        // Keep track of the current column index
        unsigned int colIdx = 0;

        VectorXd yVecNoise(num_features*2);
            
        // Extract each double
        while(ss >> val)
        {
            yVecNoise(colIdx) = val;
            
            // If the next token is a comma, ignore it and move on
            if(ss.peek() == ',') ss.ignore();
                
            // Increment the column index
            colIdx++;
        }

        //std::cout << yVecNoise << std::endl;

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

        // write to console
        std::cout << "Sent measurement: " << meas_count << std::endl;
        meas_count++;

        //-- Handling for Program Exit ---------------------------------------/
        if (finished)
        {
            // close pipe
            close(fd);

            printf("Exiting....\n");
            exit(1);
        }

        // sleep for dt sec
        usleep(dt*1000000);
    }

    input_stream_csv.close();

    //-- Close-out -----------------------------------------------------------/

    return 0;
}