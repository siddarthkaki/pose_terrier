/* Copyright (c) 2019 Siddarth Kaki
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/. */

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include "pose.pb.h"

#include "third_party/json.hpp"

using nlohmann::json;

/**
 * @function main
 * @brief main function
 */
int main(int argc, char **argv)
{
    //-- Read in pipe path params --------------------------------------------/

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

    // FIFO pipe
    int fd, rd;
    std::string pipe_path_pose = json_params["pipe_path_output"];
    const char *myfifo = pipe_path_pose.c_str();

    //-- Init ----------------------------------------------------------------/

    // ProtoPose object
    ProtoPose::Pose pose;

    // if pipe does not exist, create it
    struct stat buf;
    if (stat(myfifo, &buf) != 0)
    {
        mkfifo(myfifo, 0666); 
    }

    // Open FIFO for read only, without blocking
    fd = open(myfifo, O_RDONLY | O_NONBLOCK);

    //-- Loop ----------------------------------------------------------------/

    while (1)
    {
        // read byte size of pose object from pipe
        size_t size;
        rd = read(fd, &size, sizeof(size));

        if (rd == sizeof(size)) // if successfully received size
        {
            // allocate sufficient buffer space
            void *buffer = malloc(size);

            // read serialised pose object from pipe
            rd = read(fd, buffer, size);

            // deserialise from buffer array
            pose.ParseFromArray(buffer, size);

            // free memory
            free(buffer);

            // print received pose
            std::cout << pose.ShortDebugString() << std::endl;
        }
        else
        { // if no size message found
            std::cout << "pipe empty" << std::endl;
        }

        // sleep for 1.0 sec
        usleep(1000000);
    }

    //-- Close-out -----------------------------------------------------------/

    // close FIFO
    close(fd);

    return 0;
}