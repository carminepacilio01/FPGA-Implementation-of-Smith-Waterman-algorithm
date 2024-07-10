#ifndef _HOST_H
#define _HOST_H

#include <iostream>
#include <string.h>
#include <vector>
#include <limits>
#include <random>
#include <time.h>

#include "common/xcl2.hpp"

#define SEED 2

//////THE MAX_DIM ON THE FPGA BUILD IS 1024
#define MAX_SEQ_LEN 256
#define INPUT_SIZE 200

#define NUM_KERNEL 1

#endif // _HOST_H