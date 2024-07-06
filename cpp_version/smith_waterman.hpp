#ifndef _SMITH_WATERMAN_HPP
#define _SMITH_WATERMAN_HPP

#include <iostream>
#include <string.h>
#include <vector>
#include <limits>
#include <hls_stream.h>
#include <ap_int.h>

#define MAX_DIM 256
#define MATRIX_SIZE MAX_DIM * MAX_DIM

void computeSW(short lenA, char *seqA, short lenB, char *seqB, short wd, short ws, short gap_opening, short enlargement, short *score);
void sw_maxi(short lenA, char *seqA, short lenB, char *seqB, short wd, short ws, short gap_opening, short enlargement, short *score);

#endif // _SMITH_WATERMAN_HPP
