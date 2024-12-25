#ifndef _SMITH_WATERMAN_H
#define _SMITH_WATERMAN_H

#include <iostream>
#include <vector>
#include <limits>
#include <string.h>

#define MAX_DIM 256
#define MAX_REP MAX_DIM * 2
#define MATRIX_SIZE MAX_DIM * MAX_DIM
#define INPUT_SIZE 100

#define UP 0
#define UP_LEFT -1
#define LEFT -1

#define INFTY -32768 / 2

typedef struct conf {
	int match;
	int mismatch;
	int gap_opening;
	int gap_extension;
} conf_t;

//extern "C" {
	void sw_maxi (
		int lenT[INPUT_SIZE],
		char target[INPUT_SIZE][MAX_DIM],
		int lenD[INPUT_SIZE],
		char database[INPUT_SIZE][MAX_DIM],
		int wd, int ws, int gap_opening, int enlargement,
		int score[INPUT_SIZE],
		int offset, int input_len
	);
//}

void computeSW(int lenT, char *target, int lenD, char *database, conf_t scoring, int *score);
void readInput(int lenT[MAX_DIM], int lenT_local[MAX_DIM], char target[INPUT_SIZE][MAX_DIM], char t_local[INPUT_SIZE][MAX_DIM], int lenD[MAX_DIM], int lenD_local[MAX_DIM], char database[INPUT_SIZE][MAX_DIM], char db_local[INPUT_SIZE][MAX_DIM]);
void writeOutput(int score_l[INPUT_SIZE], int score[INPUT_SIZE], int offset, int input_len);

#endif // _SMITH_WATERMAN_H