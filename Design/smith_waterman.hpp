#ifndef _SMITH_WATERMAN_HPP
#define _SMITH_WATERMAN_HPP

#include <iostream>
#include <string.h>
#include <vector>
#include <limits>



#define MAX_DIM 1024
#define MAX_REP MAX_DIM * 2
#define MATRIX_SIZE MAX_DIM * MAX_DIM
#define INPUT_SIZE 20000

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

void computeSW(int lenT, char *target, int lenD, char *database, conf_t scoring, int *score);
void sw_maxi(int lenT[INPUT_SIZE], char target[INPUT_SIZE][MAX_DIM], int lenD[INPUT_SIZE], char database[INPUT_SIZE][MAX_DIM], int wd, int ws, int gap_opening, int enlargement, int score[INPUT_SIZE], int input_len);
void compute_diag(int num_diag, int lenT, int &diag_len, int max_diag_len, int &database_cursor, int &target_cursor, int &n_diag_repeat, int &diag_index, char target_l[MAX_DIM], char database_l[MAX_DIM], conf_t scoring, int p_buffer[3][MAX_DIM * 2], int q_buffer[3][MAX_DIM * 2], int d_buffer[3][MAX_DIM * 2], int &score_l);
void readInput(int lenT[MAX_DIM], int lenT_local[MAX_DIM], char target[INPUT_SIZE][MAX_DIM], char t_local[INPUT_SIZE][MAX_DIM], int lenD[MAX_DIM], int lenD_local[MAX_DIM], char database[INPUT_SIZE][MAX_DIM], char db_local[INPUT_SIZE][MAX_DIM]);
void writeOutput(int score_l[INPUT_SIZE], int score[INPUT_SIZE]);

#endif // _SMITH_WATERMAN_HPP
