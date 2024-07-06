// Implementation of the Smith-Waterman (Gotoh)
//
//  Author: Carmine Pacilio
//	Year: 2024
//	Version: 1.0
////////////////////////////////////////////////////
#include <limits.h>
#include <stdint.h>
#include <stdio.h>

#include "smith_waterman.h"

//  Tripcount identifier
const unsigned int out_dim = MAX_DIM;

void smithWaterman(int lenA, char *seqA, int lenB, char *seqB,
		short int *P_out, short int *D_out, short int *Q_out,
		int wd, int ws,int gap_opening, int enlargement) {

	short int P[MAX_DIM * MAX_DIM];
	short int D[MAX_DIM * MAX_DIM];
	short int Q[MAX_DIM * MAX_DIM];

    int i = 0, j = 0, itr = 0;

    const short int gap_penalty = gap_opening + enlargement * 1;

    //	Initializes first column of the matrices to help value to bound
    for (i = 0; i < lenA + 1; ++i) {
		P[i * MAX_DIM] = INT16_MIN / 2;
		D[i * MAX_DIM] = 0;
		Q[i * MAX_DIM] = INT16_MIN / 2;
	}

    //	Initializes first row of the matrices to help value to bound
    for (j = 0; j < lenB + 1; ++j) {
		P[j] = INT16_MIN / 2;
		D[j] = 0;
		Q[j] = INT16_MIN / 2;
	}

    //	Filling matrix phase
    for (i = 1; i < lenA + 1; ++i) {
        for (j = 1; j < lenB + 1; ++j) {
            //	Calculate P = min{Di-1,j + g(1); Pi-1,j}
            P[i * MAX_DIM + j] = (D[(i - 1) * MAX_DIM + j] + gap_penalty > P[(i - 1) * MAX_DIM + j] + enlargement) ?
                             D[(i - 1) * MAX_DIM + j] + gap_penalty : P[(i - 1) * MAX_DIM + j] + enlargement;

            //	Calculate Q = min{Di,j-1 + g(1); Qi,j-1}
            Q[i * MAX_DIM + j] = (D[i * MAX_DIM + (j - 1)] + gap_penalty > Q[i * MAX_DIM + (j - 1)] + enlargement) ?
                             D[i * MAX_DIM + (j - 1)] + gap_penalty : Q[i * MAX_DIM + (j - 1)] + enlargement;

            //	Calculate D = min{Di-1,j-1 + w(ai, bj); Pi,j; Qij}
            //	Calculate the score based on the match or mismatch
            int score = D[(i - 1) * MAX_DIM + (j - 1)] + ((seqA[i - 1] == seqB[j - 1]) ? wd : ws);

            //	Update D[i][j] based on the maximum score among score, P[i][j], and Q[i][j]
            if (score > 0) {
                D[i * MAX_DIM + j] = score;
            } else {
                D[i * MAX_DIM + j] = 0;
            }

            //	Check and update D[i][j] if P[i][j] or Q[i][j] is greater than D[i][j]
            if (P[i * MAX_DIM + j] > D[i * MAX_DIM + j]) {
                D[i * MAX_DIM + j] = P[i * MAX_DIM + j];
            }

            if (Q[i * MAX_DIM + j] > D[i * MAX_DIM + j]) {
                D[i * MAX_DIM + j] = Q[i * MAX_DIM + j];
            }
        }
    }

    // Burst write from output matrices to global memory
    write: for (itr = 0, i = 0, j = 0; itr < (lenA + 1) * (lenB + 1); itr++, j++) {

        if (j == lenB + 1) {
            j = 0;
            i++;
        }

        P_out[itr] = P[i * MAX_DIM + j];
        D_out[itr] = D[i * MAX_DIM + j];
        Q_out[itr] = Q[i * MAX_DIM + j];
    }
}
