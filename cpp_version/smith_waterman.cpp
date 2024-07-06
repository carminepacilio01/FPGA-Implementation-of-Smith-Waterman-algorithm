//  Implementation of the Smith-Waterman (Gotoh)
//  Author: Carmine Pacilio
//  Year: 2024


#include "smith_waterman.hpp"

void computeSW(short lenA, char *seqA, short lenB, char *seqB, short wd, short ws, short gap_opening, short enlargement, short *score) {
    short score_internal = 0;

    static short P_internal[3 * MAX_DIM];
    static short D_internal[3 * MAX_DIM];
    static short Q_internal[3 * MAX_DIM];

    static char seqA_internal[MAX_DIM];
    static char seqB_internal[MAX_DIM];

    const short row = lenA + 1;
    const short col = lenB + 1;
    const short matr_dim = row * col;

    const short gap_penalty = gap_opening + enlargement * 1;

    short i = 0, j = 0, itr = 0;

    // copy seqA and seqB to local mem
    memcpy(seqA_internal, seqA, lenA * sizeof(char));
    memcpy(seqB_internal, seqB, lenB * sizeof(char));
/*
 * For Backtrace
    //  Initializes first column of the matrices to help value to bound
    setupCol: for (i = 0; i < col; ++i) {
#pragma HLS PIPELINE II=1
        P_internal[i] = std::numeric_limits<short>::min() / 2;
        D_internal[i] = 0;
        Q_internal[i] = std::numeric_limits<short>::min() / 2;
    }

    //  Initializes first row of the matrices to help value to bound
    setupRow: for (j = 0; j < row; ++j) {
#pragma HLS PIPELINE II=1
        P_internal[j * col] = std::numeric_limits<short>::min() / 2;
        D_internal[j * col] = 0;
        Q_internal[j * col] = std::numeric_limits<short>::min() / 2;
    }
*/
    //  Compute P, D, Q
    compute: for (i = 1; i < row; ++i) {
        for (j = 1; j < col; ++j) {
            //  Calculate P = min{Di-1,j + g(1); Pi-1,j}
            P_internal[i * col + j] = std::max((D_internal[(i - 1) * col + j] + gap_penalty), (P_internal[(i - 1) * col + j] + enlargement));

            //  Calculate Q = min{Di,j-1 + g(1); Qi,j-1}
            Q_internal[i * col + j] = std::max((D_internal[i * col + (j - 1)] + gap_penalty) ,Q_internal[i * col + (j - 1)] + enlargement);

            //  Calculate D = min{Di-1,j-1 + w(ai, bj); Pi,j; Qij}
            //  Calculate the score based on the match or mismatch
            int score = D_internal[(i - 1) * col + (j - 1)] + ((seqA_internal[i - 1] == seqB_internal[j - 1]) ? wd : ws);

            //  Update D[i][j] based on the maximum score among score, P[i][j], and Q[i][j]
            //  Check and update D[i][j] if P[i][j] or Q[i][j] is greater than D[i][j]
            int tmp_max = std::max(P_internal[i * col + j], Q_internal[i * col + j]);
            tmp_max = std::max(score, tmp_max);

            if(D_internal[i * col + j] < tmp_max){
                D_internal[i * col + j] = tmp_max;
            } else {
                D_internal[i * col + j] = 0;
            }

            if(D_internal[i * col + j] > score_internal){
                score_internal = D_internal[i * col + j];
            }
        }
    }

    //  Write score
    memcpy(score, &score_internal, sizeof(short int));
}

//////////////////MASTER AXI
void sw_maxi(short lenA, char *seqA, short lenB, char *seqB, short wd, short ws, short gap_opening, short enlargement, short *score) {

#pragma HLS INTERFACE mode=s_axilite port=return bundle=control

#pragma HLS INTERFACE mode=m_axi port=seqA bundle=gmem depth=MAX_DIM offset=slave
#pragma HLS INTERFACE mode=m_axi port=seqB bundle=gmem depth=MAX_DIM offset=slave
#pragma HLS INTERFACE mode=m_axi port=score bundle=gmem offset=slave

#pragma HLS INTERFACE mode=s_axilite port=lenA bundle=control
#pragma HLS INTERFACE mode=s_axilite port=lenB bundle=control
#pragma HLS INTERFACE mode=s_axilite port=wd bundle=control
#pragma HLS INTERFACE mode=s_axilite port=ws bundle=control
#pragma HLS INTERFACE mode=s_axilite port=gap_opening bundle=control
#pragma HLS INTERFACE mode=s_axilite port=enlargement bundle=control
#pragma HLS INTERFACE mode=s_axilite port=seqA bundle=control
#pragma HLS INTERFACE mode=s_axilite port=seqB bundle=control
#pragma HLS INTERFACE mode=s_axilite port=score bundle=control

#pragma HLS DATAFLOW
    computeSW(lenA, seqA, lenB, seqB, wd, ws, gap_opening, enlargement, score);
}

