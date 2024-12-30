//	Implementation of the Smith-Waterman (Gotoh)
//  Author: Carmine Pacilio
//	Year: 2024
#include "smith_waterman.h"

void computeSW(int lenT, char *target, int lenD, char *database, conf_t scoring, int *score) {

//	initialize local buffer for the 3 dependency
	short p_buffer[3][MAX_DIM];
#pragma HLS ARRAY_PARTITION variable=p_buffer type=block factor=3 dim=1
#pragma HLS ARRAY_PARTITION variable=p_buffer type=cyclic factor=2 dim=2

	short q_buffer[3][MAX_DIM];
#pragma HLS ARRAY_PARTITION variable=q_buffer type=block factor=3 dim=1
#pragma HLS ARRAY_PARTITION variable=q_buffer type=cyclic factor=2 dim=2

	short d_buffer[3][MAX_DIM];
#pragma HLS ARRAY_PARTITION variable=d_buffer type=block factor=3 dim=1
#pragma HLS ARRAY_PARTITION variable=d_buffer type=cyclic factor=2 dim=2

    init_buffers_loop: for(int index = 0; index < MAX_DIM; index++) {
#pragma HLS PIPELINE  II=1
		p_buffer[0][index] = 0;
		q_buffer[0][index] = 0;
		d_buffer[0][index] = 0;

		p_buffer[1][index] = 0;
		q_buffer[1][index] = 0;
		d_buffer[1][index] = 0;

		p_buffer[2][index] = 0;
		q_buffer[2][index] = 0;
		d_buffer[2][index] = 0;
    }

    // Computation of the score by calling PEs
    const short max_diag_len = (lenT < lenD) ? lenT : lenD; // the maximum value diag_len can reach
    const short t_diag = lenT + lenD - 1;
    short diag_len = 0; // the number of elements on the current diagonal
    short n_diag_repeat = (lenT > lenD) ? (lenT - lenD) + 1 : (lenD - lenT) + 1;
    short max_score = 0;
    short score_l[MAX_DIM];

    int diag_index = 0;
    int database_cursor = 0;
    int target_cursor = 0;
	int score_index = 0;

    compute_score_loop: for (int num_diag = 0; num_diag <  MAX_DIM * 2 - 1; num_diag++) {
    	if(num_diag < t_diag){
			if (num_diag < lenT) {
				if (diag_len < max_diag_len) {
					diag_len++;
					database_cursor = 0;
					target_cursor = num_diag;
				} else if (n_diag_repeat > 0) {
					database_cursor = 0;
					target_cursor = num_diag;
					n_diag_repeat--;
					if (n_diag_repeat == 0) {
						diag_len--;
					}
				}
			} else if (n_diag_repeat > 0) {
				database_cursor++;
				target_cursor = lenT + database_cursor - 1;
				n_diag_repeat--;
				if (n_diag_repeat == 0) {
					diag_len--;
				}
			} else {
				database_cursor++;
				target_cursor = lenT + database_cursor - 1;
				diag_len--;
			}

			// PE
			compute_diagonal_loop: for (int j = 0; j < MAX_DIM; ++j) {
#pragma HLS DEPENDENCE variable=p_buffer array inter false
#pragma HLS DEPENDENCE variable=q_buffer array inter false
#pragma HLS DEPENDENCE variable=d_buffer array inter false
#pragma HLS DEPENDENCE variable=score_l inter false
#pragma HLS DEPENDENCE variable=score_l intra false
#pragma HLS UNROLL factor=2
				if(j < diag_len) {
					const int current_index = database_cursor + j;
					const int two_diag = (diag_index == 0) ? 1 : (diag_index == 1) ? 2 : 0;
					const int one_diag = (diag_index == 0) ? 2 : diag_index - 1;

					const short t = target[lenT - 1 - (target_cursor - current_index)];
					const short d = database[current_index];
					short match = (t == d) ? scoring.match : scoring.mismatch;

					const short tmp_dUP = d_buffer[one_diag][current_index + UP] + scoring.gap_opening;
					const short tmp_pUP = p_buffer[one_diag][current_index + UP] + scoring.gap_extension;
					const short tmp_dLEFT = d_buffer[one_diag][current_index + LEFT] + scoring.gap_opening;
					const short tmp_dUPLEFT = d_buffer[two_diag][current_index + UP_LEFT];
					const short tmp_qLEFT = q_buffer[one_diag][current_index + LEFT] + scoring.gap_extension;
					short tmp_p, tmp_q, tmp_d;

					if (t == '-') {
						tmp_d = 0;
						tmp_p = INFTY;
						tmp_q = INFTY;
						score_l[j] = 0;
					} else if (d == '-') {
						tmp_d = 0;
						tmp_p = INFTY;
						tmp_q = INFTY;
						score_l[j] = 0;
					} else {
						tmp_p = (tmp_dUP < tmp_pUP) ? tmp_pUP : tmp_dUP;
						tmp_q = (tmp_dLEFT < tmp_qLEFT) ? tmp_qLEFT : tmp_dLEFT;
						tmp_d = (tmp_dUPLEFT + match < 0) ? 0 : tmp_dUPLEFT + match;
						tmp_d = (tmp_d < tmp_p) ? tmp_p : tmp_d;
						tmp_d = (tmp_d < tmp_q) ? tmp_q : tmp_d;

						score_l[j] = tmp_d;
					}

					d_buffer[diag_index][current_index] = tmp_d;
					p_buffer[diag_index][current_index] = tmp_p;
					q_buffer[diag_index][current_index] = tmp_q;

				}
			}

			//write out scoree
			short tmp_score = 0;
			max_score_loop: for(int i = 0; i < MAX_DIM; i++){
#pragma HLS PIPELINE II=1
				if(i < diag_len){
					tmp_score = (score_l[i] > tmp_score) ? score_l[i] : tmp_score;
				}
			}

			max_score = (tmp_score > max_score) ? tmp_score : max_score;
			score_index = 0;
			diag_index = (diag_index + 1) % 3;
    	}

    }

    *score = max_score;
}

void readInput(int lenT[INPUT_SIZE], int lenT_local[INPUT_SIZE], char target[INPUT_SIZE][MAX_DIM], char t_local[INPUT_SIZE][MAX_DIM],
		int lenD[INPUT_SIZE], int lenD_local[INPUT_SIZE], char database[INPUT_SIZE][MAX_DIM], char db_local[INPUT_SIZE][MAX_DIM]){

	memcpy(lenT_local, lenT, sizeof(int) * INPUT_SIZE);
	memcpy(lenD_local, lenD, sizeof(int) * INPUT_SIZE);
	memcpy(db_local, database, sizeof(char) * INPUT_SIZE * MAX_DIM);
	memcpy(t_local, target, sizeof(char) * INPUT_SIZE * MAX_DIM);

}

void writeOutput(int score_l[INPUT_SIZE], int score[INPUT_SIZE], int offset, int input_len){

	write_output_loop: for(int iter = 0; iter < input_len; iter++){
        score[iter+offset] = score_l[iter];
    }
}

//////////////////MASTER AXI
extern "C" {
    void sw_maxi(int lenT[INPUT_SIZE], char target[INPUT_SIZE][MAX_DIM], int lenD[INPUT_SIZE], char database[INPUT_SIZE][MAX_DIM],
    int wd, int ws, int gap_opening, int enlargement, int score[INPUT_SIZE],
    int offset, int input_len) {
#pragma HLS INTERFACE s_axilite port=return bundle=control

#pragma HLS INTERFACE m_axi port=lenT bundle=gmem depth=INPUT_SIZE offset=slave
#pragma HLS INTERFACE m_axi port=target bundle=gmem depth=INPUT_SIZE*MAX_DIM offset=slave
#pragma HLS INTERFACE m_axi port=lenD bundle=gmem depth=INPUT_SIZE offset=slave
#pragma HLS INTERFACE m_axi port=database bundle=gmem depth=INPUT_SIZE*MAX_DIM offset=slave
#pragma HLS INTERFACE m_axi port=score bundle=gmem depth=INPUT_SIZE offset=slave

#pragma HLS INTERFACE s_axilite port=lenT bundle=control
#pragma HLS INTERFACE s_axilite port=target bundle=control
#pragma HLS INTERFACE s_axilite port=lenD bundle=control
#pragma HLS INTERFACE s_axilite port=database bundle=control
#pragma HLS INTERFACE s_axilite port=wd bundle=control
#pragma HLS INTERFACE s_axilite port=ws bundle=control
#pragma HLS INTERFACE s_axilite port=gap_opening bundle=control
#pragma HLS INTERFACE s_axilite port=enlargement bundle=control
#pragma HLS INTERFACE s_axilite port=score bundle=control
#pragma HLS INTERFACE s_axilite port=offset bundle=control
#pragma HLS INTERFACE s_axilite port=input_len bundle=control

    //	Defining configuration for scoring
        conf_t local_conf;
        local_conf.match 			= wd;
        local_conf.mismatch			= ws;
        local_conf.gap_opening		= gap_opening + enlargement;
        local_conf.gap_extension	= enlargement;

        int lenT_local[INPUT_SIZE];
        char t_local[INPUT_SIZE][MAX_DIM];
        int lenD_local[INPUT_SIZE];
        char db_local[INPUT_SIZE][MAX_DIM];
        int score_l[INPUT_SIZE];

#pragma HLS DATAFLOW
        readInput(lenT, lenT_local, target, t_local, lenD, lenD_local, database, db_local);

        compute_SW_loop: for (int i = 0; i < INPUT_SIZE; i++) {
			if(i < input_len){
				computeSW(lenT_local[i+offset], t_local[i+offset], lenD_local[i+offset], db_local[i+offset], local_conf, &score_l[i]);
			}
        }

        writeOutput(score_l, score, offset, input_len);
    }
}