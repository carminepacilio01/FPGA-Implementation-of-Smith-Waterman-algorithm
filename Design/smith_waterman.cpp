//	Implementation of the Smith-Waterman (Gotoh)
//  Author: Carmine Pacilio
//	Year: 2024
#include "smith_waterman.h"

void computeSW(int lenT, char *target, int lenD, char *database, conf_t scoring, int *score) {

//	copy to local buffers the input sequences
	char target_l[MAX_DIM];
	char database_l[MAX_DIM];
#pragma HLS ARRAY_PARTITION variable=target_l type=complete dim=1
#pragma HLS ARRAY_PARTITION variable=database_l type=complete dim=1

    //copy the target sequence reversed
    for (int i = 0; i < lenT; i++) {
#pragma HLS PIPELINE II=1

        target_l[i] = target[lenT - i - 1];
    }

	memcpy(database_l, database, lenD * sizeof(char));

//	initialize local buffer for the 3 dependency
	int p_buffer[3][MAX_DIM * 2];
	int q_buffer[3][MAX_DIM * 2];
	int d_buffer[3][MAX_DIM * 2];
#pragma HLS ARRAY_PARTITION variable=p_buffer complete dim=1
#pragma HLS ARRAY_PARTITION variable=q_buffer complete dim=1
#pragma HLS ARRAY_PARTITION variable=d_buffer complete dim=1

    init_buffers: for(int index = 0; index < MAX_REP; index++) {
#pragma HLS PIPELINE II=1

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
    const int max_diag_len = (lenT < lenD) ? lenT : lenD; // the maximum value diag_len can reach
    const int t_diag = lenT + lenD - 1;

    int diag_len = 0; // the number of elements on the current diagonal
    int n_diag_repeat = (lenT > lenD) ? (lenT - lenD) + 1 : (lenD - lenT) + 1;
    int diag_index = 0;
    int database_cursor = 0;
    int target_cursor = 0;

    int score_l = 0;

    compute_score: for (int num_diag = 0; num_diag < t_diag; num_diag++) {
#pragma HLS PIPELINE II=1

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
        compute_diagonal: for (int j = 0; j < diag_len; ++j) {
#pragma HLS DEPENDENCE variable=p_buffer inter false
#pragma HLS DEPENDENCE variable=q_buffer inter false
#pragma HLS DEPENDENCE variable=d_buffer inter false
#pragma HLS UNROLL factor=4

            int current_index = database_cursor + j;
            int two_diag = (diag_index == 0) ? 1 : (diag_index == 1) ? 2 : 0;
            int one_diag = (diag_index == 0) ? 2 : diag_index - 1;

            if (target_l[lenT - 1 - (target_cursor - current_index)] == '-') {
                d_buffer[diag_index][current_index] = 0;
                p_buffer[diag_index][current_index] = INFTY;
                q_buffer[diag_index][current_index] = INFTY;
            } else if (database_l[current_index] == '-') {
                d_buffer[diag_index][current_index] = 0;
                q_buffer[diag_index][current_index] = INFTY;
                p_buffer[diag_index][current_index] = INFTY;
            } else {
                p_buffer[diag_index][current_index] = ((d_buffer[one_diag][current_index + UP] + scoring.gap_opening) < (p_buffer[one_diag][current_index + UP] + scoring.gap_extension)) ? p_buffer[one_diag][current_index + UP] + scoring.gap_extension : d_buffer[one_diag][current_index + UP] + scoring.gap_opening;

                q_buffer[diag_index][current_index] = ((d_buffer[one_diag][current_index + LEFT] + scoring.gap_opening) < (q_buffer[one_diag][current_index + LEFT] + scoring.gap_extension)) ? q_buffer[one_diag][current_index + LEFT] + scoring.gap_extension : d_buffer[one_diag][current_index + LEFT] + scoring.gap_opening;

                int match = (target_l[lenT - 1 - (target_cursor - current_index)] == database_l[current_index]) ? scoring.match : scoring.mismatch;

                int tmp_D = (d_buffer[two_diag][current_index + UP_LEFT] + match < 0) ? 0 : d_buffer[two_diag][current_index + UP_LEFT] + match;
                tmp_D = (tmp_D < p_buffer[diag_index][current_index]) ? p_buffer[diag_index][current_index] : tmp_D;
                d_buffer[diag_index][current_index] = (tmp_D < q_buffer[diag_index][current_index]) ? q_buffer[diag_index][current_index] : tmp_D;

                score_l = (score_l > d_buffer[diag_index][current_index]) ? score_l : d_buffer[diag_index][current_index];
            }

            // Aggiorna diag_index per la prossima iterazione del ciclo interno
            if (j + 1 == diag_len) {
                diag_index = (diag_index + 1) % 3;
            }
        }

    }

	//write out score
	memcpy(score, &score_l, sizeof(int));
}

void readInput(int lenT[INPUT_SIZE], int lenT_local[INPUT_SIZE], char target[INPUT_SIZE][MAX_DIM], char t_local[INPUT_SIZE][MAX_DIM],
		int lenD[INPUT_SIZE], int lenD_local[INPUT_SIZE], char database[INPUT_SIZE][MAX_DIM], char db_local[INPUT_SIZE][MAX_DIM]){
#pragma HLS INLINE off

	memcpy(lenT_local, lenT, sizeof(int) * INPUT_SIZE);
	memcpy(lenD_local, lenD, sizeof(int) * INPUT_SIZE);
	memcpy(db_local, database, sizeof(char) * INPUT_SIZE * MAX_DIM);
	memcpy(t_local, target, sizeof(char) * INPUT_SIZE * MAX_DIM);

}

void writeOutput(int score_l[INPUT_SIZE], int score[INPUT_SIZE]){
#pragma HLS INLINE off

	memcpy(score, score_l, sizeof(int) * INPUT_SIZE);
}

//////////////////MASTER AXI
extern "C" {
    void sw_maxi(int lenT[INPUT_SIZE], char target[INPUT_SIZE][MAX_DIM], int lenD[INPUT_SIZE], char database[INPUT_SIZE][MAX_DIM], int wd, int ws, int gap_opening, int enlargement, int score[INPUT_SIZE], int input_len) {
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

        compute_SW_loop: for (int i = 0; i < input_len; i++) {
    #pragma HLS PIPELINE  II=1

        computeSW(lenT_local[i], t_local[i], lenD_local[i], db_local[i], local_conf, &score_l[i]);
        }

        writeOutput(score_l, score);
    }
}

