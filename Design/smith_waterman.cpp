//	Implementation of the Smith-Waterman (Gotoh)
//  Author: Carmine Pacilio
//	Year: 2024
#include "smith_waterman.h"

void computeSW(int lenT, char *target_l, int lenD, char *database_l, conf_t scoring, int *score) {
    

	//write out score
	memcpy(score, &score_l, sizeof(int));
}

void readInput(int lenT[INPUT_SIZE], int lenT_local[INPUT_SIZE], char target[INPUT_SIZE][MAX_DIM], char t_local[INPUT_SIZE][MAX_DIM],
		int lenD[INPUT_SIZE], int lenD_local[INPUT_SIZE], char database[INPUT_SIZE][MAX_DIM], char db_local[INPUT_SIZE][MAX_DIM]){
    
	memcpy(lenT_local, lenT, sizeof(int) * INPUT_SIZE);
	memcpy(lenD_local, lenD, sizeof(int) * INPUT_SIZE);
	memcpy(db_local, database, sizeof(char) * INPUT_SIZE * MAX_DIM);
	memcpy(t_local, target, sizeof(char) * INPUT_SIZE * MAX_DIM);

}

void writeOutput(int score_l[INPUT_SIZE], int score[INPUT_SIZE], int offset, int input_len){

	for(int iter = 0; iter < input_len; iter++){
#pragma HLS LOOP_TRIPCOUNT min=0 max=INPUT_SIZE
#pragma HLS PIPELINE II=1
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

        readInput(lenT, lenT_local, target, t_local, lenD, lenD_local, database, db_local);

        compute_SW_loop: for (int i = 0; i < input_len; i++) {
#pragma HLS LOOP_TRIPCOUNT min=0 max=INPUT_SIZE
#pragma HLS PIPELINE  II=1

        computeSW(lenT_local[i+offset], t_local[i+offset], lenD_local[i+offset], db_local[i+offset], local_conf, &score_l[i]);
        }

        writeOutput(score_l, score, offset, input_len);
    }
}
