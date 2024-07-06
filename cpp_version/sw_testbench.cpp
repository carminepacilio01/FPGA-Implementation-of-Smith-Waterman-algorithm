// Test Bench for the implementation of the Smith-Waterman (Gotoh)
#include <random>
#include <time.h>
#include "smith_waterman.hpp"

using namespace std;

void printConf(char *seqA, char *seqB, short ws, short wd, short gap_opening, short enlargement);
void fprintMatrix(int *P, int *D, int *Q, int lenA, int lenB);
short compute_golden(short lenA, char *seqA, short lenB, char *seqB, short wd, short ws, short gap_opening, short enlargement);
void random_seq_gen(int lenA, char *seqA, int lenB, char *seqB);
short gen_rnd(int min, int max);

int main(int argc, char const *argv[]) {
	srandom(time(NULL));

	//	match score
	const short wd = 1;
	//	mismatch score
	const short ws = -1;

	const short gap_opening = -3;
	const short enlargement = -1;

	short count_failed = 0;
	//--------------------------------------------------------------------------------//
	//	executing 10 tests each with increasing len
	for(int i = 0; i < 150; i++){
		short *score = (short *) malloc(sizeof(short));

		//	generate random sequences
		//	length of the sequences
		short lenA = (short) gen_rnd(i+5, MAX_DIM - 10);
		short lenB = (short) gen_rnd(i+5, MAX_DIM - 10);

		//	sequences
		char *seqA = (char *) malloc(sizeof(char) * (lenA + 1));
		char *seqB = (char *) malloc(sizeof(char) * (lenB + 1));

		//	generate rand sequences
		random_seq_gen(lenA, seqA, lenB, seqB);

		//	Printing current configuration
		printConf(seqA, seqB, ws, wd, gap_opening, enlargement);

		short golden_score = compute_golden(lenA, seqA, lenB, seqB, wd, ws, gap_opening, enlargement);

		//	computing result using kernel
		sw_maxi(lenA, seqA, lenB, seqB, wd, ws, gap_opening, enlargement, score);

		//	Score results
		cout << "score: " << *score << " golden_score: " << golden_score << endl;
		if(*score == golden_score)
			cout << "TEST " << i + 1 << " PASSED !" << endl;
		else{
			cout << "TEST " << i + 1 << " FAILED !!!" << endl;
			count_failed++;
		}

		free(seqA);
		free(seqB);
		free(score);
	}

	cout << "Failed tests: " << count_failed << " Passed tests: " << 100 - count_failed << endl;

	return 0;
}

//	Prints the current configuration
void printConf(char *seqA, char *seqB, short ws, short wd, short gap_opening, short enlargement) {
	cout << endl << "+++++++++++++++++++++" << endl;
	cout << "+ Sequence A [" << strlen(seqA) << "]: " << seqA << endl;
	cout << "+ Sequence B: [" << strlen(seqB) << "]: " << seqB << endl;
	cout << "+ Match Score: " << wd << endl;
	cout << "+ Mismatch Score: " << ws << endl;
	cout << "+ Gap Opening: " << gap_opening << endl;
	cout << "+ Enlargement: " << enlargement << endl;
	cout << "+++++++++++++++++++++" << endl;
}

short gen_rnd(int min, int max) {
     // Using random function to get random double value
    return (short) min + (max - min) * (random() % MAX_DIM) / MAX_DIM;
}

void random_seq_gen(int lenA, char *seqA, int lenB, char *seqB) {

	int i, j;
	for(i = 0; i < lenA; i++){
		int tmp_gen = gen_rnd(0, 3);
		seqA[i] = (tmp_gen == 0) ? 'A' :
				  (tmp_gen == 1) ? 'C' :
			      (tmp_gen == 2) ? 'G' : 'T';
	}

	for(i = 0; i < lenB; i++){
		int tmp_gen = gen_rnd(0, 3);
		seqB[i] = (tmp_gen == 0) ? 'A' :
				  (tmp_gen == 1) ? 'C' :
				  (tmp_gen == 2) ? 'G' : 'T';
	}
}

short compute_golden(short lenA, char *seqA, short lenB, char *seqB, short wd, short ws, short gap_opening, short enlargement) {
	short score_internal = 0;

	static short P_internal[MATRIX_SIZE];
	static short D_internal[MATRIX_SIZE];
	static short Q_internal[MATRIX_SIZE];

	const short row = lenA + 1;
	const short col = lenB + 1;
	const short matr_dim = row * col;

    const short gap_penalty = gap_opening + enlargement * 1;

    short i = 0, j = 0, itr = 0;

    //	Initializes first column of the matrices to help value to bound
    for (i = 0; i < col; ++i) {
		P_internal[i] = std::numeric_limits<short>::min() / 2;
		D_internal[i] = 0;
		Q_internal[i] = std::numeric_limits<short>::min() / 2;
	}

    //	Initializes first row of the matrices to help value to bound
    for (j = 0; j < row; ++j) {
		P_internal[j * col] = std::numeric_limits<short>::min() / 2;
		D_internal[j * col] = 0;
		Q_internal[j * col] = std::numeric_limits<short>::min() / 2;
	}

    //	Compute P, D, Q
    for (i = 1; i < row; ++i) {
        for (j = 1; j < col; ++j) {
            //	Calculate P = min{Di-1,j + g(1); Pi-1,j}
            P_internal[i * col + j] = std::max((D_internal[(i - 1) * col + j] + gap_penalty), (P_internal[(i - 1) * col + j] + enlargement));

            //	Calculate Q = min{Di,j-1 + g(1); Qi,j-1}
            Q_internal[i * col + j] = std::max((D_internal[i * col + (j - 1)] + gap_penalty) ,Q_internal[i * col + (j - 1)] + enlargement);

            //	Calculate D = min{Di-1,j-1 + w(ai, bj); Pi,j; Qij}
            //	Calculate the score based on the match or mismatch
            int score = D_internal[(i - 1) * col + (j - 1)] + ((seqA[i - 1] == seqB[j - 1]) ? wd : ws);

            //	Update D[i][j] based on the maximum score among score, P[i][j], and Q[i][j]
            //	Check and update D[i][j] if P[i][j] or Q[i][j] is greater than D[i][j]
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

    return score_internal;
}
