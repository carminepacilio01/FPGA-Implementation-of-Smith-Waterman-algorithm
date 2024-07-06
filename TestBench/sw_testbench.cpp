// Test Bench for the implementation of the Smith-Waterman (Gotoh)
#include <random>
#include "smith_waterman.hpp"
#include <time.h>

#define MAX_SEQ_LEN 1000

void printConf(char *seqA, char *seqB, int ws, int wd, int gap_opening, int enlargement);
void fprintMatrix(int *P, int *D, int *Q, int lenA, int lenB);
int compute_golden(int lenA, char *seqA, int lenB, char *seqB, int wd, int ws, int gap_opening, int enlargement);
void random_seq_gen(int lenA, char *seqA, int lenB, char *seqB);
int gen_rnd(int min, int max);

using namespace std;

int main(int argc, char const *argv[]) {
	srandom(2);
	clock_t start, end;

	//	match score
	const int wd = 1;
	//	mismatch score
	const int ws = -1;

	const int gap_opening = -3;
	const int enlargement = -1;

	int count_failed = 0;
	//--------------------------------------------------------------------------------//
	//	executing n tests each with increasing len
	const int n = INPUT_SIZE;
	int lenA[INPUT_SIZE];
	int lenB[INPUT_SIZE];
	char seqA[INPUT_SIZE][MAX_DIM];
	char seqB[INPUT_SIZE][MAX_DIM];

	for(int i = 0; i < n; i++){

		//	generate random sequences
		//	length of the sequences
		lenA[i] = (int)  gen_rnd(MAX_SEQ_LEN - 10, MAX_SEQ_LEN - 2);
		lenB[i] = (int)  gen_rnd(MAX_SEQ_LEN - 10, MAX_SEQ_LEN - 2);

		lenA[i] += 1;
		lenB[i] += 1;

		//	sequences
		//seqA[i] = (char *) malloc(sizeof(char) * (lenA[i]));
		//seqB[i] = (char *) malloc(sizeof(char) * (lenB[i]));

		//	generate rand sequences
		seqA[i][0] = seqB[i][0] = '-';
		seqA[i][lenA[i]] = seqB[i][lenA[i]] = '\0';
		random_seq_gen(lenA[i], seqA[i], lenB[i], seqB[i]);

		//	Printing current configuration
	}

	int score[INPUT_SIZE];
	int golden_score[INPUT_SIZE];
	int n_ops = 0;

	double mean_golden_time = 0;
	double mean_golden_gcup = 0;
	for (int golden_rep = 0; golden_rep < n; golden_rep++) {

		start = clock();
		golden_score[golden_rep] = compute_golden(lenA[golden_rep], seqA[golden_rep], lenB[golden_rep], seqB[golden_rep], wd, ws, gap_opening, enlargement);
		end = clock();

		n_ops = lenA[golden_rep]*lenB[golden_rep];
		double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
		double gcup = (double) (n_ops / cpu_time_used )  * 1e-9;

		mean_golden_time += cpu_time_used;
		mean_golden_gcup += gcup;
	}
	
    printf("Mean Execution time of Host Code: %f ms\n", mean_golden_time / n * 1e3);

	//////////// GCUP
    printf("Mean GCUP Host: %5f GCUPs\n", mean_golden_gcup / n);

	//	computing result using kernel
	sw_maxi(lenA, seqA, lenB, seqB, wd, ws, gap_opening, enlargement, score);
	
	/*	Score results
	for (int i = 0; i < INPUT_SIZE; i++)
	{
		printConf(seqA[i], seqB[i], ws, wd, gap_opening, enlargement);
		cout << "score: " << score[i] << " golden_score: " << golden_score[i] << endl;
		if(score[i] == golden_score[i])
			cout << "TEST " << i + 1 << " PASSED !" << endl;
		else{
			cout << "TEST " << i + 1 << " FAILED !!!" << endl;
			count_failed++;
		}
	}

	cout << "Failed tests: " << count_failed << " Passed tests: " << n - count_failed << endl;*/

	return 0;
}

//	Prints the current configuration
void printConf(char *seqA, char *seqB, int ws, int wd, int gap_opening, int enlargement) {
	cout << endl << "+++++++++++++++++++++" << endl;
	cout << "+ Sequence A: [" << strlen(seqA) << "]: " << seqA << endl;
	cout << "+ Sequence B: [" << strlen(seqB) << "]: " << seqB << endl;
	cout << "+ Match Score: " << wd << endl;
	cout << "+ Mismatch Score: " << ws << endl;
	cout << "+ Gap Opening: " << gap_opening << endl;
	cout << "+ Enlargement: " << enlargement << endl;
	cout << "+++++++++++++++++++++" << endl;
}

int gen_rnd(int min, int max) {
     // Using random function to get random double value
    return (int) min + rand() % (max - min + 1);
}

void random_seq_gen(int lenA, char *seqA, int lenB, char *seqB) {

	int i, j;
	for(i = 1; i < lenA; i++){
		int tmp_gen = gen_rnd(0, 3);
		seqA[i] = (tmp_gen == 0) ? 'A' :
				  (tmp_gen == 1) ? 'C' :
			      (tmp_gen == 2) ? 'G' : 'T';
	}

	for(i = 1; i < lenB; i++){
		int tmp_gen = gen_rnd(0, 3);
		seqB[i] = (tmp_gen == 0) ? 'A' :
				  (tmp_gen == 1) ? 'C' :
				  (tmp_gen == 2) ? 'G' : 'T';
	}
}

int compute_golden(int lenA, char *seqA, int lenB, char *seqB, int wd, int ws, int gap_opening, int enlargement) {
	// Inizializza le matrici D, P, Q
	    std::vector< std::vector<int> > D(lenA, std::vector<int>(lenB, 0));
	    std::vector< std::vector<int> > P(lenA, std::vector<int>(lenB, std::numeric_limits<int>::min() / 2));
	    std::vector< std::vector<int> > Q(lenA, std::vector<int>(lenB, std::numeric_limits<int>::min() / 2));

	    int gap_penalty = gap_opening + enlargement * 1;
	    int max_score = 0;

	    for (int i = 1; i < lenA; ++i) {
	        for (int j = 1; j < lenB; ++j) {
	            // Calcola P[i][j]
	            P[i][j] = std::max(P[i-1][j] + enlargement, D[i-1][j] + gap_penalty);

	            // Calcola Q[i][j]
	            Q[i][j] = std::max(Q[i][j-1] + enlargement, D[i][j-1] + gap_penalty);

	            // Calcola D[i][j]
	            int match = (seqA[i] == seqB[j]) ? wd : ws;
	            D[i][j] = std::max(0, D[i-1][j-1] + match);
	            D[i][j] = std::max(D[i][j], P[i][j]);
	            D[i][j] = std::max(D[i][j], Q[i][j]);

	            // Aggiorna lo score massimo
	            max_score = std::max(max_score, D[i][j]);
	        }
	    }

	    return max_score;
}
