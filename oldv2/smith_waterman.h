#ifndef SMITH_WATERMAN_H
#define SMITH_WATERMAN_H

#define MAX_DIM 256

void smithWaterman(int lenA, char *seqA, int lenB, char *seqB,
		short int *P_out, short int *D_out, short int *Q_out,
		int wd, int ws,int gap_opening, int enlargement);

#endif
