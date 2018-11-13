/*
 * LABS.c
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/


#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "LABS.h"
#include "random.h"


Fitness autocorrelation(Gene *seq, int N) {
	Fitness s, f = 0;

	for(int i=0; i<N-1; i++) {
		s = 0;
		for(int j=0; j<N-1-i; j ++)
			s +=  seq[j]*seq[i+1+j];
		f += (SQR(s));
	}
	return f;
}


double meritFactor(Fitness f, int N) {
	return SQR(N) / (2.0 * f);
}

double seqMeritFactor(Gene *seq, int N) {
	return meritFactor(autocorrelation(seq,N), N);
}

void LABSRandomSeq(Gene *seq, int N, Random rnd) {
	for(int i=0; i<N; i++) {
		seq[i] = RandomNextBool(rnd) ? 1 : (-1);
	}
}


void fprintSeq(FILE *fp, Gene *seq, int N) {
	for(int i=0; i<N; i++)
		fprintf(fp, "%c", seq[i] == 1 ? '+' : '-');
	fprintf(fp, " %d %.2f", autocorrelation(seq, N), seqMeritFactor(seq, N));
}

#define MAX_OPTS  	60
Fitness labsOpts[1+MAX_OPTS];


void readLABSOpts(char* fileName) {
	FILE* f;

	f = fopen(fileName, "rt");
    if(!f) {
		printf("\nCan't open %s file for reading", fileName);
        exit(1);
    }
	
	int size, opt;
	while(fscanf(f, "%d %d", &size, &opt) != EOF) {
		labsOpts[size] = opt;
	}
	
	fclose(f);
}


Fitness knownOptimum(int N) {
	return ((N <= MAX_OPTS) ? labsOpts[N] : 0);
}


