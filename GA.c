/*
 * GA.c
 *
 *  @ José E. Gallardo, Carlos Cotta & Antonio J. Fernández, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/

#include "GA.h"
#include "LABS.h"
#include "random.h"
#include "dynamicMem.h"
#include "timer.h"
#include "TabuSearch.h"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>



#define _isBetterThan_(x,y) ((x) < (y))
#define _hasBetterFitnessThan_(x,y) (_isBetterThan_((x)->fitness, (y)->fitness))






double individualMeritFactor(IndividualP indP, GAP gaP) {
	return meritFactor(indP->fitness, gaP->chromosomeLgth) ;
}

void initIndividual(IndividualP indP, GAP gaP) {
	LABSRandomSeq(indP->chromosome, gaP->chromosomeLgth, gaP->rnd);
}


Fitness evalIndividual(IndividualP indP, GAP gaP) {
	return autocorrelation(indP->chromosome, gaP->chromosomeLgth);
}


void mutateIndividual(IndividualP indP, GAP gaP) {
	for(int i=0; i<gaP->chromosomeLgth; i++) {
		if(withProb(gaP->rnd, gaP->mutationProb))
			indP->chromosome[i] = - indP->chromosome[i];
	}
}


IndividualP newIndividual(GAP gaP) {
	IndividualP indP = allocate(sizeof(Individual));
	indP->chromosome = allocate(gaP->chromosomeLgth*sizeof(Gene));
	return indP;
}

void freeIndividual(IndividualP indP) {
	deallocate(indP->chromosome);
	deallocate(indP);
}

void copyIndividualFromTo(IndividualP from, IndividualP to, GAP gaP) {
	to->fitness = from->fitness;
	memcpy(to->chromosome,from->chromosome,gaP->chromosomeLgth*sizeof(Gene));
}


void crossover(IndividualP parent1P, IndividualP parent2P, IndividualP childP, GAP gaP) {
	for(int i=0; i<gaP->chromosomeLgth; i++) 
		childP->chromosome[i] = (RandomNextBool(gaP->rnd) ? parent1P : parent2P)->chromosome[i];
}

bool isOptimal(IndividualP indP, GAP gaP) {
	return indP->fitness == knownOptimum(gaP->chromosomeLgth);
}

GAP newGA(int id, int populationSz, int chromosomeLgth) {

	GAP gaP = allocate(sizeof(GA));
	gaP->GAId = id;	
	gaP->populationSz = populationSz;
	gaP->chromosomeLgth = chromosomeLgth;
	
	gaP->population = allocate(populationSz*sizeof(IndividualP));
	for(int i=0; i<populationSz; i++) {
		gaP->population[i] = newIndividual(gaP);
	}
	
	gaP->child = newIndividual(gaP);
	gaP->aux = newIndividual(gaP);
	
	gaP->best = newIndividual(gaP);
	gaP->best->fitness = INT_MAX;
	
	return gaP;
}

void freeGA(GAP gaP) {
	for(int i=0; i<gaP->populationSz; i++) {
		freeIndividual(gaP->population[i]);
	}	
	freeIndividual(gaP->child);
	freeIndividual(gaP->best);
	freeIndividual(gaP->aux);

	deallocate(gaP);
}


// Tournament selection. Returns index into population
int tournament(GAP gaP, int tournamentSz) {
	int bestIndividual = RandomNextIntUntil(gaP->rnd,gaP->populationSz);
	for(int i=0; i<tournamentSz; i++) {
		int bestIndividual2 = RandomNextIntUntil(gaP->rnd,gaP->populationSz);
		if(_hasBetterFitnessThan_(gaP->population[bestIndividual2],gaP->population[bestIndividual])) 
			bestIndividual = bestIndividual2;
	}	
	return bestIndividual;
}

// selects an individual. Returns index into population
int randomIndividual(GAP gaP) {
	return RandomNextIntUntil(gaP->rnd, gaP->populationSz);
}


bool areEqIndividuals(GAP gaP, IndividualP indP1, IndividualP indP2) {
	if(indP1->fitness != indP2->fitness)
		return false;
	else 
		return memcmp(indP1->chromosome, indP2->chromosome, gaP->chromosomeLgth) == 0;
} 

bool isAlreadyInPopulation(GAP gaP, IndividualP indP) {
	for(int i=0; i<gaP->populationSz; i++) {
		if(areEqIndividuals(gaP, gaP->population[i], indP))
			return true;
	}
	return false;
}

void updateBestIndividual(IndividualP ind, GAP gaP) {
	if(_hasBetterFitnessThan_(ind, gaP->best)) {
	
		copyIndividualFromTo(ind, gaP->best, gaP);
		
		printf("\nGA(%2d - %2d) iter: %12d  %4d  %.2f  %.1f(secs)", gaP->GAId, gaP->seed, gaP->iter, gaP->best->fitness, individualMeritFactor(gaP->best, gaP), GetElapsedTime(gaP->timerP));
		
		if(isOptimal(gaP->best, gaP))
			gaP->stop = true; 
	}	
}

// Initializes all individuals in population
void initPopulation(GAP gaP) {
	int bestIndividualIdx = 0;
	for(int i=0; i<gaP->populationSz; i++) {
		initIndividual(gaP->population[i], gaP);
		gaP->population[i]->fitness = evalIndividual(gaP->population[i], gaP);

		if(_hasBetterFitnessThan_(gaP->population[i], gaP->population[bestIndividualIdx]))
			bestIndividualIdx = i;
	}
	updateBestIndividual(gaP->population[bestIndividualIdx], gaP);
}

void *RunGA(void *ptr) {
	GAP gaP = (GAP) ptr;

	gaP->iter = 0;	
	initPopulation(gaP);

	if(!isOptimal(gaP->best, gaP))
	  while(!(gaP->stop)) {
		// reproduction
		if(withProb(gaP->rnd, gaP->crossoverProb)) {
			int parentIdx1 = tournament(gaP, gaP->tournamentSz);
			int parentIdx2 = tournament(gaP, gaP->tournamentSz);

			crossover(gaP->population[parentIdx1], gaP->population[parentIdx2], gaP->child, gaP);
			
		} else 
			copyIndividualFromTo(gaP->population[randomIndividual(gaP)], gaP->child, gaP);
			
		// mutation
		mutateIndividual(gaP->child, gaP);
		
		//gaP->child->fitness = evalIndividual(gaP->child, gaP);
		// Local Search
		gaP->child->fitness = TabuSearch(gaP->child->chromosome, gaP->chromosomeLgth, gaP->rnd, gaP->timerP);

		// replacement	
		if(!isAlreadyInPopulation(gaP, gaP->child)) {
			int forDeathIdx = randomIndividual(gaP);
			IndividualP aux = gaP->child;
			gaP->child = gaP->population[forDeathIdx];
			gaP->population[forDeathIdx] = aux;
			updateBestIndividual(gaP->child, gaP);
		}
		(gaP->iter)++;
		
		if((GetElapsedTime(gaP->timerP) > gaP->maxExecTime))
			gaP->stop = true;
	}
}

int main(int argc, char *args[]) {

	setbuf(stdout, NULL);
	installBackTracer();

	if(argc != 5) {
		printf("usage: %s <binary seq length> <random seed> <max time (secs)> <no. of threads>\n",args[0]);
		exit(1);
	}
	
	readLABSOpts("LABSOpts.txt");
	
	int size = atoi(args[1]);
	int seed = atoi(args[2]);
	double maxExecTime = atoi(args[3]);
	int threads = atoi(args[4]);

	// Create random generators
	Random rnds[threads];
	for(int i=0; i<threads; i++) {
		rnds[i] = newRandom(i, seed);
	}
	
	printf("\nLABS %d (opt=%d,%.2f), seed=%d, maxExecTime=%.1f(secs), threads=%d", size, knownOptimum(size), meritFactor(knownOptimum(size), size),seed, maxExecTime, threads);

	TimerP timerP = newTimer();
	StartTimer(timerP);	

	GAP gasP[threads];
	for(int i=0; i<threads; i++) {
		gasP[i] = newGA(i, 100, size);
		gasP[i]->seed = seed;
		gasP[i]->rnd = rnds[i];
		gasP[i]->timerP = timerP;
		gasP[i]->crossoverProb = 0.9;
		gasP[i]->mutationProb = 1.0/(double)size;
		gasP[i]->tournamentSz = 2; // binary tournament
		gasP[i]->maxExecTime = maxExecTime;
		gasP[i]->stop = false;
	}


	int threadId[threads];
	pthread_t workerThreads[threads];
	
	
	for(int i=0; i<threads; i++) {
		threadId[i] = pthread_create(&(workerThreads[i]), NULL, RunGA, (void*) gasP[i]);
	}
	
	for(int i=0; i<threads; i++) {
		pthread_join(workerThreads[i], NULL);
	}
	
	StopTimer(timerP);

	for(int i=0; i<threads; i++) {
		printf("\nGA(%2d - %2d) BEST: %d %.2f: ", gasP[i]->GAId, gasP[i]->seed, gasP[i]->best->fitness, individualMeritFactor(gasP[i]->best,gasP[i]));
		fprintSeq(stdout, gasP[i]->best->chromosome, gasP[i]->chromosomeLgth);
	}
	printf("\n");
	
	for(int i=0; i<threads; i++) {	
		freeGA(gasP[i]);
		freeRandom(&rnds[i]);
	}
	freeTimer(&timerP);
}
