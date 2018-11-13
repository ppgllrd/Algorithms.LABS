/*
 * TabuSearch.h
 *
 *  @ Jos� E. Gallardo, Carlos Cotta & Antonio J. Fern�ndez, 2012
 *  Finding Low Autocorrelation Binary Sequences with Memetic Algorithms.
 *  Applied Soft Computing. 9(4): 1252-1262 (2009).
*/
 
#ifndef TABU_SEARCH_H_
#define TABU_SEARCH_H_

#include "types.h"
#include "random.h"
#include "timer.h"

Fitness TabuSearch(Gene *initialSeq, int N, Random rnd, Timer *timer);

#endif /* TABU_SEARCH_H_ */
