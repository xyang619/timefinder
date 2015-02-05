//
//  ModelEval.h
//  TimeFinder
//
//  Created by young on 2/5/15.
//  Copyright (c) 2015 young. All rights reserved.
//

#ifndef __TimeFinder__ModelEval__
#define __TimeFinder__ModelEval__

#include <string>
#include <vector>

/************************************************************************************
 1) iherP implement equation 1, which return a array of ancestry inherated from
 previous generation;
 2) tProp implement equation 2, which return a matrix of total ancestry at each
 generation for each ancestral population;
 3) efRec implement equation 3, which return a matrix of effective recombination rate
 for the segments from jth ancestral population at generation t;
 4) sRate implement equation 4, which return a matrix of survival rate for the
 segments from jth ancestral population at generation t;
 5) fgm calcualte the density of given x, under general model, with a specific cutoff
*************************************************************************************/

void iherP(double **m, int K, int T, double *I); //eq 1
void tTrop(double **m, double *I, int K, int T, double **M); //eq 2
void efRec(double **M, int K, int T, double **u); //eq 3
void sRate(double **m, double *I, int K, int T, double **h); //eq 4
double fgm(double x, double *h, double *u, int k, double cutoff = 0);
bool readModel(std::string & filename, double **matrix, std::vector<std::string> & labels);
double loglik(double *ht, double *ut, int T, double cutoff, const std::vector<double> & segs);
void help();
#endif /* defined(__TimeFinder__ModelEval__) */
