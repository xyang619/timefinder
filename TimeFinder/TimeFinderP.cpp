/*============================================================================
// Name        : TimeFinderP.cpp
// Author      : Young
// Version     :
// Copyright   : GNU GPL V2
// Description : Search best model and estimate parameter
// Using parallel computing
//  Created by young on 1/21/15.
//  Copyright (c) 2015 young. All rights reserved.

 1)calculate mixture proportion m;
 2)calculate the sum of likelihoods for LACS from both ancestral population: 
 L=L(pop1)+L(pop2), under HI, GA, CGFR, CGFD model;
 3)choose the one with highest likelihood as the best model, and the 
 estimated generation as the final result.
============================================================================*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
#include "omp.h"

using namespace std;

const string VERSION "0.0.6p";
void help();
int indexOfMax(double * data, int size);
double sum(vector<double> segs);
double var(vector<double> segs, double mean);
double fhi(double x, double m, int T, double cutoff = 0);
double fga(double x, double m, int T, double cutoff = 0);
double fcgfr(double x, double m, int T, double cutoff = 0);
double fcgfd(double x, double m, int T, double cutoff = 0);
double loglik(double (*fp)(double, double, int, double), vector<double> segs, double m, int T, double cutoff);

int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        cerr << "More arguments required, please use -h or --help to get help" << endl;
        exit(1);
    }
    string fname;
    string outfile = "loglik.txt";
    int maxT = 500;
    double cutoff = 0;
    //dealing with command line arguments
    for (int i = 1; i < argc; ++i)
    {
        string agv(argv[i]);
        if (agv == "-h" || agv == "--help")
        {
            help();
            exit(0);
        }
        else if (agv == "-f" || agv == "--file")
        {
            fname = string(argv[++i]);
        }
        else if (agv == "-t" || agv == "--max")
        {
            maxT = atoi(argv[++i]);
        }
        else if (agv == "-c" || agv == "--cutoff")
        {
            cutoff = atof(argv[++i]);
        }
        else if (agv == "-o" || agv == "--output")
        {
            outfile = argv[++i];
        }
    }

    ifstream fin(fname.c_str());
    if (!fin.is_open())
    {
        cerr << "Can't open file " << fname << endl;
        exit(1);
    }
    //reading the data
    map<string, vector<double> > segs;
    map<string, int> counts;
    map<string, double> sumLens;
    vector<string> labels;

    string line;
    double start, end;
    string label;
    while (getline(fin, line))
    {
        stringstream ss(line);
        ss >> start;
        ss >> end;
        ss >> label;
        double len = end - start;
        //no key not found, then create an new one
        if (segs.find(label) == segs.end())
        {
            vector<double> tmp;
            segs[label] = tmp;
            counts[label] = 0;
            sumLens[label] = 0.0;
            labels.push_back(label);
        }
        //more than 2 labels, close file and exit;
        if (labels.size() > 2)
        {
            cerr << "There are more than two labels in your data, please check" << endl;
            fin.close();
            exit(1);
        }
        counts[label]++;
        sumLens[label] += len;
        if (len >= cutoff)
        {
            segs[label].push_back(len);
        }
    }
    fin.close();

    //calculate the admixture proportions, mean length and its variance
    if (labels.size() == 1)
    {
        cerr << "Only segments from one population found, please check the data" <<endl;
        exit(1);
    }
    double m = sumLens[labels.at(0)] / (sumLens[labels.at(0)] + sumLens[labels.at(1)]);
    string pop1 = "";
    string pop2 = "";
    if (m < 0.5)
    {
        pop1 = labels.at(0);
        pop2 = labels.at(1);
    }
    else
    {
        pop1 = labels.at(1);
        pop2 = labels.at(0);
        m = 1 - m;
    }
    double meanLen1 = sumLens[pop1] / counts[pop1];
    double meanLen2 = sumLens[pop2] / counts[pop2];
    cout << "Populations: " << pop1 << "; " << pop2 << endl;
    cout << "Proportions: " << m << "; " << (1 - m) << endl;
    cout << "Means: " << meanLen1 << "; " << meanLen2 << endl;
    cout << "Variances: " << var(segs[pop1], meanLen1) << "; " << var(segs[pop2], meanLen2) << endl;


    //calculate the likelihood of pop1 as HI, GA, CGFR and CGFD
    double ** llk = new double *[4];
    for (int i = 0; i < 4; ++i)
    {
        llk[i] = new double[maxT];
    }

    int t = 1;
	#pragma omp parallel for
	for (t = 1; t <= maxT; ++t)
	{
		llk[0][t - 1] = loglik(fhi, segs[pop1], m, t, cutoff) + loglik(fhi, segs[pop2], 1 - m, t, cutoff);
		llk[1][t - 1] = loglik(fga, segs[pop1], m, t, cutoff) + loglik(fga, segs[pop2], 1 - m, t, cutoff);
		llk[2][t - 1] = loglik(fcgfr, segs[pop1], m, t, cutoff) + loglik(fcgfd, segs[pop2], m, t, cutoff); //1 as recepient
		llk[3][t - 1] = loglik(fcgfd, segs[pop1], 1 - m, t, cutoff) + loglik(fcgfr, segs[pop2], 1 - m, t, cutoff); //1 as donor
	}

    double *globalMax = new double[4];
    int *index = new int[4];
    string models[4] = {"HI", "GA", "CGFR", "CGFD"};

    for (int i = 0; i < 4; ++i)
    {
        index[i] = indexOfMax(llk[i], maxT);
        globalMax[i] = llk[i][index[i]];
    }
    //print the max likelihood of each column and corresponding generation
    for (int i = 0; i < 4; ++i)
    {
        cout << "Assumed-model: "<< models[i] << "; generation: " << index[i] + 1 << "; likelihood: " << globalMax[i] << endl;
    }
    int mod = indexOfMax(globalMax, 4);
	cout << "Best-model: "<<models[mod] << "; generation: " << index[mod] + 1 << "; likelihood: " << globalMax[mod] << endl; //updated
    
	//save the likelihoods to files
	string outHead = "Generation\tHI\tGA\tCGFR\tCGFD";
    ofstream fout(outfile.c_str());
    if (!fout.good())
    {
        cerr << "Can't open file " << outfile << endl;
    }
    fout << outHead << endl;
    for (int t = 0; t < maxT; ++t)
    {
        fout << setw(6) << t + 1;
        fout << fixed << setprecision(6);
        fout << setw(16) << llk[0][t];
        fout << setw(16) << llk[1][t];
        fout << setw(16) << llk[2][t];
        fout << setw(16) << llk[3][t];
        fout << endl;
    }
    fout.close();

    //clean stuff
    delete[] index;
    delete[] globalMax;
    for (int i = 0; i < 4; ++i)
    {
        delete[] llk[i];
    }
    delete[] llk;
    //test();
    return 0;
}

void help()
{
    cout << "TimeFinder version " << VERSION << endl;
    cout << "TimeFinder is designed to search for the time T when the admixture started," << endl;
    cout << "and figure out the model best fit the data using maximum likelihood estimation" << endl;
    cout << "Arguments" << endl;
    cout << "	-h/--help	print help message [optional]" << endl;
    cout << "	-f/--file	file name for the ancestral segments [required]" << endl;
    cout << "	-t/--max	maximum generation to search for [optional, default=500]" << endl;
    cout << "	-c/--cutoff	cutoff to throw away small segment [optional, default=0]" << endl;
    cout << "	-o/--output	output file for likelihood [optional, default=loglik.txt]" << endl;
}

int indexOfMax(double * data, int size)
{
    double maxVal = data[0];
    int index = 0;
    for (int i = 1; i < size; ++i)
    {
        if (data[i] > maxVal)
        {
            maxVal = data[i];
            index = i;
        }
    }
    return index;
}

double sum(vector<double> segs)
{
    double sum = 0;
    size_t size = segs.size();
    for (size_t i = 0; i < size; ++i)
    {
        sum += segs.at(i);
    }
    return sum;
}

double var(vector<double> segs, double mean)
{
    double sum = 0;
    size_t size = segs.size();
    for (size_t i = 0; i < size; ++i)
    {
        sum += (segs.at(i) - mean) * (segs.at(i) - mean);
    }
    return sum / (size - 1);
}

double fhi(double x, double m, int T, double cutoff)
{
    return (1 - m) * T * exp(-(1 - m) * T * (x-cutoff));
}

double fga(double x, double m, int T, double cutoff)
{
    if(T==1)
    {
        return fhi(x, m, T, cutoff);
    }

    double num = 0;
    double den = 0;
    double tmp = 0;

    for (int t = 2; t <= T; ++t)
    {
        tmp = (T - t + 1) * pow(1 - 1.0 / T, 1 - t);
        num += tmp * (T - t + 1) * exp(-(1 - m) * (T - t + 1) * x);
        den += tmp * exp(-(1 - m) * (T - t + 1) * cutoff);
    }

    num = (1 - m) * (T * T * exp(-(1 - m) * T * x) + num / T);
    den = T * exp(-(1 - m) * T * cutoff) + den / T;

    return num / den;
}

double fcgfr(double x, double m, int T, double cutoff)
{
    double u = pow(m, 1.0 / T);
    u = T - (1 - m) * u / (1 - u);
    return u * exp(-u * (x - cutoff));
}

double fcgfd(double x, double m, int T, double cutoff)
{
    double num = 0;
    double den = 0;
    double tmp1 = 0;

    double tmp2 = pow(m, (T + 1.0) / T);
    double tmp3 = 1 - pow(m, 1.0 / T);
    for (int t = 1; t <= T; ++t)
    {
        tmp1 = pow(m, t * 1.0 / T);
        num += (tmp1 - tmp2) * (tmp1 - tmp2) * exp(-(tmp1 - tmp2) * x / tmp3) / tmp1;
        den += tmp3 * (tmp1 - tmp2) * exp(-(tmp1 - tmp2) * cutoff / tmp3) / tmp1;
    }

    return num / den;
}

double loglik(double (*fp)(double, double, int, double), vector<double> segs, double m, int T, double cutoff)
{
    double ll = 0;
    size_t size = segs.size();
    for (size_t i = 0; i < size; ++i)
    {
        ll += log(fp(segs.at(i), m, T, cutoff));
    }
    return ll;
}
