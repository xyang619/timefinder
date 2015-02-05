/*
//  TimeFinder4r3.cpp
//  AdmSimulator
//
//  Created by young on 1/21/15.
//  Copyright (c) 2015 young. All rights reserved.
 
 1) calculate the admixture proportions of each ancestral population, let the ancestral population has minor
 proportion as population 1, and it's mixture proportion is m;
 2) find the generation has mamimum likelihood of LACS from population 1 under assumption HI, GA, CGFR and CGFD
 and record the generations T1(HI), T1(GA), T1(CGFR) and T1(CGFD); likelihoods L1(HI), L1(GA), L1(CGFR) and L1(CGFD);
 3) find the maxmium of L1, and it's corresponding model;
 4) if model in {HI, CGFR}, then find maximum likelihood of population 2 under HI and CGFD, record the generations
 T2(HI) and T2(CGFD), likelihood L2(HI) and L2(CGFD); a）if L2(HI)>L2(CGFD), model is HI; b) else model is CGF, population 1
 is CGFR and 2 is CGFD;
 5) otherwise model in {GA, CGFD}, then find maximum likelihood of population 2 under GA and CGFR, record the generations
 T2(GA) and T2(CGFR), likelihood L2(GA) and L2(CGFR); a) if L2(GA)>L2(CGFR), model is GA; b) else model is CGF, population 1
 is CGFD and 2 is CGFR;
 
 update r3: change the final generation estimation strategy: for HI and GA, output the generation of population 1; for CGF
 model, always output the generation estimated from CGFR.

*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>
//#include "omp.h"
#define VERSION "0.0.4r2"
using namespace std;

//void test();
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
    //vector<double> segs1;
    //vector<double> segs2;
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
        //cout << start <<"; "<<end <<"; " <<label <<"; "<<len<<endl;
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
        //cout << counts[label]<<"; "<<sumLens[label]<<endl;
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
//    cout << labels.at(0) <<"; "<< labels.at(1) << endl;
//    cout << sumLens[labels.at(0)] <<"; "<< sumLens[labels.at(1)] <<endl;
//    cout <<counts[labels.at(0)] <<"; "<<counts[labels.at(1)] <<endl;
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
    cout <<"Pop1: " << pop1 << "; Pop2: " << pop2 << endl;
    cout << "Proportions: " << m << "; " << (1 - m) << endl;
    cout << "Means: " << meanLen1 << "; " << meanLen2 << endl;
    cout << "Variances: " << var(segs[pop1], meanLen1) << "; " << var(segs[pop2], meanLen2) << endl;
    
    
    //calculate the likelihood of pop1 as HI, GA, CGFR and CGFD
    double ** llk = new double *[6];
    for (int i = 0; i < 6; ++i)
    {
        llk[i] = new double[maxT];
    }
    
    int t = 1;
    #pragma omp parallel for
    for (t = 1; t <= maxT; ++t)
    {
        llk[0][t - 1] = loglik(fhi, segs[pop1], m, t, cutoff);
        //llk[1][t - 1] = loglik(fhi, segs[pop2], 1 - m, t, cutoff);
        llk[1][t - 1] = loglik(fga, segs[pop1], m, t, cutoff);
        //llk[3][t - 1] = loglik(fga, segs[pop2], 1 - m, t, cutoff);
        llk[2][t - 1] = loglik(fcgfr, segs[pop1], m, t, cutoff); //1 as recepient
        //llk[5][t - 1] = loglik(fcgfd, segs[pop2], m, t, cutoff); //proportion here refer to receipient m !!! 2 as donor
        //llk[6][t - 1] = loglik(fcgfr, segs[pop2], m, t, cutoff);
        llk[3][t - 1] = loglik(fcgfd, segs[pop1], 1 - m, t, cutoff); //1 as donor
    }
    //Important notes:
    //as the m in CGF model always refer to the recepient, and m here refer to the small one
    //therefore, if 1 is CGFR, then m=m; if 1 is CGFD, then m =1-m
    //find the maximum likelihood for each column, and corresponding generation
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
        //cout << "Generation: " << (index[i] + 1) << "; Likelihood: " << globalMax[i] << endl;
        cout << "Pop1-" << models[i] << ":generation: " << index[i] + 1 << "; Likelihood: " << globalMax[i] << endl;
    }
    
    /*comparison and choose the best model, two steps
    1) compare the likelihoods of population 1 under HI, GA, CGFR, CGFD; 
     if GA has the highest likelihood, output the model and generation;
     moved, find them on the top
    */
    string outHead = "Generation\tHI_Pop1\tGA_Pop1\tCGFR_Pop1\tCGFD_Pop1";

    int mod = indexOfMax(globalMax, 4);
    if (mod == 1 || mod == 3) //GA or CGFD is the candidate
    {
        #pragma omp parallel for
        for (t = 1; t <= maxT; ++t)
        {
            llk[4][t - 1] = loglik(fga, segs[pop2], 1 - m, t, cutoff);
            llk[5][t - 1] = loglik(fcgfr, segs[pop2], 1 - m, t, cutoff);
        }
        outHead = outHead + "\tGA_Pop2\tCGFR_Pop2";
        int t1 = indexOfMax(llk[4], maxT);
        int t2 = indexOfMax(llk[5], maxT);
        double llik1 = llk[4][t1];
        double llik2 = llk[5][t2];
        cout << "Pop2-GA:generation: " << t1 + 1 << "; Likelihood: " << llik1 << endl;
        cout << "Pop2-CGFR:generation: " << t2 + 1 << "; Likelihood: " << llik2 << endl;
        //double llik1 = loglik(fga, segs[pop2], 1 - m, t1, cutoff);
        //double llik2 = loglik(fcgfr, segs[pop2], 1 - m, t2, cutoff);
        if (llik1 > llik2)
        {
            cout << "model: GA; generation: " << index[1] + 1 << "; likelihood: " << globalMax[1] << endl;
        }
        else
        {
            //cout << "model: CGFD; generation: " << index[3] + 1 << "; likelihood: " << globalMax[3] << endl;
            cout << "model: CGFD; generation: " << t2 + 1 << "; likelihood: " << llik2 << endl; //updated
        }
        
    }
    else //HI or CGFR is the candidate
    {
        #pragma omp parallel for
        for (t = 1; t <= maxT; ++t)
        {
            llk[4][t - 1] = loglik(fhi, segs[pop2], 1 - m, t, cutoff);
            llk[5][t - 1] = loglik(fcgfd, segs[pop2], m, t, cutoff);
        }
        outHead = outHead + "\tHI_Pop2\tCGFD_Pop2";
        int t1 = indexOfMax(llk[4], maxT);
        int t2 = indexOfMax(llk[5], maxT);
        double llik1 = llk[4][t1];
        double llik2 = llk[5][t2];
        cout << "Pop2-HI:generation: " << t1 + 1 << "; Likelihood: " << llik1 << endl;
        cout << "Pop2-CGFD:generation: " << t2 + 1 << "; Likelihood: " << llik2 << endl;
        //double llik1 = loglik(fhi, segs[pop2], 1 - m, t1, cutoff);
        //double llik2 = loglik(fcgfd, segs[pop2], m, t2, cutoff);
        if (llik1 > llik2)
        {
            cout << "model: HI, generation: " << index[0] + 1 << ", likelihood: " << globalMax[0] << endl;
        }
        else
        {
            cout << "model: CGFR, generation: " << index[2] + 1 << ", likelihood: " << globalMax[2] << endl;
        }
    }
//    if (mod == 1) {
//        int t = index[mod];
//        cout << "model: GA, generation: " << (t + 1) << ", likelihood: " << llk[mod][t] << endl;
//    }
//    else
//    {
//        if(globalMax[2] >= globalMax [3])
//        { //1 as CGFR is better
//            int t1 = index[0];
//            int t2 = index[2];
//            double llik1 = loglik(fhi, segs[pop2], 1 - m, t1, cutoff);
//            double llik2 = loglik(fcgfd, segs[pop2], m, t2, cutoff);
//            if (llik1 > llik2)
//            {
//                cout << "model: HI, generation: " << (t1 + 1) << ", likelihood: " << llk[0][t1] << endl;
//            }
//            else
//            {
//                cout << "model: CGFR, generation: " << (t2 + 1) << ", likelihood: " << llk[2][t2] << endl;
//            }
//        }
//        else
//        { //1 as CGFD is better
//            if (globalMax[0] > globalMax[3])
//            {
//                int t1 = index[0];
//                cout << "model: HI, generation: " << (t1 + 1) << ", likelihood: " << llk[0][t1] << endl;
//            }
//            else
//            {
//                int t2= index[3];
//                cout << "model: CGFD, generation: " << (t2 + 1) << ", likelihood: " << llk[3][t2] << endl;
//            }
//        }
//    }
    
    //write the likelihood to disk file
    ofstream fout(outfile.c_str());
    if (!fout.is_open())
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
        fout << setw(16) << llk[4][t];
        fout << setw(16) << llk[5][t];
//        fout << setw(16) << llk[6][t];
//        fout << setw(16) << llk[7][t];
        fout << endl;
    }
    fout.close();
    
    //clean stuff
    delete[] index;
    delete[] globalMax;
    for (int i = 0; i < 6; ++i)
    {
        delete[] llk[i];
    }
    delete[] llk;
    //test();
    return 0;
}

//void test() {
//	double m = 0.1;
//	int T = 50;
//	for (int i = 1; i < 100; ++i) {
//		double x = i / 100.0;
//		cout << x << "\t" << fhi(x, m, T) << "\t" << fhi(x, 1 - m, T) << "\t"
//				<< fga(x, m, T) << "\t" << fga(x, 1 - m, T) << "\t"
//				<< fcgf1(x, m, T) << "\t" << fcgf2(x, m, T) << endl;
//	}
//}

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