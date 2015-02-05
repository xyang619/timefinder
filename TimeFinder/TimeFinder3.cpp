/*
 * TimeFinder3.cpp
 *
 *  Created on: Dec 16, 2014
 *      Author: young
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include "omp.h"
#define VERSION "0.0.3"
using namespace std;

//void test();
void help();
int indexOfMax(double * data, int size);
double sum(vector<double> segs);
double var(vector<double> segs, double mean);
double fhi(double x, double m, int T);
double fga(double x, double m, int T);
double fcgf1(double x, double m, int T);
double fcgf2(double x, double m, int T);
double loglik(double (*fp)(double, double, int), vector<double> segs, double m, int T);

int main(int argc, char ** argv) {
    if (argc < 2) {
        cerr << "More arguments required, please use -h or --help to get help" << endl;
        exit(1);
    }
    char *fname;
    string outfile = "loglik.txt";
    int maxT = 1000;
    double cutoff = 0;
    //dealing with command line arguments
    for (int i = 1; i < argc; ++i) {
        string agv(argv[i]);
        if (agv == "-h" || agv == "--help") {
            help();
            exit(0);
        } else if (agv == "-f" || agv == "--file") {
            fname = argv[++i];
        } else if (agv == "-t" || agv == "--max") {
            maxT = atoi(argv[++i]);
        } else if (agv == "-c" || agv == "--cutoff") {
            cutoff = atof(argv[++i]);
        } else if (agv == "-o" || agv == "--output") {
            outfile = argv[++i];
        }
    }
    
    ifstream fin(fname);
    if (!fin.is_open()) {
        cerr << "Can't open file " << fname << endl;
        exit(1);
    }
    //reading the data
    vector<double> segs1;
    vector<double> segs2;
    int size1 = 0;
    int size2 = 0;
    double sum1 = 0;
    double sum2 = 0;
    
    string line;
    double start, end;
    int label;
    while (getline(fin, line)) {
        stringstream ss(line);
        ss >> start;
        ss >> end;
        ss >> label;
        double len = end - start;
        if (label == 1) {
            ++size1;
            sum1 += len;
            if (len >= cutoff) {
                segs1.push_back(end - start);
            }
        } else {
            ++size2;
            sum2 += len;
            if (len >= cutoff) {
                segs2.push_back(end - start);
            }
        }
    }
    fin.close();
    //calculate the admixture proportions, mean length and its variance
    double m = sum1 / (sum1 + sum2);
    double mLen1 = sum1 / size1;
    double mLen2 = sum2 / size2;
    cout << "Proportions of admixture are: " << m << "; " << (1 - m) << endl;
    cout << "Means of Length are: " << mLen1 << "; " << mLen2 << endl;
    //cout << "Variances of Length are: " << var(segs1, mLen1) << "; " << var(segs2, mLen2) << endl;
    
    //calculate the likelihood
    
    double ** llk = new double *[6];
    for (int i = 0; i < 6; ++i) {
        llk[i] = new double[maxT];
    }
    
    int t=1;
#pragma omp parallel for
    for (t = 1; t <= maxT; ++t) {
        llk[0][t - 1] = loglik(fhi, segs1, m, t);
        llk[1][t - 1] = loglik(fhi, segs2, 1 - m, t);
        llk[2][t - 1] = loglik(fga, segs1, m, t);
        llk[3][t - 1] = loglik(fga, segs2, 1 - m, t);
        llk[4][t - 1] = loglik(fcgf1, segs1, m, t);
        llk[5][t - 1] = loglik(fcgf2, segs2, m, t); //proportion here refer to pop1 m !!!
    }
    //write the likelihood to disk file
    ofstream fout(outfile.c_str());
    if (!fout.is_open()) {
        cerr << "Can't open file " << outfile << endl;
    }
    fout << "Generation\tHI_Pop1\tHI_Pop2\tGA_Pop1\tGA_Pop2\tCGF_Pop1\tCGF_Pop2" << endl;
    for (int t = 0; t < maxT; ++t) {
        fout << setw(6) << t + 1;
        fout << fixed << setprecision(6);
        fout << setw(16) << llk[0][t];
        fout << setw(16) << llk[1][t];
        fout << setw(16) << llk[2][t];
        fout << setw(16) << llk[3][t];
        fout << setw(16) << llk[4][t];
        fout << setw(16) << llk[5][t];
        fout << endl;
    }
    fout.close();
    
    //find the maximum likelihood for each column, and corresponding generation
    double *globalMax = new double[6];
    int *index = new int[6];
    for (int i = 0; i < 6; ++i) {
        index[i] = indexOfMax(llk[i], maxT);
        globalMax[i] = llk[i][index[i]];
    }
    //print the max likelihood of each column and corresponding generation
    //for (int i = 0; i < 6; ++i) {
    //cout << "Generation: " << (index[i] + 1) << "; Likelihood: " << globalMax[i] << endl;
    //}
    
    //comparison and choose the best model, two steps
    //1) compare the likelihoods of population 1 under three model,
    //if GA has the highest likelihood, then output the model and generation;
    //2) if is not GA, then compare the
    double *byPop1 = new double[3];
    for (int i = 0; i < 3; ++i) {
        byPop1[i] = globalMax[2 * i];
    }
    int mod = indexOfMax(byPop1, 3);
    delete[] byPop1;
    if (mod == 1) {
        int t = index[2];
        cout << "model is GA, generation is " << (t + 1) << ", likelihood is " << llk[2][t] << endl;
    } else {
        int t1 = index[0];
        int t2 = index[4];
        if (llk[1][t1] >= llk[5][t2]) {
            cout << "model is HI, generation is " << (t1 + 1) << ", likelihood is " << llk[1][t1] << endl;
        } else {
            cout << "model is CGF, generation is " << (t2 + 1) << ", likelihood is " << llk[5][t2] << endl;
        }
    }
    
    //clean stuff
    delete[] index;
    delete[] globalMax;
    for (int i = 0; i < 6; ++i) {
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

void help() {
    cout << "TimeFinder version " << VERSION << endl;
    cout << "TimeFinder is designed to search for the time T when the admixture started," << endl;
    cout << "and figure out the model best fit the data using maximum likelihood estimation" << endl;
    cout << "Arguments" << endl;
    cout << "	-h/--help	print help message [optional]" << endl;
    cout << "	-f/--file	file name for the ancestral segments [required]" << endl;
    cout << "	-t/--max	maximum generation to search for [optional, default=1000]" << endl;
    cout << "	-c/--cutoff	cutoff to throw away small segment [optional, default=0]" << endl;
    cout << "	-o/--output	output file for likelihood [optional, default=loglik.txt]" << endl;
}

int indexOfMax(double * data, int size) {
    double maxVal = data[0];
    int index = 0;
    for (int i = 1; i < size; ++i) {
        if (data[i] > maxVal) {
            maxVal = data[i];
            index = i;
        }
    }
    return index;
}

double sum(vector<double> segs) {
    double sum = 0;
    size_t size = segs.size();
    for (size_t i = 0; i < size; ++i) {
        sum += segs.at(i);
    }
    return sum;
}

double var(vector<double> segs, double mean) {
    double sum = 0;
    size_t size = segs.size();
    for (size_t i = 0; i < size; ++i) {
        sum += (segs.at(i) - mean) * (segs.at(i) - mean);
    }
    return sum / (size - 1);
}

double fhi(double x, double m, int T) {
    return (1 - m) * T * exp(-(1 - m) * T * x);
}

double fga(double x, double m, int T) {
    if(T==1){
        return fhi(x, m, T);
    }
    double tmp = pow(1 - 1.0 / T, T - 1);
    double num = 0;
    for (int t = 2; t <= T; ++t) {
        num += (T - t + 1) * (T - t + 1) * pow(1 - 1.0 / T, 1 - t) * exp(-(1 - m) * (T - t + 1) * x);
    }
    num = T * T * exp(-(1 - m) * T * x) + num / T;
    num *= (1 - m) * tmp;
    
    double den = 0;
    for (int t = 2; t <= T; ++t) {
        den += (T - t + 1) * pow(1 - 1.0 / T, 1 - t);
    }
    den = T + den / T;
    den *= tmp;
    return num / den;
}

double fcgf1(double x, double m, int T) {
    double u = pow(m, 1.0 / T);
    u = (T - (T + 1 - m) * u) / (1 - u);
    return u * exp(-u * x);
}

double fcgf2(double x, double m, int T) {
    double num = 0;
    double tmp1 = 0;
    double tmp2 = pow(m, (T + 1.0) / T);
    double tmp3 = 1 - pow(m, 1.0 / T);
    for (int t = 1; t <= T; ++t) {
        tmp1 = pow(m, t * 1.0 / T);
        num += (tmp1 - tmp2) * (tmp1 - tmp2) * exp(-(tmp1 - tmp2) * x / tmp3) / tmp1;
    }
    double den = 0;
    for (int t = 1; t <= T; ++t) {
        tmp1 = pow(m, t * 1.0 / T);
        den += (tmp1 - tmp2) / tmp1;
    }
    den *= tmp3;
    return num / den;
}

double loglik(double (*fp)(double, double, int), vector<double> segs, double m, int T) {
    double ll = 0;
    size_t size = segs.size();
    for (size_t i = 0; i < size; ++i) {
        ll += log(fp(segs.at(i), m, T));
    }
    return ll;
}