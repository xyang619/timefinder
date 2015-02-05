/*
 * TimeFinder.cpp
 *
 *  Created on: Dec 16, 2014
 *      Author: young
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#define VERSION "0.0.1"
using namespace std;

//void test();
void help();
int indexOfMax(double * data, int size);
double sum(vector<double> segs);
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
	//dealing with arguments
	for (int i = 1; i < argc; ++i) {
		string agv(argv[i]);
		if (agv == "-h" || agv == "--help") {
			help();
			exit(0);
		} else if (agv == "-f" || agv == "--file") {
			fname = argv[++i];
		} else if (agv == "-t" || agv == "max") {
			maxT = atoi(argv[++i]);
		} else if (agv == "-o" || agv == "output") {
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
	string line;
	double start, end;
	int label;
	while (getline(fin, line)) {
		stringstream ss(line);
		ss >> start;
		ss >> end;
		ss >> label;
		if (label == 1) {
			segs1.push_back(end - start);
		} else {
			segs2.push_back(end - start);
		}
		//cout << start << " " << end << " " << label;
	}
	fin.close();
	//calculate the mixture proportions
	double sum1 = sum(segs1);
	double sum2 = sum(segs2);
	double m = sum1 / (sum1 + sum2);
	double meanLen1 = sum1 / segs1.size();
	double meanLen2 = sum2 / segs2.size();
	cout << "Average Length is: " << meanLen1 << "; " << meanLen2 << endl;
	cout << "Admixture proportion is: " << m << "; " << (1 - m) << endl;
	ofstream fout(outfile.c_str());
	if (!fout.is_open()) {
		cerr << "Can't open file " << outfile << endl;
	}
	double ** llk = new double *[6];
	for (int i = 0; i < 6; ++i) {
		llk[i] = new double[maxT];
	}
	fout << "Generation\tllik_1(HI)\tllik_2(HI)\tllik_1(GA)\tllik_2(GA)\tllik_1(CGF)\tllik_2(CGF)" << endl;
	for (int t = 1; t <= maxT; ++t) {
		fout << t << "\t";
		llk[0][t - 1] = loglik(fhi, segs1, m, t);
		fout << llk[0][t - 1] << "\t";
		llk[1][t - 1] = loglik(fhi, segs2, 1 - m, t);
		fout << llk[1][t - 1] << "\t";
		llk[2][t - 1] = loglik(fga, segs1, m, t);
		fout << llk[2][t - 1] << "\t";
		llk[3][t - 1] = loglik(fga, segs2, 1 - m, t);
		fout << llk[3][t - 1] << "\t";
		llk[4][t - 1] = loglik(fcgf1, segs1, m, t);
		fout << llk[4][t - 1] << "\t";
		llk[5][t - 1] = loglik(fcgf2, segs2, m, t); //proportion always be m !!!
		fout << llk[5][t - 1] << "\t";
		fout << endl;
	}
	fout.close();
	double *globalMax = new double[6];
	int *index = new int[6];
	for (int i = 0; i < 6; ++i) {
		index[i] = indexOfMax(llk[i], maxT);
		globalMax[i] = llk[i][index[i]];
	}
	//print the max likelihood of each column and corresponding generation
	for (int i = 0; i < 6; ++i) {
		cout << "Generation: " << index[i] +1 << "; Likelihood: " << globalMax[i] << endl;
	}
	int col = indexOfMax(globalMax, 6);
	string model;
	switch (col) {
	case 0:
	case 1:
		model = "HI";
		break;
	case 2:
	case 3:
		model = "GA";
		break;
	default:
		model = "CGF";
	}
	cout << "max likelihood found at " << index[col] + 1 << "; value: " << globalMax[col] << "; model: " << model << endl;
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

double fhi(double x, double m, int T) {
	return (1 - m) * T * exp(-(1 - m) * T * x);
}

double fga(double x, double m, int T) {
	if(T==1){
		return fhi(x, m, T);
	}
	double num = 0;
    double den = 0;
    double tmp = 0;
	for (int t = 2; t <= T; ++t) {
        tmp = (T - t + 1) * pow(1 - 1.0 / T, 1 - t);
		num += tmp * (T - t + 1) * exp(-(1 - m) * (T - t + 1) * x);
        den += tmp;
	}
    
	num = T * T * exp(-(1 - m) * T * x) + num / T;
	den = T + den / T;
    
	return (1 - m) * num / den;
}

double fcgf1(double x, double m, int T) {
	double u = pow(m, 1.0 / T);
	u = T - (1 - m) * u / (1 - u);
	return u * exp(-u * x);
}

double fcgf2(double x, double m, int T) {
	double num = 0;
    double den = 0;
	double tmp1 = 0;
	double tmp2 = pow(m, (T + 1.0) / T);
	double tmp3 = 1 - pow(m, 1.0 / T);
	for (int t = 1; t <= T; ++t) {
		tmp1 = pow(m, t * 1.0 / T);
		num += (tmp1 - tmp2) * (tmp1 - tmp2) * exp(-(tmp1 - tmp2) * x / tmp3) / tmp1;
        den += (tmp1 - tmp2) / tmp1;
	}
    den = den * tmp3;
	//den = tmp3 * (T - tmp1);
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
