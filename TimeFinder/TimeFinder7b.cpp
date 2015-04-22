/*============================================================================
 // Name        : TimeFinderB.cpp
 // Author      : Young
 // Version     :
 // Copyright   : GNU GPL V2
 // Description : Search best model and estimate parameter
 // Using binary search
 //  Created by young on 1/21/15.
 //  Copyright (c) 2015 young. All rights reserved.
 
 1)calculate mixture proportion m;
 2)calculate the sum of likelihoods for LACS from both ancestral population:
 L=L(pop1)+L(pop2), under HI, GA, CGFR, CGFD model;
 3)choose the one with highest likelihood as the best model, and the
 estimated generation as the final result.
 Update 7: 2015-4-22
 Consider chromosome length L, distributions and likelihood are all modified
 ============================================================================*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <cstdlib>

using namespace std;

const string VERSION = "0.0.7b";
const double kInf = exp(1000);
//void test();
void help();
int indexOfMax(double * data, int size);
double sum(vector<double> segs);
double var(vector<double> segs, double mean);
double fhi(double x, double m, int T, double cutoff = 0, double length = kInf);
double fga(double x, double m, int T, double cutoff = 0, double length = kInf);
double fcgfr(double x, double m, int T, double cutoff = 0, double length = kInf);
double fcgfd(double x, double m, int T, double cutoff = 0, double length = kInf);
double loglik(double (*fp)(double, double, int, double, double), vector<double> segs, double m, int T, double cutoff, double length);
int bSearchT(double (*fp1)(double, double, int, double, double), double (*fp2)(double, double, int, double, double), vector<double> segs1, vector<double> segs2, double m1, double m2, double cutoff, double length, int maxT, double &llk);

int main(int argc, char ** argv)
{
    if (argc < 2)
    {
        cerr << "More arguments required, please use -h or --help to get help" << endl;
        exit(1);
    }
    
    string fname = "";
    int maxT = 500;
    double cutoff = 0;
    double length = kInf;
    
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
        else if(agv == "-l" || agv == "--length"){
            length = atof(argv[++i]);
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
        if (len > cutoff && len < length)
        {
            segs[label].push_back(len);
        }
    }
    fin.close();
    
    //calculate the admixture proportions, average and variance of length
    if (labels.size() == 1)
    {
        cerr << "Only segments from one population found, please check the data" << endl;
        exit(1);
    }
    double m = sumLens[labels.at(0)] / (sumLens[labels.at(0)] + sumLens[labels.at(1)]);
    string pop1 = "";
    string pop2 = "";
    //decide which one is population 1
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
    cout <<"Populations: " << pop1 << "; " << pop2 << endl;
    cout << "Proportions: " << m << "; " << (1 - m) << endl;
    cout << "Means: " << meanLen1 << "; " << meanLen2 << endl;
    cout << "Variances: " << var(segs[pop1], meanLen1) << "; " << var(segs[pop2], meanLen2) << endl;
    
    
    //calculate the likelihood of pop1 under assumption of HI, GA, CGFR and CGFD
    double *globalMax = new double[4];
    int *index = new int[4];
    string models[4] = {"HI", "GA", "CGFR", "CGFD"};
    
    index[0] = bSearchT(fhi, fhi, segs[pop1], segs[pop2], m, 1 - m, cutoff, length, maxT, globalMax[0]);
    index[1] = bSearchT(fga, fga, segs[pop1], segs[pop2], m, 1 - m, cutoff, length, maxT, globalMax[1]);
    index[2] = bSearchT(fcgfr, fcgfd, segs[pop1], segs[pop2], m, m, cutoff, length, maxT, globalMax[2]);
    index[3] = bSearchT(fcgfd, fcgfr, segs[pop1], segs[pop2], 1 - m, 1 - m, cutoff, length, maxT, globalMax[3]);
    
    //print the max likelihood under each assumption and corresponding generation
    for (int i = 0; i < 4; ++i)
    {
        cout << "Assumed-model: " << models[i] << "; generation: " << index[i] << "; likelihood: " << globalMax[i] << endl;
    }
    
    //comparison and choose the best model
    
    int mod = indexOfMax(globalMax, 4);
    cout << "Best-model: " << models[mod] << "; generation: " << index[mod] << "; likelihood: " << globalMax[mod] << endl;
    //clean stuff
    delete[] index;
    delete[] globalMax;
    
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
    cout << "	-l/--length	chromosome length [optional, default=Infinity]" << endl;
}

//find the max
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

//Summation
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

//calculate variance
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

//HI distribution
double fhi(double x, double m, int T, double cutoff, double length)
{
    double u = (1-m)*T;
    return u* exp(-u * x)/(exp(-u*cutoff)-exp(-u*length));
}

//GA distribution
double fga(double x, double m, int T, double cutoff, double length)
{
    if(T==1)
    {
        return fhi(x, m, T, cutoff, length);
    }
    
    double num = 0;
    double den = 0;
    double tmp = 0;
    
    for (int t = 2; t <= T; ++t)
    {
        tmp = (T - t + 1) * pow(1 - 1.0 / T, 1 - t);
        num += tmp * (T - t + 1) * exp(-(1 - m) * (T - t + 1) * x);
        den += tmp * (exp(-(1 - m) * (T - t + 1) * cutoff)-exp(-(1 - m) * (T - t + 1) * length));
    }
    
    num = (1 - m) * (T * T * exp(-(1 - m) * T * x) + num / T);
    den = T * (exp(-(1 - m) * T * cutoff)-exp(-(1 - m) * T * length)) + den / T;
    
    return num / den;
}

//CGFR distribution
double fcgfr(double x, double m, int T, double cutoff, double length)
{
    double u = pow(m, 1.0 / T);
    u = T - (1 - m) * u / (1 - u);
    return u * exp(-u * x)/(exp(-u*cutoff)-exp(-u*length));
}

//CGFD distribution
double fcgfd(double x, double m, int T, double cutoff, double length)
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
        den += tmp3 * (tmp1 - tmp2) * (exp(-(tmp1 - tmp2) * cutoff / tmp3)-exp(-(tmp1 - tmp2) * length / tmp3)) / tmp1;
    }
    
    return num / den;
}

//Likelihood function
double loglik(double (*fp)(double, double, int, double, double), vector<double> segs, double m, int T, double cutoff, double length)
{
    double ll = 0;
    size_t size = segs.size();
    for (size_t i = 0; i < size; ++i)
    {
        ll += log(fp(segs.at(i), m, T, cutoff,length));
    }
    return ll;
}

//MLE by binary search
int bSearchT(double (*fp1)(double, double, int, double, double), double (*fp2)(double, double, int, double, double), vector<double> segs1, vector<double> segs2, double m1, double m2, double cutoff, double length, int maxT, double &llk)
{
    if (maxT < 1)
    {
        cerr << " Max search range should larger than 1" << endl;
        exit(1);
    }
    if (maxT == 1){
        cerr << " Only one generation to compare, don't ask me, you know it" << endl;
        llk = 0;
        return maxT;
    }
    if (maxT == 2)
    {
        double llk1 = loglik(fp1, segs1, m1, 1, cutoff, length) + loglik(fp2, segs2, m2, 1, cutoff, length);
        double llk2 = loglik(fp1, segs1, m1, 2, cutoff, length) + loglik(fp2, segs2, m2, 2, cutoff, length);
        if (llk1 > llk2)
        {
            llk = llk1;
            return 1;
        }
        else
        {
            cerr << "Warning: likelihood still increasing, please increase maxT" << endl;
            llk = llk2;
            return maxT;
        }
    }
    
    int tleft = 1;
    int tright = maxT;
    while (1)
    {
        double lleft, lmiddle, lright;
        int tmiddle = (tleft + tright) / 2;
        lmiddle = loglik(fp1, segs1, m1, tmiddle, cutoff, length) + loglik(fp2, segs2, m2, tmiddle, cutoff, length);
        lleft = loglik(fp1, segs1, m1, tmiddle - 1, cutoff, length) + loglik(fp2, segs2, m2, tmiddle - 1, cutoff, length);
        lright = loglik(fp1, segs1, m1, tmiddle + 1, cutoff, length) + loglik(fp2, segs2, m2, tmiddle + 1, cutoff, length);
        //cerr << tmiddle << ": " << lleft << "-" << lmiddle << "-" << lright << endl;
        if (lmiddle >= lleft && lmiddle >= lright)
        {
            llk = lmiddle;
            return tmiddle;
        }
        else if (lleft > lmiddle)
        {
            if (tmiddle == 2)
            {
                llk = lleft;
                return 1;
            }
            else
            {
                tright = tmiddle;
            }
        }
        else if (lmiddle < lright)
        {
            if (tmiddle == (maxT - 1)) {
                cerr << "Warning: ikelihood still increasing, please increase maxT" << endl;
                llk = lright;
                return maxT;
            }
            else
            {
                tleft = tmiddle;
            }
        }
    }
}
