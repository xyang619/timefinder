//
//  WaveScanner.cpp
//  TimeFinder
//
//  Created by young on 2/4/15.
//  Copyright (c) 2015 young. All rights reserved.
//

#include "WaveScanner.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;
const double kDelta = 0.000001; //converge condition
const double kminP = 0.01;      //minimum proportion for a wave
const string VERSION = "0.0.1";

//implementation of class ParamExp
//constructor
ParamExp::ParamExp(int K) : K(K)
{
    lambda = new double[K];
    for (int j = 0; j < K; ++j)
    {
        lambda[j] = 1.0 * rand() / RAND_MAX;
    }
    prop = new double[K];
    double tmp = 0;
    for (int j = 0; j < (K - 1); ++j)
    {
        prop[j] = rand() / (1.0 * K * RAND_MAX);
        tmp += prop[j];
    }
    prop[K - 1] = 1 - tmp;
}

//another constructor
ParamExp::ParamExp(int K, double *l, double *p) : K(K)
{
    lambda = new double[K];
    prop = new double[K];
    for (int j = 0; j < K; ++j)
    {
        lambda[j] = l[j];
        prop[j] = p[j];
    }
}

//copy constructor
ParamExp::ParamExp(const ParamExp & rhs)
{
    K = rhs.K;
    lambda = new double[K];
    prop = new double[K];
    for (int j = 0; j < K; ++j)
    {
        lambda[j] = rhs.lambda[j];
        prop[j] = rhs.prop[j];
    }
}

//overloading operator assignment
ParamExp & ParamExp::operator=(const ParamExp & rhs)
{
    if (this != &rhs) {
        K = rhs.K;
        if (lambda != NULL)
            delete[] lambda;
        if (prop != NULL)
            delete[] prop;
        lambda = new double[K];
        prop = new double[K];
        for (int j = 0; j < K; ++j)
        {
            lambda[j] = rhs.lambda[j];
            prop[j] = rhs.prop[j];
        }
    }
    return *this;
}

int ParamExp::getK() const
{
    return K;
}

double ParamExp::getLambda(int index) const
{
    return lambda[index];
}

double ParamExp::getProp(int index) const
{
    return prop[index];
}
//TODO
//bool ParamExp::isDistinct() const
//{
//    return false;
//}

//test converge or not
bool ParamExp::isConverge(const ParamExp & par)
{
    bool converge = 1;
    for (int i = 0; i < K; ++i)
    {
        if (abs(lambda[i] - par.getLambda(i)) > kDelta
            || abs(prop[i] - par.getProp(i)) > kDelta)
        {
            converge = 0;
            break;
        }
    }
    return converge;
}

//output the results
void ParamExp::print()
{
    cout << "par=(";
    for (int j = 0; j < K; ++j)
    {
        cout << lambda[j] << ", " << prop[j] << "; ";
    }
    cout << ")" << endl;
}

//deconstructor
ParamExp::~ParamExp()
{
    if (lambda != NULL)
        delete[] lambda;
    if (prop != NULL)
        delete[] prop;
}
//end of class ParamExp

//implementation of class EMExp
//constructor
EMExp::EMExp(const ParamExp & par, const std::vector<double> & observ) : par(par), observ(observ)
{
    updateLik();
}

double EMExp::getLik() const
{
    return lik;
}

ParamExp & EMExp::getPar()
{
    return par;
}

//update likelihood
void EMExp::updateLik()
{
    double tmp = 0;
    for (size_t i = 0; i < observ.size(); ++i)
    {
        int k = par.getK();
        double fval = 0;
        for (int j = 0; j < k; ++j)
        {
            double m = par.getProp(j);
            double l = par.getLambda(j);
            fval += m * l * exp(-l * observ.at(i));
        }
        tmp += log(fval);
    }
    lik = tmp;
}

//EM iteration, converge or reach max iteration time will terminate
void EMExp::iterate(int maxIter)
{
    int it = 0;
    int kval = par.getK();
    int size = static_cast<int>(observ.size());
    double *nlambda;
    double *nprop;
    double **pval;
    double **pxval;
    nlambda = new double[kval];
    nprop = new double[kval];
    pval = new double *[kval];
    pxval = new double *[kval];
    for (int j = 0; j < kval; ++j)
    {
        pval[j] = new double[size];
        pxval[j] = new double[size];
    }
    while (it < maxIter)
    {
        //E-step
        for (int i = 0; i < size; ++i)
        {
            double denormator = 0;
            for (int j = 0; j < kval; ++j)
            {
                double mj = par.getProp(j);
                double lj = par.getLambda(j);
                pval[j][i] = mj * lj * exp(-lj * observ.at(i));
                denormator += pval[j][i];
            }
            for (int j = 0; j < kval; ++j)
            {
                pval[j][i] /= denormator;
                pxval[j][i] = pval[j][i] * observ.at(i);
            }
        }
        //M-step
        for (int j = 0; j < kval; ++j)
        {
            double sump = sum(pval[j], size);
            nprop[j] = sump / size;
            nlambda[j] = sump / sum(pxval[j], size);
        }
        ParamExp updatedPar(kval, nlambda, nprop);
        if (par.isConverge(updatedPar))
        {
            par = updatedPar;
            updateLik();
            break;
        }
        else
        {
            par = updatedPar;
            updateLik();
        }
        //cout << "Iteration " << ++it << " --> llk: " << getLik() << "; ";
        //getPar().print();
    }
    
    //clean stuff
    delete[] nlambda;
    delete[] nprop;
    for (int j = 0; j < kval; ++j)
    {
        delete[] pval[j];
        delete[] pxval[j];
    }
    delete[] pval;
    delete[] pxval;
}

EMExp::~EMExp()
{
}
//end of class EMExp

//summation over data
double sum(double * data, int size)
{
    double tmp = 0;
    for (int i = 0; i < size; ++i)
    {
        tmp += data[i];
    }
    return tmp;
}

inline double max(double a, double b){
    return  a > b ? a : b ; 
}

void help() {
    cout << "WaveScanner Version " << VERSION << endl;
    cout << "WaveScanner is used to scan the number of waves and estimate parameters" << endl;
    cout << "Arguments:" << endl;
    cout << "   -h/--help   print help message [optional]" << endl;
    cout << "   -f/--file   file name of observations [required]" << endl;
    cout << "   -m/--iter   max iteration to perform EM [optional, 10000]" << endl;
    cout << "   -c/--cutoff cutoff of minimum length [optional, default, 0.0]" << endl;
}

int main(int argc, char ** argv)
{
    if(argc < 2)
    {
        cerr << "Need more arguments than provided, please use -h/--help to get help" << endl;
        exit(1);
    }
    string filename = "";
    int maxIter = 10000;
    double cutoff = 0;
    //dealing with command line arguments
    for (int i = 1; i < argc; ++i)
    {
        string arg(argv[i]);
        if (arg == "-h" || arg == "--help")
        {
            help();
            exit(0);
        }
        else if (arg == "-f" || arg == "--file")
        {
            filename = string(argv[++i]);
        }
        else if (arg == "-c" || arg == "--cutoff")
        {
            cutoff = atof(argv[++i]);
        }
        else if (arg == "-m" || arg == "--iter")
        {
            maxIter = atoi(argv[++i]);
        }
    }
    if(filename.size()==0){
        cerr <<"File name required, please check help" <<endl;
        exit(1);
    }
    srand(time(NULL));
    vector<double> observ;
    //Reading read and store as a map
    ifstream fin(filename.c_str());
    if (!fin.is_open()) {
        cout << "Can't open file " << filename << endl;
        exit(1);
    }
    
    map<string, vector<double> > segs;
    map<string, double> sumLens;
    vector<string> labels;
    double tLen = 0;
    string line;
    while (getline(fin, line))
    {
        stringstream ss(line);
        double start, end;
        string label;
        ss >> start;
        ss >> end;
        ss >> label;
        double len = end - start;
        //no key not found, then create an new one
        if (segs.find(label) == segs.end())
        {
            vector<double> tmp;
            segs[label] = tmp;
            sumLens[label] = 0.0;
            labels.push_back(label);
        }
        sumLens.at(label) += len;
        tLen += len;
        //greater than or equal to cutoff
        if (len >= cutoff)
        {
            segs.at(label).push_back(len);
        }
    }
    fin.close();
    map<string, double> props;
    int numLabel = static_cast<int>(labels.size());
    #pragma omp parallel for
    for (int i = 0; i < numLabel; ++i)
    {
        string label = labels.at(i);
        props[label] = sumLens.at(label) / tLen;
        cout << "Label: " << label << "; Proption: " << props.at(label) <<endl;
        int maxK = static_cast<int>(props.at(label) / kminP);
        bool findBest = false;
        int k = 1;
        //ParamExp parPrev(1);
        double llkPrev = 0;
        while (!findBest)
        {
            ParamExp parCur(k);
            EMExp em(parCur, segs.at(label));
            //em.updateLik();
            em.iterate(maxIter);
            cout << "K=" << k <<"; Likelihood=" << setprecision(8) << em.getLik() <<"; ";
            em.getPar().print();
            if (abs((em.getLik() - llkPrev) / em.getLik()) < kDelta)
            {
                findBest = true;
            }
            else
            {
                ++k;
                llkPrev = em.getLik();
            }
            if (k > maxK)
            {
                cout << "No Optimal found" <<endl;
                break;
            }
        }
    }

    return 0;
}
