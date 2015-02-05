//
//  ModelEval.cpp
//  TimeFinder
//
//  Created by young on 2/5/15.
//  Copyright (c) 2015 young. All rights reserved.
//

#include <cmath>
#include <cstdlib>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include "ModelEval.h"

using namespace std;

//equation 1
void iherP(double **mat, int K, int T, double *it)
{
    for (int t = 0; t < T; ++t)
    {
        it[t] = 0;
        for (int i = 0; i < K; ++i)
        {
            it[t] += mat[i][t];
        }
        it[t] = 1 - it[t];
    }
}

//equation 2
void tProp(double **mat, double *it, int K, int T, double **ap)
{
    for (int i = 0; i < K; ++i)
    {
        ap[i][0] = mat[i][0];
        for (int t = 1; t < T; ++t)
        {
            ap[i][t] = ap[i][t - 1] * it[t] + mat[i][t];
        }
    }
}

//equation 3
void efRec(double **ap, int K, int T, double **u)
{
    for (int i = 0; i < K; ++i)
    {
        u[i][T - 1] = 1 - ap[i][T - 1];
        for (int t = T - 2; t >= 0; --t)
        {
            u[i][t] = 1 - ap[i][t] + u[i][t + 1];
        }
    }
}

//equation 4
void sRate(double **mat, double *it, int K, int T, double **h)
{
    double tmp = 1;
    for (int t = T - 1; t >= 0; --t)
    {
        for (int i = 0; i < K; ++i)
        {
            h[i][t] = mat[i][t] * tmp;
        }
        tmp *= it[t];
    }
}

//equation 8, distribution of general model
double fgm(double x, double *h, double *u, int T, double cutoff)
{
    double num = 0;
    double den = 0;
    for (int t = 0; t < T; ++t)
    {
        double tmp = h[t] * u[t];
        num += tmp * u[t] * exp(-u[t] * x);
        den += tmp * exp(-u[t] * cutoff);
        //cout << h[t] <<" "<<u[t] <<" "<<num <<" "<<den <<endl;
    }
    return num / den;
}

//read model from the file, in which model s
bool readModel(string & filename, int K, int T, double **matrix, vector<string> & labels)
{
    fstream fin(filename.c_str());
    if (!fin.is_open())
    {
        cout << "Can't open file " << filename <<endl;
        return false;
    }
    //read population label
    for (int i = 0; i < K; ++i)
    {
        string label;
        fin >> label;
        labels.push_back(label);
    }
    //read proportions
    for (int t = 0; t < T; ++t)
    {
        for (int i = 0; i < K; ++i)
        {
            fin >> matrix[i][t];
        }
    }
    fin.close();
    return true;
}

double loglik(double *ht, double *ut, int T, double cutoff, const vector<double> & segs)
{
    double llk = 0;
    int size = static_cast<int>(segs.size());
    for (int i = 0; i < size; ++i)
    {
        llk +=log(fgm(segs.at(i), ht, ut, T, cutoff));
    }
    return llk;
}

void help()
{
    cout << "ModelEvaluator" << endl;
    cout << "ModelEvaluator is designed to evaluate the fittness of model and data" << endl;
    cout << "Arguments" << endl;
    cout << "   -h/--help   print help message [optional]" << endl;
    cout << "   -m/--model  parameter file of the alternative model [required]" << endl;
    cout << "   -d/--data   data file of the observation [required]" << endl;
    cout << "   -k/--numb   number of ancestral populations [required]" <<endl;
    cout << "   -g/--gen    generation since admixture [required]" <<endl;
    cout << "   -c/--cutoff cutoff to throw away small segment [optional, default=0]" << endl;
}

int main(int argc, char **argv)
{
    if (argc == 2 && (string(argv[1]) == "-h" || string(argv[1]) == "--help"))
    {
        help();
        exit(0);
    }
    if (argc < 10)
    {
        cerr << "More arguments required, please use -h or --help to get help" << endl;
        exit(1);
    }
    
    string modelfile = "";
    string datafile = "";
    int K = 2;
    int T = 2;
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
        else if (agv == "-m" || agv == "--model")
        {
            modelfile = string(argv[++i]);
        }
        else if (agv == "-d" || agv == "--data")
        {
            datafile = string(argv[++i]);
        }
        else if (agv == "-k" || agv == "--numb")
        {
            K = atoi(argv[++i]);
        }
        else if (agv == "-g" || agv == "--gen")
        {
            T = atoi(argv[++i]);
        }
        else if (agv == "-c" || agv == "--cutoff")
        {
            cutoff = atof(argv[++i]);
        }
    }

    vector<string> labels;
    //read and initialize model
    double **m = new double *[K];
    for (int i = 0; i < K; ++i)
    {
        m[i] = new double[T];
    }
    bool ret = readModel(modelfile, K, T, m, labels);
    if (!ret)
    {
        for (int i = 0; i < K; ++i)
        {
            delete [] m[i];
        }
        delete [] m;
        exit(1);
    }
    double *I = new double[T];
    double **M = new double *[K];
    double **u = new double *[K];
    double **h = new double *[K];
    for (int i = 0; i < K; ++i)
    {
        M[i] = new double[T];
        u[i] = new double[T];
        h[i] = new double[T];
    }
    iherP(m, K, T, I);      //calculate I(t)
    tProp(m, I, K, T, M);   //calcualte M_i(t)
    efRec(M, K, T, u);      //calculate u_i(t)
    sRate(m, I, K, T, h);   //calculate h_i(t)
    
    //reading data
    fstream fin(datafile.c_str());
    if (!fin.is_open())
    {
        cout << "Can't open file " << datafile <<endl;
        exit(1);
    }
    
    map<string, vector<double> > segs;
    map<string, double> lliks;
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
            lliks[label] = 0;
        }
        //greater than or equal to cutoff
        if (len >= cutoff)
        {
            segs.at(label).push_back(len);
        }
    }
    fin.close();
    double totLlik = 0;
    for (size_t i = 0; i < labels.size(); ++i)
    {
        lliks.at(labels.at(i)) = loglik(h[i], u[i], T, cutoff, segs.at(labels.at(i)));
        totLlik += lliks.at(labels.at(i));
    }
    for (size_t i = 0; i < labels.size(); ++i)
    {
        cout << labels.at(i) << "\t";
    }
    cout << "Sum" <<endl;
    for (size_t i = 0; i < labels.size(); ++i)
    {
        cout << lliks.at(labels.at(i)) << "\t";
    }
    cout << totLlik <<endl;
    for (int i = 0; i < K; ++i)
    {
        delete m[i];
        delete M[i];
        delete u[i];
        delete h[i];
    }
    delete [] I;
    delete [] m;
    delete [] M;
    delete [] u;
    delete [] h;

    return 0;
}
