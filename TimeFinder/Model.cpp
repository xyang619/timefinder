//
//  Model.cpp
//  TimeFinder
//
//  Created by young on 2/6/15.
//  Copyright (c) 2015 young. All rights reserved.
//

#include "Model.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace std;

//implement of class Model
double Model::loglik(vector<double> observ, double cutoff) const
{
    double llk = 0;
    for (size_t i = 0; i < observ.size(); ++i) {
        llk += log(dist(observ.at(i), cutoff));
    }
    return llk;
}

//implement of class SpecModel
SpecModel::SpecModel(int T, double m): T_(T), m_(m) { }

void SpecModel::setT(int T)
{
    T_ = T;
}

int SpecModel::getT() const
{
    return T_;
}

double SpecModel::getM() const
{
    return m_;
}

//implement of class HIModel
HIModel::HIModel(int T, double m): SpecModel(T, m) { }

double HIModel::dist(double x, double cutoff) const
{
    double m = getM();
    int T = getT();
    return (1 - m) * T * exp(-(1 - m) * T * (x - cutoff));
}

//implement of class GAModel
GAModel::GAModel(int T, double m): SpecModel(T, m) { }

double GAModel::dist(double x, double cutoff) const
{
    double m = getM();
    int T = getT();
    if(T == 1)
    {
        return (1 - m) * T * exp(-(1 - m) * T * (x - cutoff));
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

////implement of class CGFRModel
CGFRModel::CGFRModel(int T, double m): SpecModel(T, m) { }

double CGFRModel::dist(double x, double cutoff) const
{
    double m = getM();
    int T = getT();
    double u = pow(m, 1.0 / T);
    u = T - (1 - m) * u / (1 - u);
    return u * exp(-u * (x - cutoff));
    
}

//implement of class CGFDModel
CGFDModel::CGFDModel(int T, double m): SpecModel(T, m) { }

double CGFDModel::dist(double x, double cutoff) const
{
    double m = getM();
    int T = getT();
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


GeneralModel::GeneralModel(int K, int T, const std::string & modelfile):K_(K),T_(T)
{
    readModel(modelfile);
    init();
}

GeneralModel::GeneralModel(const GeneralModel & gm)
{
    
}

GeneralModel & GeneralModel::operator=(const GeneralModel & gm){
    return *this;
}

GeneralModel::~GeneralModel(){
    
}

double GeneralModel::dist(double x, double cutoff) const{
    return 0;
}

void GeneralModel::readModel(const std::string & modelfile){
    fstream fin(modelfile.c_str());
    if (!fin.is_open()){
        cout << "Can't read model from file " << modelfile << endl;
        exit(1);
    }
    //read labels
    string label;
    for (int i = 0; i < K_; ++i)
    {
        fin >> label;
        labels.push_back(label);
    }
    //read proportions m_ij
    m_ = new double *[K_];
    for (int i = 0; i < K_; ++i)
    {
        m_[i] = new double[T_];
    }
    for (int t = 0; t < T_; ++t)
    {
        for (int i = 0; i < K_; ++i)
        {
            fin >> m_[i][t];
        }
    }
    fin.close();
}
//equation 1
void GeneralModel::calIt()
{
    I_ = new double [T_];
    for (int t = 0; t < T_; ++t)
    {
        I_[t] = 0;
        for (int i = 0; i < K_; ++i)
        {
            I_[t] += m_[i][t];
        }
        I_[t] = 1 - I_[t];
    }
}

//equation 2
void GeneralModel::calMt(){
    M_ = new double *[K_];
    for (int i = 0; i < K_; ++i) {
        M_[i] = new double[T_];
    }
    for (int i = 0; i < K_; ++i)
    {
        M_[i][0] = M_[i][0];
        for (int t = 1; t < T_; ++t)
        {
            M_[i][t] = M_[i][t - 1] * I_[t] + m_[i][t];
        }
    }
}

//equation 4
void GeneralModel::calHt(){
    h_ = new double *[K_];
    for (int i = 0; i < K_; ++i) {
        h_[i] = new double[T_];
    }
    double tmp = 1;
    for (int t = T_ - 1; t >= 0; --t)
    {
        for (int i = 0; i < K_; ++i)
        {
            h_[i][t] = m_[i][t] * tmp;
        }
        tmp *= I_[t];
    }
}

//equation 3
void GeneralModel::calUt(){
    u_ = new double *[K_];
    for (int i = 0; i < K_; ++i) {
        u_[i] = new double[T_];
    }
    for (int i = 0; i < K_; ++i)
    {
        u_[i][T_ - 1] = 1 - M_[i][T_ - 1];
        for (int t = T_ - 2; t >= 0; --t)
        {
            u_[i][t] = 1 - M_[i][t] + u_[i][t + 1];
        }
    }
}

void GeneralModel::init(){
    calIt();
    calMt();
    calHt();
    calUt();
}

//test cases
int main(int argc, char **argv){
    vector<double> observ;
    observ.push_back(0.1793446);
    observ.push_back(1.9395135);
    observ.push_back(0.1398489);
    observ.push_back(0.2446893);
    observ.push_back(1.7759924);
    observ.push_back(3.7637582);
    observ.push_back(0.8054408);
    observ.push_back(2.2669937);
    observ.push_back(0.5962381);
    observ.push_back(0.4876046);
    Model *hi= new HIModel(5, 0.8);
    Model *ga = new GAModel(5, 0.8);
    Model *cgfr = new CGFRModel(5, 0.8);
    Model *cgfd = new CGFDModel(5, 0.8);
    cout << "x = 0.2; HI f(0.2)= " << hi->dist( 0.2, 0) << "; GA f(0.2)=" << ga->dist(0.2, 0) ;
    cout << "; CGFR f(0.2)= " << cgfr->dist( 0.2, 0) << "; CGFD f(0.2)=" << cgfd->dist(0.2, 0)<<endl;
    cout << "HI llk=" << hi->loglik(observ) << "; GA llk=" << ga->loglik(observ) ;
    cout << "CGFR llk=" << cgfr->loglik(observ) << "; CGFD llk=" << cgfd->loglik(observ) <<endl;
    delete hi;
    delete ga;
    delete cgfr;
    delete cgfd;
}
