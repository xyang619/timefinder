//
//  WaveScanner.h
//  TimeFinder
//
//  Created by young on 2/4/15.
//  Copyright (c) 2015 young. All rights reserved.
//

#ifndef __TimeFinder__WaveScanner__
#define __TimeFinder__WaveScanner__

#include <vector>

/*********************************************************************************
 **Class ParamExp to store the parameters of EM
 **Here assume the observations are from the combination of K sets of numbers
 **which following exponential distributions, with parameters lambda[j], and the
 **mix proportions prop[j], the purpose of EM is used to estimate the parameters
 **lambda[j] and mix proportions[j], where j in [1,2,...,K]
 ********************************************************************************/
class ParamExp
{
public:
    ParamExp(int K);
    ParamExp(int K, double *lambda, double *prop);
    ParamExp(const ParamExp & rhs);
    ParamExp & operator=(const ParamExp & rhs);
    int getK() const;
    double getLambda(int index) const;
    double getProp(int index) const;
    //void setK(int index);
    //void setLambda(int, double);
    //void setProp(int, double);
    //bool isDistinct() const;
    bool isConverge(const ParamExp & par);
    void print();
    virtual ~ParamExp();
private:
    int K;
    double * lambda;
    double * prop;
};

/*********************************************************************************
 Class EMExp to perform EM algorithm in parameters estimation with hidden variable
 The EMExp class has three attributes:
 1) parameters to be estimated;
 2) and a sequences of observations, which used to estimate the parameters;
 3) the likelihood corresponding to the observations and parameters.
 
 The method iterate is used to perform EM algorithm, including two steps:
 E-Step:
 1) calculate p_j=prob(z_i=j|x_i, theta_t)
 =prob(x_i|z_i=j,theta_t)*prob(z_i=j,theta_t)/(sum_j{from 1 to K}(.)
 . denotes the numerator
 2) calculate p_j*xi ; j=1,...,K
 theta_t: parameters at the t-th iteration
 M-Step:
 1) update prop, refer as m, m_j=sum(p_j)/sum(p_1+p_2+...+p_K)=sum(p_j)/n
 2) update lambda, lambda_j=sum(p_j)/sum(p_j*x_i)
 note: update the parameters for t+1 times iteration
 *********************************************************************************/

class EMExp
{
public:
    EMExp(const ParamExp & par, const std::vector<double> & observ);
    double getLik() const;
    ParamExp & getPar();
    void updateLik();
    void iterate(int maxIter);
    virtual ~EMExp();
private:
    double lik;
    ParamExp par;
    std::vector<double> observ;
};

/*********************************************************************************
 Auxiliary functions
 summation over a set of data
*********************************************************************************/
double sum(double * data, int size);
inline double max(double a, double b);
void help();

#endif /* defined(__TimeFinder__WaveScanner__) */
