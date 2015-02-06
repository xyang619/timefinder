//
//  Model.h
//  TimeFinder
//
//  Created by young on 2/6/15.
//  Copyright (c) 2015 young. All rights reserved.
//

#ifndef __TimeFinder__Model__
#define __TimeFinder__Model__

#include <vector>
#include <string>

//base class for all models
class Model
{
public:
    virtual ~Model() {};
    virtual double dist(double x, double cutoff) const = 0;
    double loglik(std::vector<double> observ, double cutoff = 0) const;
};

//base class for specific model HI, GA, CGFR and CGFD
class SpecModel: public Model
{
public:
    SpecModel(int T, double m);
    virtual ~SpecModel() {};
    void setT(int T);
    int getT() const;
    double getM() const;
private:
    int T_;
    double m_;
};

//HI model
class HIModel: public SpecModel
{
public:
    HIModel(int T = 1, double m = 0.5);
    ~HIModel() {};
    double dist(double x, double cutoff = 0) const;
};

//GA model
class GAModel: public SpecModel
{
public:
    GAModel(int T = 1, double m = 0.5);
    ~GAModel() {}
    double dist(double x, double cutoff = 0) const;
};

//CGFR model
class CGFRModel: public SpecModel
{
public:
    CGFRModel(int T = 1, double m = 0.5);
    ~CGFRModel() {}
    double dist(double x, double cutoff = 0) const;
};

//CGFD model
class CGFDModel: public SpecModel
{
public:
    CGFDModel(int T = 1, double m = 0.5);
    ~CGFDModel() {}
    double dist(double x, double cutoff = 0) const;
};

//General model
class GeneralModel: public Model
{
public:
    GeneralModel(int K, int T, const std::string & modelfile);
    GeneralModel(const GeneralModel & gm);
    GeneralModel & operator=(const GeneralModel & gm);
    ~GeneralModel();
    double dist(double x, double cutoff = 0) const;
private:
    void init();
    void calIt();
    void calMt();
    void calHt();
    void calUt();
    void readModel(const std::string & modelfile);
    int K_;
    int T_;
    std::vector<std::string> labels;
    double *I_;
    double **m_;
    double **M_;
    double **h_;
    double **u_;
};
#endif /* defined(__TimeFinder__Model__) */
