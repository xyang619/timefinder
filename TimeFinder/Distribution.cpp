#include <cmath>
#include <iostream> 
#include <cstdlib>
using namespace std;

double fhi(double x, double m, int T, double cutoff = 0);
double fga(double x, double m, int T, double cutoff = 0);
double fcgfr(double x, double m, int T, double cutoff = 0);
double fcgfd(double x, double m, int T, double cutoff = 0);
double fgm(double x, double *h, double *u, int k, double cutoff = 0);
void iherP(double **mat, int K, int T, double *it); //eq 1
void tProp(double **mat, double *it, int K, int T, double **ap); //eq 2
void efRec(double **ap, int K, int T, double **u); //eq 3
void sRate(double **mat, double *it, int K, int T, double **h); //eq 4

//HI distribution
double fhi(double x, double m, int T, double cutoff)
{
    return (1 - m) * T * exp(-(1 - m) * T * (x-cutoff));
}

//GA distribution
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

//CGFR distribution
double fcgfr(double x, double m, int T, double cutoff)
{
    double u = pow(m, 1.0 / T);
    u = T - (1 - m) * u / (1 - u);
    return u * exp(-u * (x - cutoff));
}

//CGFD distribution
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

//distribution of general model
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

int main(int argc, char **argv)
{
    if (argc < 4) {
        cout << "Usage: dist <m> <T> <c>" <<endl;
        exit(1);
    }
    int K = 2;
    double m = atof(argv[1]);
    int T = atoi(argv[2]);
    double c = atof(argv[3]);
    double **mat = new double *[K];
    double **ap = new double *[K];
    double **u = new double *[K];
    double **h = new double *[K];
    for (int i = 0; i < K; ++i)
    {
        mat[i] = new double[T];
        ap[i] = new double[T];
        u[i] = new double[T];
        h[i] = new double[T];
    }
    mat[0][0] = m;
    mat[1][0] = 1 - m;
    for (int i = 0; i < K; ++i)
    {
        for (int t = 1; t < T; ++t)
        {
            mat[i][t] = 0;
        }
    }
    double *it = new double[T];
    
    iherP(mat, K, T, it);
    tProp(mat, it, K, T, ap);
    efRec(ap, K, T, u);
    sRate(mat, it, K, T, h);
//    cout <<" m:" <<endl;
//    for (int i=0; i<K; ++i) {
//        for(int t=0; t<T;++t){
//            cout <<mat[i][t] <<" ";
//        }
//        cout <<endl;
//    }
//    cout <<"I(t):";
//    for (int t=0; t<T; ++t) {
//        cout << it[t] << " ";
//    }
//    cout <<endl;
//    cout <<" M:" <<endl;
//    for (int i=0; i<K; ++i) {
//        for(int t=0; t<T;++t){
//            cout <<ap[i][t] <<" ";
//        }
//        cout <<endl;
//    }
//    cout <<" u:" <<endl;
//    for (int i=0; i<K; ++i) {
//        for(int t=0; t<T;++t){
//            cout <<u[i][t] <<" ";
//        }
//        cout <<endl;
//    }
//    cout <<" h:" <<endl;
//    for (int i=0; i<K; ++i) {
//        for(int t=0; t<T;++t){
//            cout <<h[i][t] <<" ";
//        }
//        cout <<endl;
//    }
    cout << "x\tHI\tGA\tCGFR\tCGFD\tGM_HI" <<endl;
    double step = 0.1 / ((1 - m) * T);
    //double x=0.05;
    for (int i = 1; i <= 100; ++i)
    {
        double x = step * i;
        cout << x << "\t";
        cout << fhi(x, m, T, c) << "\t";
        cout << fga(x, m, T, c) << "\t";
        cout << fcgfr(x, m, T, c) << "\t";
        cout << fcgfd(x, 1 - m, T, c) << "\t";
        cout << fgm(x, h[0], u[0], T, c) << endl;
    }
    for (int i = 0; i < K; ++i)
    {
        delete mat[i];
        delete ap[i];
        delete u[i];
        delete h[i];
    }
    delete [] it;
    delete [] mat;
    delete [] ap;
    delete [] u;
    delete [] h;
    return 0;
}
