#include "PolStr.h"
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
class mandatoryMethods
{
    public:
    mandatoryMethods()
    {}
    double LeftRectangularMethod(int n, double a, double b, double* fxNode)
    {
        double h = (b - a) / n;
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            sum += fxNode[i];
        }
        return sum * h;
    }
    double LeftRectangularMethod(int n, double* xNode, double* fxNode)
    {
        double h;
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            h = xNode[i + 1] - xNode[i];
            sum += fxNode[i] * h;
        }
        return sum;
    }
    pair<double, pair<int, double>> LeftRectangularMethod(int n, double a, double b, double eps, char* fx)
    {
        int currN = n;
        double sum = 0;
        double sumPrev = 0;
        char *pstr = CreatePolStr(fx, 0);
        if (GetError() != 0)
            throw invalid_argument("Polstr error type " + to_string(GetError()));
        int iter = 0;
        while (true)
        {
            sum = 0;
            double h = (b - a) / currN;

            double x = a;
            for (int i = 0; i < currN; i++)
            {
                sum += EvalPolStr(pstr, x);
                x += h;
            }
            sum *= h;
            if (fabs(sum - sumPrev) < eps && iter != 0)
                break;
            sumPrev = sum;
            currN *= 2;
            iter++;
        }
        double finEps = fabs(sum - sumPrev);
        pair<int, double> subPair(iter + 1, finEps);
        pair<double, pair<int, double>> result(sum, subPair);
        return result;            
    }
    double RightRectangularMethod(int n, double a, double b, double* fxNode)
    {
        double h = (b - a) / n;
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            sum += fxNode[i + 1];
        }
        return sum * h;
    }
    double RightRectangularMethod(int n, double* xNode, double* fxNode)
    {
        double h;
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            h = xNode[i + 1] - xNode[i];
            sum += fxNode[i + 1] * h;
        }
        return sum;
    }
    pair<double, pair<int, double>> RightRectangularMethod(int n, double a, double b, double eps, char* fx)
    {
        double sum = 0;
        double sumPrev = 0;
        int currN = n;
        char *pstr = CreatePolStr(fx, 0);
        if (GetError() != 0)
            throw invalid_argument("Polstr error type " + to_string(GetError()));
        int iter = 0;
        while (true)
        {
            sum = 0;
            double h = (b - a) / currN;
            double x = a + h;
            for (int i = 0; i < currN; i++)
            {
                sum += EvalPolStr(pstr, x);
                x += h;
            }
            sum *= h;
            if (fabs(sum - sumPrev) < eps && iter != 0)
                break;
            currN *= 2;
            sumPrev = sum;
            iter++;
        }
        double finEps = fabs(sum - sumPrev);
        pair<int, double> subPair(iter + 1, finEps);
        pair<double, pair<int, double>> result(sum, subPair);
        return result;   
    }
    double trapezoidMethod(int n, double a, double b, double* fxNode)
    {
        double h = (b - a) / n;
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            sum += fxNode[i] + fxNode[i + 1];
        }
        return sum * h / 2;
    }
    double trapezoidMethod(int n, double* xNode, double* fxNode)
    {
        double h;
        double sum = 0;
        for (int i = 0; i < n; i++)
        {
            h = xNode[i + 1] - xNode[i];
            sum += (fxNode[i] + fxNode[i + 1]) * h / 2;
        }
        return sum;
    }
    pair<double, pair<int, double>> trapezoidMethod(int n, double a, double b, double eps, char* fx)
    {
        int currN = n;
        double sum = 0;
        double sumPrev = 0;
        char *pstr = CreatePolStr(fx, 0);
        if (GetError() != 0)
            throw invalid_argument("Polstr error type " + to_string(GetError()));
        int iter = 0;
        while (true)
        {
            sum = 0;
            double h = (b - a) / currN;
            double x = a;
            for (int i = 0; i < currN; i++)
            {
                sum += EvalPolStr(pstr, x) + EvalPolStr(pstr, x + h);
                x += h;
            }
            sum *= h / 2;
            if (fabs(sum - sumPrev) < eps && iter != 0)
                break;
            currN *= 2;
            sumPrev = sum;
            iter++;
        }
        int p = 2;
        double finEps = fabs(sum - sumPrev) / (pow(2, p) - 1);
        pair<int, double> subPair(iter + 1, finEps);
        pair<double, pair<int, double>> result(sum, subPair);
        return result;   
    }
    double SimpsonMethod(int n, double a, double b, double* fxNode)
    {
        double h = (b - a) / n;
        double sum = 0;
        for (int i = 0; i < n - 1; i += 2)
        {
            sum += fxNode[i] + 4 * fxNode[i + 1] + fxNode[i + 2];
        }
        return sum * h / 3;        
    }
    double SimpsonMethod(int n, double* xNode, double* fxNode)
    {
        double sum = 0;
        for (int i = 0; i < n - 1; i += 2)
        {
            double h1 = xNode[i+1] - xNode[i];
            double h2 = xNode[i+2] - xNode[i+1];
            double hSum = h1 + h2;
            double coeff0 = 2 - h2 / h1;
            double coeff1 = (hSum * hSum) / (h1 * h2);
            double coeff2 = 2 - h1 / h2;
            sum += (hSum / 6.0) * (coeff0 * fxNode[i] + coeff1 * fxNode[i+1] + coeff2 * fxNode[i+2]);
        }
        return sum;   
    }
    pair<double, pair<int, double>> SimpsonMethod(int n, double a, double b, double eps, char* fx)
    {
        double sum = 0;
        double sumPrev = 0;
        int currN = n;
        double h = 0;
        char *pstr = CreatePolStr(fx, 0);
        if (GetError() != 0)
            throw invalid_argument("Polstr error type " + to_string(GetError()));
        int iter = 0;
        while (true)
        {
            sum = 0;
            h = (b - a) / currN;
            double x = a;
            for (int i = 0; i < currN; i += 2)
            {
                sum +=  EvalPolStr(pstr, x) + 4 * EvalPolStr(pstr, x + h) + EvalPolStr(pstr, x + 2 * h);
                x += 2 * h;
            }
            sum *= h / 3;
            if (fabs(sum - sumPrev) < eps && iter != 0)
                break;
            currN *= 2;
            sumPrev = sum;
            iter++;
        }
        int p = 4;
        double finEps = fabs(sum - sumPrev) / (pow(2, p) - 1);
        pair<int, double> subPair(iter + 1, finEps);
        pair<double, pair<int, double>> result(sum, subPair);
        return result;  
    }
};
