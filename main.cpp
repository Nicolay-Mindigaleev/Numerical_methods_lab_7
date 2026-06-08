#include <iostream>
#include <fstream>
#include "mandatoryMethods.cpp"
#include "ChebyshevMethod.cpp"
#include <iomanip>
using namespace std;
int main()
{
    int method;
    string grid = "";
    int n;
    int k;
    double a, b;
    string s = "";
    char fx[128] = "";
    double eps;
    ifstream input("input.txt");
    input >> method;
    double* xNode = nullptr;
    double* fxNode = nullptr;
    if (method != 5)
    {
        input >> grid;
        input >> n;
        if (grid != "dynamic") 
        {
            xNode = new double[n + 1];
            fxNode = new double[n + 1];
        }
    }
    if (method == 4 && n % 2 != 0)
        throw runtime_error("For Simpson method n must be even");
    if (method == 5)
        input >> k;
    if (grid == "uniform" || grid == "dynamic" || method == 5)
        input >> a >> b;
    else if (grid == "nonuniform")
    {
        for (int i = 0; i <= n; i++)
        {
            input >> xNode[i];
        }
    }
    if (method != 5)
        input >> s;
    if (s == "table")
    {
        for (int i = 0; i <= n; i++)
            input >> fxNode[i];
    }
    else if (s == "analitical" || method == 5)
    {
        input >> fx;
        if (grid == "dynamic")
            input >> eps;
    }
    else
    {
        throw runtime_error("s must be table or analitical. Check txt file");
    }
    mandatoryMethods methods;
    double result = 0;
    int iter = -1;
    double FinEps = -1;
    double* xArr = nullptr;
    switch(method)
    {
        case 1:
            if (grid == "uniform") {
                if (s == "table")
                    result = methods.LeftRectangularMethod(n, a, b, fxNode);
                else
                    throw runtime_error("s must be table");
            }
            else if (grid == "nonuniform") {
                if (s == "table")
                    result = methods.LeftRectangularMethod(n, xNode, fxNode);
                else
                    throw runtime_error("s must be table");
            }
            else if (grid == "dynamic") {
                auto res = methods.LeftRectangularMethod(n, a, b, eps, fx);
                result = res.first;
                iter = res.second.first;
                FinEps = res.second.second;
            }
            break;
            
        case 2:
            if (grid == "uniform") {
                if (s == "table")
                    result = methods.RightRectangularMethod(n, a, b, fxNode);
                else
                    throw runtime_error("s must be table");
            }
            else if (grid == "nonuniform") {
                if (s == "table")
                    result = methods.RightRectangularMethod(n, xNode, fxNode);
                else
                    throw runtime_error("s must be table");
            }
            else if (grid == "dynamic") {
                auto res = methods.RightRectangularMethod(n, a, b, eps, fx);
                result = res.first;
                iter = res.second.first;
                FinEps = res.second.second;
            }
            break;
            
        case 3:
            if (grid == "uniform") {
                if (s == "table")
                    result = methods.trapezoidMethod(n, a, b, fxNode);
                else
                    throw runtime_error("s must be table");
            }
            else if (grid == "nonuniform") {
                if (s == "table")
                    result = methods.trapezoidMethod(n, xNode, fxNode);
                else
                    throw runtime_error("s must be table");
            }
            else if (grid == "dynamic") {
                auto res = methods.trapezoidMethod(n, a, b, eps, fx);
                result = res.first;
                iter = res.second.first;
                FinEps = res.second.second;
            }
            break;
            
        case 4:
            if (grid == "uniform") {
                if (s == "table")
                    result = methods.SimpsonMethod(n, a, b, fxNode);
                else
                    throw runtime_error("s must be table");
            }
            else if (grid == "nonuniform") {
                if (s == "table")
                    result = methods.SimpsonMethod(n, xNode, fxNode);
                else
                    throw runtime_error("s must be table");
            }
            else if (grid == "dynamic") {
                auto res = methods.SimpsonMethod(n, a, b, eps, fx);
                result = res.first;
                iter = res.second.first;
                FinEps = res.second.second;
            }
            break;
        case 5:
            {
                auto res = ChebyshevMethod(k, a, b, fx);
                result = res.first;
                xArr = res.second;
            }
            break;
        default:
            throw invalid_argument("Case error: invalid method value");
    }  
    if (grid != "dynamic") 
    {
        delete[] xNode;
        delete[] fxNode;
    }
    ofstream output("output.txt");
    if (grid == "dynamic")
    {
        double digits = (int)-log10(eps) + 0.5;
        if (digits < 0) digits = 0;
        double formatedResult = round(result * pow(10, digits)) / pow(10, digits);
        output <<fixed << setprecision(digits)<<formatedResult<<endl;
    }
    else
    {
        output << setprecision(10) << result << endl;
    }
    if (iter != -1)
        output<<iter<<endl; 
    if (FinEps != -1)
    {
        output<<scientific<<FinEps<<endl;   
    }
    if (method == 5)
    {
        for (int i = 0; i < k; i++)
        {
            output << xArr[i] << "\t";
        }
        delete[] xArr;
    }  
}

