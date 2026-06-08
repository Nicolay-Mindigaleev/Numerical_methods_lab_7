#include "PolStr.h"
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
pair<double, double*> ChebyshevMethod(int k, double a, double b, char* fx)
{
    if (k == 8 || k >= 10)
        throw runtime_error("System hasn't valid solves");
    double t[9][9] = {
        {0.0},                                                          // k=1
        {0.577350, -0.577350},                                          // k=2
        {0.707107, 0.0, -0.707107},                                    // k=3
        {0.794654, 0.187592, -0.187592, -0.794654},                    // k=4
        {0.832497, 0.374541, 0.0, -0.374541, -0.832497},               // k=5
        {0.866246, 0.422518, 0.266635, -0.266635, -0.422518, -0.866246}, // k=6
        {0.883861, 0.529656, 0.323913, 0.0, -0.323913, -0.529656, -0.883861}, // k=7
        {},
        {0.911579, 0.601029, 0.528761, 0.167905, 0.0, -0.167905, -0.528761, -0.601018, -0.911579} // k=9
    };
    double c = (b - a) / k;
    double sum = 0;
    char *pstr = CreatePolStr(fx, 0);
    if (GetError() != 0)
        throw invalid_argument("Polstr error type " + to_string(GetError()));
    double* xArr = new double[k];
    for (int i = 0; i < k; i++)
    {
        double tNode = t[k - 1][i];
        xArr[i] = (b + a) / 2 + ((b - a) / 2) * tNode;
    }
    for (int i = 0; i < k; i++)
    {
        sum += EvalPolStr(pstr, xArr[i]);
    }
    pair<double, double*> result(sum * c, xArr);
    return result;
}