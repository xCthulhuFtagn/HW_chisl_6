#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

void GoThrough(const vector<vector<double>>& array, const vector<double>& x){
    unsigned n = array.size();
    vector <double> P(n + 1), Q(n + 1), ans;
    P[0] = 0;
    Q[0] = 0;
    unsigned step = -1;
    double A, B, C, D;
    for (int i = 0; i < n; i++)
    {
        if (i == n - 1)
        {
            A = array[n - 1][n - 2];
            B = array[n - 1][n - 1];
            C = 0;
        }
        else if (i == 0)
        {
            A = 0;
            B = array[0][0];
            C = array[0][1];
        }
        else
        {
            A = array[i][step];
            B = array[i][step + 1];
            C = array[i][step + 2];
        }
        step++;
        D = array[i][n-1];
        Q[i + 1] = (D - A * Q[i]) / (B + A * P[i]);
        P[i + 1] = (-C) / (B + A * P[i]);
    }
    for (int i = n - 1; i > 0; --i) {
        Q[i] = Q[i] + P[i] * Q[i + 1];
    }
    for (int i = 1; i < n; ++i) {
        cout << fixed << setprecision(3) << "y(" << x[i - 1] << ") = ";
        cout<<setprecision(10) << Q[i] << endl;
    }
}

void FiniteDifferenceMethod(double a, double b, unsigned n,
    double (*K)(double), double (*L)(double), double (*M)(double), double (*F)(double),
    double R, double S, double T, double V, double W, double Z) {
    double h = (b - a) / (n - 1);
    vector<double> x;
    for (unsigned i = 0; i < n; ++i) {
        x.push_back(a + i * h);
    }
    vector<vector<double>> matrix;
    matrix.resize(n+1);
    for (auto& line : matrix) {
        line.assign(n+1, 0);
    }
    matrix[0][0] = -R / h + S;
    matrix[0][1] = R / h;
    matrix[0][n] = T;
    for (unsigned i = 1; i < n; ++i) {
        double x_i = x[i];
        matrix[i][i - 1] = K(x_i) / pow(h, 2) - L(x_i) / (2 * h);
        matrix[i][i] = -2 * K(x_i) + M(x_i);
        matrix[i][i + 1] = K(x_i) / pow(h, 2) + L(x_i) / (2 * h);
        matrix[i][n] = F(x_i);
    }
    matrix[n][n - 2] = V / h;
    matrix[n][n-1] = -(V / h + W);
    matrix[n][n] = -Z;
    GoThrough(matrix, x);
}

double K(double x) { return 1; }
double L(double x) { return 2*pow(x,2); }
double M(double x) { return 1; }
double F(double x) { return x; }

int main()
{
    double a = 0.5, b = 0.8, step = 0.1, R = -1, S = 2, T = 1, V = 0, W = 1, Z = 3;
        FiniteDifferenceMethod(a, b, (b - a) / step,
            K, L, M, F,
            R, S, T, V, W, Z);
}
