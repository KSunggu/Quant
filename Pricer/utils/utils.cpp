#include <ql/quantlib.hpp>

#include <stdio.h>

#include <iostream>

#include <ql/utilities/dataparsers.hpp >



#include <fstream>

#include <sstream>

#include <filesystem>



#include <numeric>

#include <algorithm>



#include "utils.h"





void Dataprogress::SchedulePrint(Schedule S) {



    int i = 0;

    /* string str = GET_NAME(S);

    cout << str << endl;

    */



    for (i = 0; i < S.size(); ++i) {

        cout << S[i] << endl;

    };

    cout << "\n";



}

void Dataprogress::SchedulePrint(vector <Date> S) {



    int i = 0;

    /* string str = GET_NAME(S);

    cout << str << endl;

    */



    for (i = 0; i < S.size(); ++i) {

        cout << S[i] << endl;

    };

    cout << "\n";



}



void Dataprogress::luDecomposition(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& U) {

    // Original martrix A, Lower matrix L, Upper matrix U

    int n = A.size();



    L.resize(n, std::vector<double>(n, 0.0));

    U.resize(n, std::vector<double>(n, 0.0));



    for (int i = 0; i < n; i++) {

        for (int j = 0; j < n; j++) {

            U[i][j] = A[i][j];

        }

    }



    for (int k = 0; k < n; k++) {

        int pivotRow = k;

        double maxPivot = std::abs(U[k][k]);



        for (int i = k + 1; i < n; i++) {

            if (std::abs(U[i][k]) > maxPivot) {

                maxPivot = std::abs(U[i][k]);

                pivotRow = i;

            }

        }



        std::swap(U[k], U[pivotRow]);

        std::swap(L[k], L[pivotRow]);



        L[k][k] = 1.0;



        for (int i = k + 1; i < n; i++) {

            L[i][k] = U[i][k] / U[k][k];

            for (int j = k; j < n; j++) {

                U[i][j] -= L[i][k] * U[k][j];

            }

        }

    }



    cout << "LU decomposition Result " << endl;



    cout << "Lower Matrix" << endl;

    for (int i = 0; i < n; ++i) {

        for (int j = 0; j < n; ++j) {

            cout << L[i][j] << " ";

        }

        cout << endl;

    }



    cout << "Upper Matrix" << endl;

    cout << endl;



    for (int i = 0; i < n; ++i) {

        for (int j = 0; j < n; ++j) {

            cout << U[i][j] << " ";

        }

        cout << "\n" << endl;



    }



}



vector <double> Dataprogress::product_MatVec(const vector <vector <double>>& A, vector <double>& v) {



    int rown = A.size();

    int coln = A[0].size();

    double memory;



    vector <double> output;

    // output.resize(rown, vector <double>(coln, 0.0));

    output.resize(rown);



    for (int i = 0; i < rown; ++i) {

        memory = 0;

        for (int j = 0; j < coln; ++j) {

            memory += A[i][j] * v[j];

        }

        output[i] = memory;

    }



    return output;

}
