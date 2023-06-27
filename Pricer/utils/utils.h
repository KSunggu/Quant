#pragma once

#include <ql/quantlib.hpp>

#include <stdio.h>

#include <iostream>

#include <ql/utilities/dataparsers.hpp >



#include <fstream>

#include <sstream>

#include <filesystem>



#include <numeric>

#include <algorithm>



using namespace QuantLib;

using namespace std;

namespace fs = std::filesystem;



class Dataprogress {

public:

    void SchedulePrint(Schedule);

    void SchedulePrint(vector <Date>);



    void luDecomposition(const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);

    vector <double> product_MatVec(const vector <vector <double>>&, vector <double>&);

};

#pragma once
