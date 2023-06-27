/* C++ 17++ Based code*/ // Visual Basic 속성 변경을 통해 수정할 것.

// FT * (S - K) 꼴이므로 S 는 해외 Money account가 risk neutral인 확률분포에 대하여 진행. 즉, 기초자산에는 콴토조정이 없고, 환율만 상관관계를 입혀서 path Genration 할 것

// Doestic CCY : USD



#include <stdio.h>

#include <iostream>

#include <ql/quantlib.hpp>

#include <ctime>

#include <fstream>

#include <filesystem>

#include <sstream>

#include <numeric>

#include <algorithm>

#include <ql/utilities/dataparsers.hpp >



#define GET_NAME(name) #name



using namespace QuantLib;

using namespace std;

namespace fs = std::filesystem;



bool PrintSwitch = false;

string Inputpath;

string Outputpath;

vector <Date> daily; vector <double> daily_tenors;

vector <Rate> kr_spot; vector <double> kr_fwds;

vector <double> us_spot; vector <double> us_fwds;



vector <string> Open_getline(string);

vector <string> Open_getline(string, bool);

vector <double> product_MatVec(const vector <vector <double>>&, vector <double>&);

void SchedulePrint(vector <Date>);

void SchdeulePrint(Schedule);

void luDecomposition(const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);



int main() {

    // ************************* Path 지정

    fs::path(fs::current_path());

    if (fs::is_directory("Inputs")) {

        Inputpath = "./Inputs/";

    }



    if (fs::is_directory("Outputs")) { Outputpath = "./Outputs/"; }

    else { fs::create_directories("Outputs"); Outputpath = "./Outputs/"; }





    Date Base = Date(31, Mar, 2023);

    Date Trade = Date(23, Nov, 2019);

    Date Term = Date(17, May, 2023);



    for (int i = Base.serialNumber(); i <= Term.serialNumber(); ++i) {

        daily.push_back(Date(i));

        if (i == Base.serialNumber() || i == Term.serialNumber()) { std::cout << Date(i) << endl; }

    }





    double Notional = 60000000000;

    double fxsig = 0.0983329798956531;

    double FX0 = 1303.8;

    double S0 = 322.03;

    double sig = 0.186385084757187;

    double cor = -0.690181395519357;

    double q = 2.21 * 0.01;









    Calendar FXcalendar = SouthKorea();

    Calendar calendar = JointCalendar(SouthKorea(), UnitedStates(UnitedStates::NYSE));

    BusinessDayConvention BCC = Following;

    DayCounter PayoffDCC = Actual365Fixed();



    Settings::instance().evaluationDate() = Base;



    // *****************  Curve **********************

    // Curve for USDSOFR

    std::cout << "Read Curve Data..." << endl;

    ifstream file1(Inputpath + "USD_SOFR.txt");

    vector <int> USD_Day;

    vector <string> USD_STD;

    vector <Rate> USD_Rates;

    if (file1.is_open()) {

        vector <int> colum_idx; // Example: Read column number

        colum_idx.push_back(1); colum_idx.push_back(2); colum_idx.push_back(5);



        int flag = 0;

        string line;

        while (std::getline(file1, line)) {

            std::istringstream iss(line);

            std::string column;

            std::vector<std::string> columns;



            while (std::getline(iss, column, '\t')) {

                columns.push_back(column);

            }



            if (flag > 0) {

                for (int i = 0; i < colum_idx.size(); ++i) {

                    if (colum_idx[i] <= columns.size()) {

                        // Do something with the desired column

                        int idx;

                        idx = colum_idx[i] - 1;

                        if (i == 0) {

                            int value = stod(columns[idx]);

                            USD_Day.push_back(value);

                        }

                        else if (i == 1) { USD_STD.push_back(columns[idx]); }

                        else {

                            Rate value = stod(columns[idx]);

                            value *= 0.01;

                            USD_Rates.push_back(value);

                        }

                        //cout << "Column " << colum_idx[i] << ": " << columns[idx] << "\t";



                    }

                    else {

                        std::cout << "Invalid column number." << std::endl;

                        break;

                    }

                }

                // cout << "\n";

            }

            flag++;

        }

        file1.close();

    }

    else {

        std::cout << "Failed to open the file." << std::endl;

        return 1; // Return from the program indicating an error.

    }

    vector <Period> USD_Tenor;

    for (int i = 0; i < USD_Day.size(); ++i) {

        vector <string> STD;

        STD.push_back("WK"); STD.push_back("MO"); STD.push_back("YR");

        if (USD_STD[i].compare(STD[0]) == 0) { USD_Tenor.push_back(Period(USD_Day[i], Weeks)); }

        else if (USD_STD[i].compare(STD[1]) == 0) { USD_Tenor.push_back(Period(USD_Day[i], Months)); }

        else if (USD_STD[i].compare(STD[2]) == 0) { USD_Tenor.push_back(Period(USD_Day[i], Years)); }

        // cout << USD_Tenor[i] << endl;

    }



    Integer settlement = 2;

    Calendar USDcalendar = UnitedStates(UnitedStates::FederalReserve);

    Integer pay_delay = 2;

    DayCounter DCC = Actual360();

    // Period coupon_period = Period(6, Months);

    BusinessDayConvention BCC1 = ModifiedFollowing; // Fixed

    BusinessDayConvention BCC2 = ModifiedFollowing; // Float



    /*

    helpers1.push_back(boost::make_shared<OISRateHelper>(settlement, USD_Tenor[i],

             Handle <Quote>(boost::make_shared <SimpleQuote>(USD_Rates[i])),

             Index1,

             Handle<YieldTermStructure> {},

             false,

             pay_delay,

             BCC2,

             Annual,

             USDcalendar))

    */



    boost::shared_ptr Index1 = boost::make_shared <OvernightIndex>(OvernightIndex("USD SOFR", settlement, USDCurrency(), USDcalendar, DCC));

    vector <boost::shared_ptr <RateHelper>> helpers1;

    for (int i = 0; i < USD_Tenor.size(); ++i) {

        helpers1.push_back(boost::make_shared<OISRateHelper>(settlement, USD_Tenor[i],

            Handle <Quote>(boost::make_shared <SimpleQuote>(USD_Rates[i])),

            Index1,

            Handle<YieldTermStructure> {},

            false,

            pay_delay,

            BCC2,

            Annual,

            USDcalendar));

    }

    std::cout << "Curve check " << endl;

    boost::shared_ptr <YieldTermStructure> USDSOFR(new PiecewiseYieldCurve <ZeroYield, Cubic, IterativeBootstrap>(Base, helpers1, DCC));

    USDSOFR->enableExtrapolation();



    for (const auto& date : daily) {

        us_spot.push_back(USDSOFR->zeroRate(date, DCC, Continuous));

        us_fwds.push_back(USDSOFR->forwardRate(date, NullCalendar().advance(date, Period(1, Days)), DCC, Continuous));

    }

    us_fwds.pop_back();



    // Curve for KRWCD91

    ifstream file2(Inputpath + "KRW_CD.txt");

    vector <int> KRW_Day;

    vector <string> KRW_STD;

    vector <Rate> KRW_Rates;

    if (file2.is_open()) {

        vector <int> colum_idx; // Example: Read column number

        colum_idx.push_back(1); colum_idx.push_back(2); colum_idx.push_back(5);







        int flag = 0;

        string line;

        while (std::getline(file2, line)) {

            std::istringstream iss(line);

            std::string column;

            std::vector<std::string> columns;



            while (std::getline(iss, column, '\t')) {

                columns.push_back(column);

            }



            if (flag > 0) {

                for (int i = 0; i < colum_idx.size(); ++i) {

                    if (colum_idx[i] <= columns.size()) {

                        int idx;

                        idx = colum_idx[i] - 1;

                        if (i == 0) {

                            int value = stod(columns[idx]);

                            KRW_Day.push_back(value);

                        }

                        else if (i == 1) { KRW_STD.push_back(columns[idx]); }

                        else {

                            Rate value = stod(columns[idx]);

                            value *= 0.01;

                            KRW_Rates.push_back(value);

                        }

                        //cout << "Column " << colum_idx[i] << ": " << columns[idx] << "\t";



                    }

                    else {

                        std::cout << "Invalid column number." << std::endl;

                        break;

                    }

                }

                //cout << "\n";

            }

            flag++;

        }

        file2.close();

    }

    else {

        std::cout << "Failed to open the file." << std::endl;

        return 1; // Return from the program indicating an error.

    }

    vector <Period> KRW_Tenor;

    for (int i = 0; i < KRW_Day.size(); ++i) {

        vector <string> STD;

        STD.push_back("WK"); STD.push_back("MO"); STD.push_back("YR");

        if (KRW_STD[i].compare(STD[0]) == 0) { KRW_Tenor.push_back(Period(KRW_Day[i], Weeks)); }

        else if (KRW_STD[i].compare(STD[1]) == 0) { KRW_Tenor.push_back(Period(KRW_Day[i], Months)); }

        else if (KRW_STD[i].compare(STD[2]) == 0) { KRW_Tenor.push_back(Period(KRW_Day[i], Years)); }

        // cout << KRW_Tenor[i] << endl;

    }



    settlement = 1;

    Calendar KRWcalendar = SouthKorea();

    DCC = Actual365Fixed();

    Period coupon_period = Period(3, Months);

    BCC1 = ModifiedFollowing; // Fixed

    BCC2 = ModifiedFollowing; // Float



    boost::shared_ptr Index2 = boost::make_shared<IborIndex>(IborIndex("KRW91CD_onshore", coupon_period, settlement, KRWCurrency(), KRWcalendar, BCC2, true, DCC));

    vector <boost::shared_ptr <RateHelper>> helpers2;

    helpers2.push_back(boost::make_shared<DepositRateHelper>(Handle <Quote>(boost::make_shared<SimpleQuote>(KRW_Rates[0])), KRW_Tenor[0], settlement, KRWcalendar, BCC1, true, DCC));

    for (int i = 1; i < KRW_Tenor.size(); ++i) {

        helpers2.push_back(boost::make_shared<SwapRateHelper>(Handle <Quote>(boost::make_shared <SimpleQuote>(KRW_Rates[i])), KRW_Tenor[i], KRWcalendar, Quarterly, BCC1, DCC, Index2, Handle <Quote>(boost::make_shared<SimpleQuote>(0.0))));

    }

    std::cout << "Curve check " << endl;

    boost::shared_ptr <YieldTermStructure> KRWCD91_Onshore(new PiecewiseYieldCurve <Discount, Cubic, IterativeBootstrap>(Base, helpers2, DCC));

    KRWCD91_Onshore->enableExtrapolation();





    for (const auto& date : daily) {

        kr_spot.push_back(KRWCD91_Onshore->zeroRate(date, DCC, Continuous));

        kr_fwds.push_back(KRWCD91_Onshore->forwardRate(date, NullCalendar().advance(date, Period(1, Days)), DCC, Continuous));

    }

    kr_fwds.pop_back();



    /*

    vector <vector <double>> sampleM = { {1.0, 2.0},

             {3.0, 4.0},

             {5.0, 6.0} };



    std::cout << "Row size : " << sampleM.size() << "\t Column size : " << sampleM[0].size() << endl;



    vector <vector <double>> sampleM2;



    sampleM2.resize(sampleM.size(), vector<double>(sampleM[0].size(), 1.0));



    std::cout << "Row size : " << sampleM2.size() << "\t Column size : " << sampleM2[1].size() << endl;

    */



    // ************************** Simulation ************************

    int ns = 100000;

    int nt = daily.size() - 1;

    int jj = 0;

    int seed = 10125;



    daily_tenors.push_back(PayoffDCC.yearFraction(Base, daily[0]));

    vector <double> dt;

    for (int i = 1; i < nt + 1; ++i) {

        daily_tenors.push_back(PayoffDCC.yearFraction(Base, daily[i]));

        int ii = i - 1;

        dt.push_back(daily_tenors[i] - daily_tenors[ii]);

    }



    // random_device seed; mt19937 gen(seed()); // Random seed



    mt19937 gen(seed);

    normal_distribution<double> dist(0.0, 1.0);



    vector <vector <double>> corrmatrix = { { 1.0, cor },

                                                                               {cor, 1.0} };



    vector <vector<double>> L;

    vector <vector<double>> U;

    luDecomposition(corrmatrix, L, U);



    /*

    for (int i = 0; i < w.size(); ++i) {

             cout << z[i] << " ";

    }

    cout << endl;

    for (int i = 0; i < w.size(); ++i) {

             cout << w[i] << " ";

    }

    */





    std::cout << "\n Generate Underlying Path..." << endl;

    vector <vector<double>> FX(ns, vector<double>(nt + 1));

    vector <vector <double>> S(ns, vector<double>(nt + 1));

    vector <vector <double>> lnS(ns, vector<double>(nt + 1)); // R = ln(Perf/Perf)

    vector <vector <double>> barlnS(ns, vector<double>(nt + 1));

    vector <vector <double>> sigt(ns, vector<double>(nt + 1)); // Realized Vol

    vector <vector <double>> alpha(ns, vector<double>(nt + 1)); // Realized Participation

    vector <vector <double>> volidx(ns, vector<double>(nt + 1)); // Vol_Index

    vector <vector<double>> Payoff(ns, vector <double>(nt + 1));

    int VP = 30;



    vector <int> A = { 1,2,3,4,5,6,7,8,9,10 };



    /*

    //cout << total << endl;

    for (int i = 0; i < A.size(); ++i) {

             if (i <= A.size() -2) {

             barlnS[0].insert(barlnS[0].begin(), accumulate(A.end() - i - 2, A.end() - i, 0));

             }

             else {

                     barlnS[0].insert(barlnS[0].begin(), 0);

             }

    }

    for (int i = 0; i < A.size(); ++i) {

             cout << barlnS[0].at(i) << "\t" << barlnS[1].at(i) << endl;

    }

    */







    for (int i = 0; i < ns; ++i) {

        FX[i][0] = FX0; S[i][0] = S0; lnS[i][0] = 0.0;

        barlnS[i][0] = 0.0; sigt[i][0] = 0.0; alpha[i][0] = 0.273944728; volidx[i][0] = 96.48967361;

        for (int j = 0; j < nt; ++j) {

            vector <double> z = { dist(gen), dist(gen) };

            vector <double> w;



            w.resize(z.size());

            w = product_MatVec(L, z);



            jj = j + 1;

            S[i][jj] = S[i][j] * exp((kr_fwds[j] - 0.5 * sig * sig) * dt[j] + sig * sqrt(dt[j]) * w[0]);

            FX[i][jj] = FX[i][j] * exp((kr_fwds[j] - us_fwds[j] + 0.5 * fxsig * fxsig) * dt[j] + fxsig * sqrt(dt[j]) * w[1]); // Domestic currency : USD

            lnS[i][jj] = log(S[i][jj] / S[i][j]);



            vector <double> tempsigt(30);

            if (j >= 29) {

                barlnS[i][jj] = accumulate(lnS[i].begin() + j - 29, lnS[i].begin() + j + 1, (double)0.0) / 30;

                for (int k = 0; k < 30; ++k) { tempsigt[k] = pow(lnS[i][jj + (k - 29)] - barlnS[i][jj], 2); }

                sigt[i][jj] = sqrt(252 * accumulate(tempsigt.begin(), tempsigt.end(), (double)0.0) / 30);

                alpha[i][jj] = min(1.0, 0.05 / sigt[i][jj]);



            }

            else { barlnS[i][jj] = 0.0; sigt[i][jj] = 0.0; alpha[i][jj] = alpha[i][0]; volidx[i][jj] = volidx[i][0]; }

            volidx[i][jj] = volidx[i][j] * (1.0 + alpha[i][j] * (S[i][jj] / S[i][j] - 1.0));

        }

        Payoff[i][nt] = Notional * 1.0 * max(volidx[i][nt] / 100 - 1, 0.0) / FX[i][nt];

    }



    cout << "Checking ..." << endl;

    for (int j = 0; j < nt + 1; ++j) {

        cout << "j: " << j << "\t barlnS: " << barlnS[0][j] << "\t sigt: " << sigt[0][j] << "\t alpha: " << alpha[0][j] << "\t vol_Idx: " << volidx[0][j] << endl;

    }





    // cout << "Payoff Calcualting..." << endl;

    /*

    for (int i = 0; i < ns; ++i) {

             idx1 = afteridx1; idx2 = afteridx2;

             for (int j = 0; j < nt + 1; ++j) {

                     if (daily[j] == USDFixedDate[idx1]) {

                             USDflow[i][j] = USDNotional * USDFixedRate * USDFixedDCC.yearFraction(USDFixedDate[idx1 - 1], USDFixedDate[idx1]);

                             idx1++;

                             if (idx1 == USDFixedDate.size()) { idx1 = 0; }

                     }



                     if (daily[j] == KRWFixedDate[idx2]) {

                             KRWflow[i][j] = -KRWNotional * KRWFixedRate * KRWFixedDCC.yearFraction(KRWFixedDate[idx2 - 1], KRWFixedDate[idx2]) / S[i][j];

                             idx2++;

                             if (idx2 == KRWFixedDate.size()) { idx2 = 0; }

                     }



                     if (daily[j] == FinalExchangeDate) {

                             if (S[i][j] <= Cap) { Payoff[i][j] = -USDNotional * (K - Cap) / S[i][j]; cout << i << " " << j << ": " << "Case1" << endl; }

                             else if (S[i][j] > Cap && S[i][j] <= K) { Payoff[i][j] = -USDNotional * (K - S[i][j]) / S[i][j]; cout << i << " " << j << ": " << "Case2" << endl; }

                             else if (S[i][j] > K) { Payoff[i][j] = USDNotional * (S[i][j] - K) / S[i][j]; cout << i << " " << j << ": " << "Case3" << endl; }

                     }

             }

    }





    double KRWPrice = 0;

    double USDPrice = 0;

    */



    double Price = 0;



    cout << "\n Discounting..." << endl;

    for (int i = 0; i < ns; ++i) {

        for (int j = nt; j >= 1; --j) {

            jj = j - 1;

            Payoff[i][jj] += Payoff[i][j] / (1 + us_fwds[jj] * dt[jj]);

        }

        Price += Payoff[i][0];

        // cout << Payoff[i][0] << endl;

    }

    Price /= ns;





    cout << "\n Process is done. Option Value is " << Price << endl;





    // Outputs

    /*

    ofstream filename1(Outputpath + "Price.txt");

    if (filename1.is_open()) {

             filename1 << Price;

    }

    filename1.close();

    */



    ofstream filename2(Outputpath + "KRW Rates.txt");

    if (filename2.is_open()) {

        cout << "Write data: " << "Rates" << endl;

        filename2 << "Index" << "\t" << "Date" << "\t" << "Spot" << "\t" << "Fwd" << "\n";

        filename2 << "1" << "\t" << daily[0] << "\t" << kr_spot[0] << "\t" << "0.0" << "\n";

        for (int i = 1; i < kr_spot.size(); ++i) {

            filename2 << i + 1 << "\t" << daily[i] << "\t" << kr_spot[i] << "\t" << kr_fwds[i - 1] << "\n";

        }

    }

    filename2.close();





    ofstream filename3(Outputpath + "USD Rates.txt");

    if (filename3.is_open()) {

        cout << "Write data: " << "Rates" << endl;

        filename3 << "Index" << "\t" << "Date" << "\t" << "Spot" << "\t" << "Fwd" << "\n";

        filename3 << "1" << "\t" << daily[0] << "\t" << us_spot[0] << "\t" << "0.0" << "\n";

        for (int i = 1; i < us_spot.size(); ++i) {

            filename3 << i + 1 << "\t" << daily[i] << "\t" << us_spot[i] << "\t" << us_fwds[i - 1] << "\n";

        }

    }

    filename3.close();

    /*



    ofstream filename4(Outputpath + "Underlying path.txt");

    if (filename4.is_open()) {

             cout << "Write data: " << "Underlying path" << endl;

             for (const auto& t: daily_tenors){

                     filename4 << "\t" << t ;

             }

             filename4 << "\n";

             for (int i = 0; i < ns; ++i) {

                     filename4 << i + 1;

                     for (int j = 0; j < nt+1; ++j) {

                     filename4 << "\t" << S[i][j];

                     }

                     filename4 << "\n";

             }

    }

    filename4.close();

    */

    /*

    ofstream filename5(Outputpath + "KRW Flow.txt");

    if (filename5.is_open()) {

             cout << "Write data: " << "Flow" << endl;

             for (const auto& t : daily_tenors) {

                     filename5 << "\t" << t;

             }

             for (int i = 0; i < ns; ++i) {

                     filename5 << i + 1;

                     for (int j = 0; j < nt + 1; ++j) {

                             filename5 << "\t" << KRWflow[i][j];

                     }

                     filename5 << "\n";

             }

    }

    filename5.close();



    ofstream filename6(Outputpath + "USD Flow.txt");

    if (filename6.is_open()) {

             cout << "Write data: " << "Flow" << endl;

             for (const auto& t : daily_tenors) {

                     filename6 << "\t" << t;

             }

             for (int i = 0; i < ns; ++i) {

                     filename6 << i + 1;

                     for (int j = 0; j < nt + 1; ++j) {

                             filename6 << "\t" << USDflow[i][j];

                     }

                     filename6 << "\n";

             }

    }

    filename6.close();





    ofstream filename7(Outputpath + "Payoff.txt");

    if (filename7.is_open()) {

             cout << "Write data: " << "Payoff" << endl;

             for (const auto& t : daily_tenors) {

                     filename7 << "\t" << t;

             }

             for (int i = 0; i < ns; ++i) {

                     filename7 << i + 1;

                     for (int j = 0; j < nt + 1; ++j) {

                             filename7 << "\t" << Payoff[i][j];

                     }

                     filename7 << "\n";

             }

    }

    filename7.close();





    */



    return 0;

}



void SchedulePrint(Schedule S) {



    int i = 0;

    /* string str = GET_NAME(S);

    cout << str << endl;

    */



    for (i = 0; i < S.size(); ++i) {

        cout << S[i] << endl;

    };

    cout << "\n";



}



void SchedulePrint(vector <Date> S) {



    int i = 0;

    /* string str = GET_NAME(S);

    cout << str << endl;

    */



    for (i = 0; i < S.size(); ++i) {

        cout << S[i] << endl;

    };

    cout << "\n";



}



void luDecomposition(const std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& L, std::vector<std::vector<double>>& U) {

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



vector <string> Open_getline(string filename) {

    /*

    filename은 확장자까지 포함해야합니다.

    output은 vector <string> 형식입니다. 각 목적에 맞게 데이터를 변형하세요.

    */

    vector <string> output;



    if (fs::exists(Inputpath + filename)) {

        ifstream localfilename(Inputpath + filename);

        if (localfilename.is_open()) {

            cout << "Open Inputs: " << filename << endl;

            string line;



            while (getline(localfilename, line)) {

                output.push_back(line);

            }

        }

        localfilename.close();

    }



    return output;

}



vector <string> Open_getline(string filename, bool localPrintSwitch) {

    /*

    filename은 확장자까지 포함해야합니다.

    output은 vector <string> 형식입니다. 각 목적에 맞게 데이터를 변형하세요.

    */

    vector <string> output;



    if (fs::exists(Inputpath + filename)) {

        ifstream localfilename(Inputpath + filename);

        if (localfilename.is_open()) {

            cout << "Open Inputs: " << filename << endl;

            string line;



            if (localPrintSwitch == true) {

                int i = 0;

                while (getline(localfilename, line)) {

                    output.push_back(line);

                    cout << output[i] << endl;

                    i += 1;

                }

            }

            else {

                while (getline(localfilename, line)) {

                    output.push_back(line);

                }



            }

        }

        localfilename.close();

    }

    return output;

}



vector <double> product_MatVec(const vector <vector <double>>& A, vector <double>& v) {



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