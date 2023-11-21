#include <iostream>

#include <stdio.h>

#include <ql/quantlib.hpp>

#include <ql/utilities/dataparsers.hpp >



#include <fstream>

#include <sstream>

#include <filesystem>



#include <numeric>

#include <algorithm>



#include "Bootstrapping.h"

#include "Dataprogress.h"



using namespace QuantLib;

using namespace std;

namespace fs = std::filesystem;



int main() {

    Dataprogress datatreat;



    Date base = Date(31, Mar, 2023);

    Settings::instance().evaluationDate() = base;



    enum Parties {

        PartyA = -1,

        PartyB = 1,

        Nomura = -1,

        Hana = 1

    };



    Date trade = Date(24, Feb, 2021);

    Date effect = Date(3, Mar, 2021);

    Date term = Date(29, Feb, 2024); // or early Termination Date

    Date value_date = Date(23, Feb, 2024);







    Calendar calendar = JointCalendar(vector <Calendar> {SouthKorea(), UnitedStates(UnitedStates::NYSE)}, JoinHolidays);

    Calendar calendar_CCY = SouthKorea();



    vector <Date> knockout_date{ Date(24, Aug, 2021), Date(24,Feb, 2022),Date(24, Aug, 2022), Date(24, Feb, 2023), Date(24, Aug, 2023), value_date };

    vector <Date> fxfixing_date;

    for (const auto& d : knockout_date) {

        fxfixing_date.push_back(calendar_CCY.advance(d, Period(2, Days), Following));

    }

    Date cashsettle = calendar_CCY.advance(fxfixing_date.back(), Period(2, Days), Following);

    vector <Date> earlyterm{ Date(24, Aug, 2021), Date(24,Feb, 2022),Date(24, Aug, 2022), Date(24, Feb, 2023), Date(24, Aug, 2023), cashsettle };



    vector <double> trigger_rate{ 0.95, 0.90, 0.85, 0.85, 0.80, 0.80 };



    vector <double> ko_per{ 1.0, 2.0, 3.0, 4.0, 5.0 ,6.0 };





    while (1)

    {

        if (earlyterm[0] >= base) { break; }

        knockout_date.erase(knockout_date.begin());

        earlyterm.erase(earlyterm.begin());

        fxfixing_date.erase(fxfixing_date.begin());

        trigger_rate.erase(trigger_rate.begin());

        ko_per.erase(ko_per.begin());

    }

    int exernum = trigger_rate.size();



    Schedule float_schedule = MakeSchedule().from(effect)

        .to(Date(29, Feb, 2024))

        .withRule(DateGeneration::Forward)

        .withCalendar(calendar)

        .withConvention(ModifiedFollowing)

        .withTenor(Period(3, Months));



    vector <Date> float_dates;

    for (int i = 0; i < float_schedule.size(); ++i) {

        float_dates.push_back(float_schedule[i]);

    }



    while (1)

    {

        if (float_dates[1] >= base) { break; }

        float_dates.erase(float_dates.begin());

    }

    datatreat.SchedulePrint(float_dates);



    double notional = 10000000000;

    DayCounter float_dc = Actual365Fixed();

    double spread = 0.15 / 100;

    BootstrapKRWCD91 curveKRW(base, term, "./Inputs/KRWCD91.txt", Bootstraps::linear);

    curveKRW.setting(float_dc, calendar);



    vector <double> float_flow; // Reset Date 등 세부 디테일은 없으므로 주의 (fixing rate 확인할것)

    double fixing_rate = 3.76 / 100;

    for (int i = 2; i < float_dates.size(); ++i) {

        int ii = i - 1;

        float_flow.push_back(PartyB * notional * float_dc.yearFraction(float_dates[ii], float_dates[i]) * (curveKRW.forward(float_dates[ii], float_dates[i]) + spread));

    }

    float_flow.insert(float_flow.begin(), (double)PartyB * notional * float_dc.yearFraction(float_dates[0], float_dates[1]) * fixing_rate);

    for (int i = 0; i < float_flow.size(); ++i)

    {

        std::cout << float_flow[i] << endl;

    }



    // SPX(0), HSCEI(1), SX5E(2)

    double initial_SPX = 3925.43; // 0

    double initial_HSCEI = 11509.73; // 1

    double initial_SX5E = 3705.99; // 2

    vector <double> initials{ initial_SPX , initial_HSCEI, initial_SX5E };



    double knockout_SPX = 1962.715;

    double knockout_HSCEI = 5754.865;

    double knockout_SX5E = 1852.995;

    vector <double> knockouts{ knockout_SPX , knockout_HSCEI , knockout_SX5E };



    int basketnum = initials.size();



    vector <vector<double>> triggers;

    triggers.resize(basketnum, vector(exernum, 0.0));



    for (int i = 0; i < basketnum; ++i) {

        for (int j = 0; j < exernum; ++j) { triggers[i][j] = (trigger_rate[j] * initials[i]); }

    }



    DayCounter dc = Actual365Fixed();

    BootstrapUSDSofr curveUSD(base, term, "./Inputs/USDSofr.txt");

    curveUSD.setting(dc, calendar);

    BootstrapHKDOIS curveHKD(base, term, "./Inputs/HKDOis.txt");

    curveHKD.setting(dc, calendar);

    BootstrapESTR curveEUR(base, term, "./Inputs/EUROis.txt");

    curveHKD.setting(dc, calendar);



    vector <Date> daily_dates = curveUSD.daily();

    vector <double> tenors = curveUSD.tenors();



    vector <Rate> krw_spots = curveUSD.spots();

    vector <Rate> krw_fwds = curveUSD.forwards();

    vector <Rate> usd_spots = curveUSD.spots();

    vector <Rate> usd_fwds = curveUSD.forwards();

    vector <Rate> hkd_spots = curveUSD.spots();

    vector <Rate> hkd_fwds = curveUSD.forwards();

    vector <Rate> eur_spots = curveUSD.spots();

    vector <Rate> eur_fwds = curveUSD.forwards();



    curveKRW.saveTxt("./Outputs/KRW Curve.txt", daily_dates, tenors, krw_spots, krw_fwds);

    curveUSD.saveTxt("./Outputs/USD Curve.txt", daily_dates, tenors, usd_spots, usd_fwds);

    curveHKD.saveTxt("./Outputs/HKD Curve.txt", daily_dates, tenors, hkd_spots, hkd_fwds);

    curveEUR.saveTxt("./Outputs/EUR Curve.txt", daily_dates, tenors, eur_spots, eur_fwds);







    // Domestic : KRW, Measure : Domestic

    double s0_SPX = 4109.31;

    double s0_HSCEI = 6968.86;

    double s0_SX5E = 4315.05;

    vector <double> s0{ s0_SPX, s0_HSCEI, s0_SX5E };



    double q_SPX = 0.01 * (1.57 + 1.27 + 1.60) / 3;

    double q_HSCEI = 0.01 * (2.81 + 2.82 + 3.77) / 3;

    double q_SX5E = 0.01 * (3.53 + 2.38 + 2.18) / 3;

    vector <double> q{ q_SPX, q_HSCEI, q_SX5E };



    double sig_SPX = 0.23306426064253;

    double sig_HSCEI = 0.334137084868717;

    double sig_SX5E = 0.200183489979534;



    // sig 조정 (변동성 range 조정 후 범위안에 값이 들어오는 지 확인

    sig_SPX *= 3;

    sig_HSCEI *= 3;

    sig_SX5E *= 3;



    vector <double> sig{ sig_SPX, sig_HSCEI, sig_SX5E };



    double fxsig_USD = 0.116792399209799;

    double fxsig_HKD = 0.116375036528614;

    double fxsig_EUR = 0.0897414831073909;





    vector <double> fxsig{ fxsig_USD ,fxsig_HKD ,fxsig_EUR };



    // USD, HKD, EUR

    vector <double> fxcor{ -0.093, -0.566, 0.162 };



    vector <double> vols;

    for (int i = 0; i < basketnum; ++i) {

        vols.push_back(sqrt(pow(sig[i], 2) + pow(fxsig[i], 2) + 2.0 * fxcor[i] * sig[i] * fxsig[i]));

    }







    // MCS

    int ns = 10000;

    int nt = daily_dates.size() - 1;

    int seed = 10125;



    vector <double> dt;

    for (int i = 1; i < nt + 1; ++i) {

        int ii = i - 1;

        dt.push_back(tenors[i] - tenors[ii]);

    }



    mt19937 gen(seed);

    normal_distribution<double> dist(0.0, 1.0);



    vector <vector <double>> corrmatrix = { { 1.0, 0.484042851727314, 0.402730242519576 },

                                                                               {0.484042851727314, 1.0, 0.200183489979534},

                                                                               {0.402730242519576, 0.200183489979534, 1.0} };



    vector <vector<double>> L;

    vector <vector<double>> U;

    Dataprogress().luDecomposition(corrmatrix, L, U);



    vector <vector <double>> s1(ns, vector<double>(nt + 1)); // SPX

    vector <vector <double>> s2(ns, vector<double>(nt + 1)); // HSCEI

    vector <vector <double>> s3(ns, vector<double>(nt + 1)); // SX5E

    vector <vector<double>> payoff(ns, vector <double>(nt + 1));





    vector <double> knockin;

    for (int i = 0; i < basketnum; ++i) { knockin.push_back(initials[i] * 0.5); }

    vector <int> knockin_flag(ns, 0);



    int value_idx = distance(daily_dates.begin(), find(daily_dates.begin(), daily_dates.end(), value_date)); //distance는 거래를 계산하는 method (daily_dates.begin()과 index사이의 차이)

    vector <int> knockout_idx;

    for (const auto& d : knockout_date) {

        knockout_idx.push_back(distance(daily_dates.begin(), find(daily_dates.begin(), daily_dates.end(), d)));

    }

    vector<int> knockout_flag(ns, 0);



    std::cout << daily_dates[knockout_idx[0]] << " " << daily_dates[knockout_idx[1]] << endl;



    std::cout << "\n Generate Underlying Path..." << endl;

    for (int i = 0; i < ns; ++i) {

        s1[i][0] = s0[0]; s2[i][0] = s0[1]; s3[i][0] = s0[1];

        int flag = 0;

        for (int j = 0; j < nt; ++j) {

            vector <double> z = { dist(gen), dist(gen), dist(gen) };

            vector <double> w;



            w.resize(z.size());

            w = Dataprogress().product_MatVec(L, z);



            int jj = j + 1;

            s1[i][jj] = s1[i][j] * exp((krw_fwds[jj] - q[0] - 0.5 * pow(vols[0], 2)) * dt[j] + vols[0] * sqrt(dt[j]) * w[0]);

            s2[i][jj] = s2[i][j] * exp((krw_fwds[jj] - q[1] - 0.5 * pow(vols[1], 2)) * dt[j] + vols[1] * sqrt(dt[j]) * w[1]);

            s3[i][jj] = s3[i][j] * exp((krw_fwds[jj] - q[2] - 0.5 * pow(vols[2], 2)) * dt[j] + vols[2] * sqrt(dt[j]) * w[2]);



            if (flag < exernum) {

                if (daily_dates[knockout_idx[flag]] == daily_dates[j]) {

                    if (knockout_flag[i] == 0) { // 기존에 knockout이 발생하지 않은 경우에만 확인

                        if (s1[i][jj] >= knockouts[0] && s2[i][jj] >= knockouts[1] && s3[i][jj] >= knockouts[2]) {

                            payoff[i][jj] = PartyA * notional * ko_per[flag] / 100 * 0.0371 + PartyB * 742207439;

                            knockout_flag[i] = jj; // flag는 일자까지 알려줍니다.

                        }

                    }

                    flag += 1;

                }

            }

            if (daily_dates[value_idx] == daily_dates[j]) {

                if (knockout_flag[i] == 0) {

                    if (s1[i][jj] < knockin[0] || s2[i][jj] < knockin[1] || s3[i][jj] < knockin[2]) {

                        double worst = min(min(s1[i][jj] / initials[0], s2[i][jj] / initials[1]), s3[i][jj] / initials[2]);

                        payoff[i][jj] = PartyB * max(0.0, notional * (1 - worst));

                        knockin_flag[i] = jj;

                    }

                    else {

                        payoff[i][jj] = PartyA * notional * 0.2226;

                    }

                }

            }



        }

    }





    std::cout << "Checking ..." << endl;

    for (int i = 0; i < ns; ++i) {

        std::cout << "knock out flag " << knockout_flag[i] << "knock in flag" << knockin_flag[i] << endl;



    }

    double price = 0;



    std::cout << "\n Discounting..." << endl;

    for (int i = 0; i < ns; ++i) {

        for (int j = nt; j >= 1; --j) {

            int jj = j - 1;

            payoff[i][jj] += payoff[i][j] / (1 + krw_fwds[j] * dt[jj]);

        }

        price += payoff[i][0];

        // cout << Payoff[i][0] << endl;

    }

    price /= ns;





    std::cout << "\n Process is done. Option Value is " << price << endl;





    return 0;

}