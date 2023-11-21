#include "Bootstrapping.h"



#include <iostream>

#include <ql/quantlib.hpp>

#include <stdio.h>

#include <filesystem>

#include <sstream>

#include <ql/utilities/dataparsers.hpp>

#include <fstream>



using namespace std;

using namespace QuantLib;

namespace fs = std::filesystem;





// Bootsrappers

void Bootstraps::ReadTxt_bloomcurve(vector <Period>& tenor, vector <Rate>& rate, string filepath, vector <int> ReadIndex) {

    //Data(ReadIndex)[Type]

    //Tenor(1)[int] STD(2)[string] Market(3)[rate] Shift(4)[rate] Shifted(5)[rate] Zero(6)[rate] Discount(7)[rate]





    ifstream file(filepath);

    vector <int> column_idx;



    for (const auto& i : ReadIndex) {

        column_idx.push_back(i);

    }

    vector <int> day;

    vector <string> std;        vector <Rate> market;     vector <Rate> shift;

    vector <Rate> Zero;                   vector <Rate> Discount;





    if (file.is_open()) {



        int flag = 0;

        string line;

        while (std::getline(file, line)) {

            std::istringstream iss(line);

            std::string column;

            std::vector<std::string> columns;



            while (std::getline(iss, column, '\t')) {

                columns.push_back(column);

            }



            if (flag > 0) {

                for (int i = 0; i < column_idx.size(); ++i) {

                    if (column_idx[i] <= columns.size()) {

                        int idx = column_idx[i] - 1;

                        if (i == 0) { int value = stod(columns[idx]); day.push_back(value); }

                        else if (i == 1) { std.push_back(columns[idx]); }

                        else if (i == 2) { Rate value = stod(columns[idx]); value *= 0.01; market.push_back(value); }

                        else if (i == 3) { Rate value = stod(columns[idx]); value *= 0.01; shift.push_back(value); }

                        else if (i == 4) { Rate value = stod(columns[idx]); value *= 0.01; rate.push_back(value); }

                        else if (i == 5) { Rate value = stod(columns[idx]); value *= 0.01; Zero.push_back(value); }

                        else if (i == 6) { Rate value = stod(columns[idx]); value *= 0.01; Discount.push_back(value); }

                    }

                    else { std::cout << "Invalid column number." << std::endl; break; }

                }

            }

            flag++;

        }

        file.close();



        for (int i = 0; i < day.size(); ++i) {

            vector <string> STD;

            STD.push_back("WK"); STD.push_back("MO"); STD.push_back("YR");

            if (std[i].compare(STD[0]) == 0) { tenor.push_back(Period(day[i], Weeks)); }

            else if (std[i].compare(STD[1]) == 0) { tenor.push_back(Period(day[i], Months)); }

            else if (std[i].compare(STD[2]) == 0) { tenor.push_back(Period(day[i], Years)); }

            // cout << KRW_Tenor[i] << endl;

        }

    }

    else {

        throw "Failed to open the file. | ReadTxt_bloomcurve";

        // Return from the program indicating an error.

    }

}

void Bootstraps::savetoTxt(string filename, vector <Date> daily, vector <double> tenors, vector <Rate> spot, vector <Rate> fwds) {

    // needs all outputs

    ofstream file(filename);

    if (file.is_open()) {

        std::cout << "Write data: " << "Rates" << endl;

        file << "Index" << "\t" << "Date" << "\t" << "Tenors" << "\t" << "Spot" << "\t" << "Fwd" << "\n";

        for (int i = 0; i < spot.size(); ++i) {

            file << i + 1 << "\t" << daily[i] << "\t" << tenors[i] << "\t" << spot[i] << "\t" << fwds[i] << "\n";

        }

    }

    file.close();

}



void Bootstraps::setBase(Date d) { base = d; }

void Bootstraps::setTerm(Date d) { term = d; }

void Bootstraps::setCurve(boost::shared_ptr <YieldTermStructure> c) { curve = c; }

void Bootstraps::setDayCounter(DayCounter d) { dc = d; }

void Bootstraps::setCalendar(Calendar c) { calendar = c; }

void Bootstraps::setCompounding(Compounding c) { compound = c; }

void Bootstraps::setFreq(Frequency f) { freq = f; }

void Bootstraps::setBCC(BusinessDayConvention b) { bcc = b; }



Date Bootstraps::getBase() { return base; }

Date Bootstraps::getTerm() { return term; }

DayCounter Bootstraps::getDayCounter() { return dc; }

Calendar Bootstraps::getCalendar() { return calendar; }

Compounding Bootstraps::getCompounding() { return compound; }

Frequency Bootstraps::getFreq() { return freq; }

BusinessDayConvention Bootstraps::getBCC() { return bcc; }



vector <Date> Bootstraps::getdaily() {

    Calendar calendar = getCalendar();

    BusinessDayConvention bcc = getBCC();



    vector <Date> days;



    int i = base.serialNumber();

    if (calendar.adjust(Date(i), bcc) == term) { days.push_back(calendar.adjust(Date(i), bcc)); }

    else { days.push_back(calendar.adjust(Date(i), bcc)); }

    for (int i = base.serialNumber() + 1; i <= term.serialNumber(); ++i) {

        if (calendar.adjust(Date(i), bcc) == term) {

            if (calendar.adjust(Date(i), bcc) == calendar.adjust(Date(i - 1), bcc)) { ; }

            else { days.push_back(calendar.adjust(Date(i), bcc)); break; }

        }

        else {

            if (calendar.adjust(Date(i), bcc) == calendar.adjust(Date(i - 1), bcc)) { ; }

            else { days.push_back(calendar.adjust(Date(i), bcc)); }

        }

    }

    return days;

}



vector <Rate> Bootstraps::getspots() {

    DayCounter dcc = getDayCounter();

    Calendar calendar = getCalendar();

    Compounding compound = getCompounding();

    Frequency freq = getFreq();

    BusinessDayConvention bcc = getBCC();

    vector <Date> days = getdaily();



    vector <Rate> spot;

    for (const auto& date : days) {

        spot.push_back(curve->zeroRate(date, dcc, compound, freq));

    }

    return spot;

}



vector <Rate> Bootstraps::getforwards() {

    DayCounter dcc = getDayCounter();

    Calendar calendar = getCalendar();

    Compounding compound = getCompounding();

    Frequency freq = getFreq();

    BusinessDayConvention bcc = getBCC();



    vector <Rate> fwd;

    vector <Date> days = getdaily();



    for (const auto& date : days) {

        fwd.push_back(curve->forwardRate(date, calendar.advance(date, Period(1, Days)), dcc, compound, freq));

    }

    fwd.pop_back();

    fwd.insert(fwd.begin(), (Rate)0.0);



    return fwd;

}



Rate Bootstraps::getforward(Date d1, Date d2) {

    DayCounter dcc = getDayCounter();

    Calendar calendar = getCalendar();

    Compounding compound = getCompounding();

    Frequency freq = getFreq();

    BusinessDayConvention bcc = getBCC();



    Rate fwd = curve->forwardRate(d1, d2, dcc, compound, freq);



    return fwd;



}



vector <double> Bootstraps::gettenors() {

    DayCounter dcc = getDayCounter();

    Calendar calendar = getCalendar();

    Compounding compound = getCompounding();

    Frequency freq = getFreq();

    BusinessDayConvention bcc = getBCC();



    vector <double> daily_tenors;

    vector <Date> days = getdaily();



    for (int i = 0; i < days.size() - 1; ++i) {

        daily_tenors.push_back(dcc.yearFraction(base, days[i + 1]));

    }

    daily_tenors.insert(daily_tenors.begin(), (double)0.0);



    return daily_tenors;

}



boost::shared_ptr <YieldTermStructure> Bootstraps::getcurve() { return curve; }





// KRWCD91 curve

BootstrapKRWCD91::BootstrapKRWCD91(Date Base, Date Term, string filename, Bootstraps::Interpolationtype interp) {

    file = filename;

    setBase(Base); setTerm(Term);

    base = getBase();

    term = getTerm();



    ReadTxt_bloomcurve(tenor, rate, file, Readindex);

    boost::shared_ptr <IborIndex> Index = boost::make_shared<IborIndex>(IborIndex("KRWCD91 Oneshore", coupon_period, settlement, KRWCurrency(), innercalendar, bcc2, true, innerdcc));

    vector <boost::shared_ptr <RateHelper>> helpers;

    helpers.push_back(boost::make_shared<DepositRateHelper>(Handle <Quote>(boost::make_shared<SimpleQuote>(rate[0])), tenor[0], settlement, innercalendar, bcc1, true, innerdcc)); // settlement for first rate is 1?? (for Deposite rate)

    for (int i = 1; i < tenor.size(); ++i) {

        helpers.push_back(boost::make_shared<SwapRateHelper>(Handle <Quote>(boost::make_shared<SimpleQuote>(rate[i])), tenor[i], innercalendar, Quarterly, bcc1, innerdcc, Index, Handle <Quote>(boost::make_shared<SimpleQuote>(0.0))));

    }



    switch (interp) {

    case linear: {

        curve = boost::shared_ptr <YieldTermStructure>(new PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap>(base, helpers, innerdcc));

        break; }

    case cubic: {

        // Cubic logic이 현재 좀 이상함 이건 Quantlib 모듈 문제인거 같은데

        boost::shared_ptr <YieldTermStructure> curve(new PiecewiseYieldCurve<ZeroYield, Cubic, IterativeBootstrap>(base, helpers, innerdcc));

        break;  }

    default: {

        curve = boost::shared_ptr <YieldTermStructure>(new PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap>(base, helpers, innerdcc));

        break; }

    }

    curve->enableExtrapolation();



    setCurve(curve);

}



void BootstrapKRWCD91::setting(DayCounter dc, Calendar cal, Compounding com, Frequency f, BusinessDayConvention bc) {

    setDayCounter(dc);         setCalendar(cal); setCompounding(com); setFreq(f); setBCC(bc);

}



vector <Date> BootstrapKRWCD91::daily() { return getdaily(); }



vector <Rate> BootstrapKRWCD91::spots() { return getspots(); }



vector <Rate> BootstrapKRWCD91::forwards() { return getforwards(); }



vector <double> BootstrapKRWCD91::tenors() { return gettenors(); }



Rate BootstrapKRWCD91::forward(Date d1, Date d2) { return getforward(d1, d2); }



boost::shared_ptr <YieldTermStructure> BootstrapKRWCD91::curveout() { return getcurve(); }



void BootstrapKRWCD91::saveTxt(string filename, vector <Date> daily, vector <double> tenors, vector <Rate> spot, vector <Rate> fwds) {

    savetoTxt(filename, daily, tenors, spot, fwds);

}



BootstrapKRWCD91::~BootstrapKRWCD91() {

    std::cout << "Using Curve : KRW CD91 OnShore" << endl;

}







// USDSofr Curve

BootstrapUSDSofr::BootstrapUSDSofr(Date Base, Date Term, string filename, Bootstraps::Interpolationtype interp) {

    file = filename;

    setBase(Base); setTerm(Term);

    base = getBase();

    term = getTerm();



    ReadTxt_bloomcurve(tenor, rate, file, Readindex);



    boost::shared_ptr Index = boost::make_shared<OvernightIndex>(OvernightIndex("USD Sofr", settlement, USDCurrency(), innercalendar, innerdcc));

    vector <boost::shared_ptr <RateHelper>> helpers;

    for (int i = 0; i < tenor.size(); ++i) {

        helpers.push_back(boost::make_shared<OISRateHelper>(settlement,

            tenor[i],

            Handle <Quote>(boost::make_shared <SimpleQuote>(rate[i])),

            Index,

            Handle<YieldTermStructure> {},

            false,

            pay_delay,

            bcc2,

            Annual,

            innercalendar));

    }



    switch (interp) {

    case linear: {

        curve = boost::shared_ptr <YieldTermStructure>(new PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap>(base, helpers, innerdcc));

        break; }

    case cubic: {

        boost::shared_ptr <YieldTermStructure> curve(new PiecewiseYieldCurve<ZeroYield, Cubic, IterativeBootstrap>(base, helpers, innerdcc));

        break;  }

    default: {

        curve = boost::shared_ptr <YieldTermStructure>(new PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap>(base, helpers, innerdcc));

        break; }

    }

    curve->enableExtrapolation();



    setCurve(curve);

}



void BootstrapUSDSofr::setting(DayCounter dc, Calendar cal, Compounding com, Frequency f, BusinessDayConvention bc) {

    setDayCounter(dc);         setCalendar(cal); setCompounding(com); setFreq(f); setBCC(bc);

}



vector <Date> BootstrapUSDSofr::daily() { return getdaily(); }



vector <Rate> BootstrapUSDSofr::spots() { return getspots(); }



vector <Rate> BootstrapUSDSofr::forwards() { return getforwards(); }



vector <double> BootstrapUSDSofr::tenors() { return gettenors(); }



Rate BootstrapUSDSofr::forward(Date d1, Date d2) { return getforward(d1, d2); }



void BootstrapUSDSofr::saveTxt(string filename, vector <Date> daily, vector <double> tenors, vector <Rate> spot, vector <Rate> fwds) {

    savetoTxt(filename, daily, tenors, spot, fwds);

}



BootstrapUSDSofr::~BootstrapUSDSofr() {

    std::cout << "Using Curve : USD Sofr (vs. Fixed)" << endl;

}

// HKDOIS Curve

BootstrapHKDOIS::BootstrapHKDOIS(Date Base, Date Term, string filename, Bootstraps::Interpolationtype interp) {

    file = filename;

    setBase(Base); setTerm(Term);

    base = getBase();

    term = getTerm();



    ReadTxt_bloomcurve(tenor, rate, file, Readindex);



    boost::shared_ptr Index = boost::make_shared<OvernightIndex>(OvernightIndex("HKD OIS", settlement, HKDCurrency(), innercalendar, innerdcc));

    vector <boost::shared_ptr <RateHelper>> helpers;

    for (int i = 0; i < tenor.size(); ++i) {

        helpers.push_back(boost::make_shared<OISRateHelper>(settlement,

            tenor[i],

            Handle <Quote>(boost::make_shared <SimpleQuote>(rate[i])),

            Index,

            Handle<YieldTermStructure> {},

            false,

            pay_delay,

            bcc2,

            Quarterly,

            innercalendar));

    }



    switch (interp) {

    case linear: {

        curve = boost::shared_ptr <YieldTermStructure>(new PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap>(base, helpers, innerdcc));

        break; }

    case cubic: {

        boost::shared_ptr <YieldTermStructure> curve(new PiecewiseYieldCurve<ZeroYield, Cubic, IterativeBootstrap>(base, helpers, innerdcc));

        break;  }

    default: {

        curve = boost::shared_ptr <YieldTermStructure>(new PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap>(base, helpers, innerdcc));

        break; }

    }

    curve->enableExtrapolation();



    setCurve(curve);

}



void BootstrapHKDOIS::setting(DayCounter dc, Calendar cal, Compounding com, Frequency f, BusinessDayConvention bc) {

    setDayCounter(dc);         setCalendar(cal); setCompounding(com); setFreq(f); setBCC(bc);

}



vector <Date> BootstrapHKDOIS::daily() { return getdaily(); }



vector <Rate> BootstrapHKDOIS::spots() { return getspots(); }



vector <Rate> BootstrapHKDOIS::forwards() { return getforwards(); }



vector <double> BootstrapHKDOIS::tenors() { return gettenors(); }



Rate BootstrapHKDOIS::forward(Date d1, Date d2) { return getforward(d1, d2); }



void BootstrapHKDOIS::saveTxt(string filename, vector <Date> daily, vector <double> tenors, vector <Rate> spot, vector <Rate> fwds) {

    savetoTxt(filename, daily, tenors, spot, fwds);

}



BootstrapHKDOIS::~BootstrapHKDOIS() {

    std::cout << "Using Curve : HKD OIS" << endl;

}



// ESTR Curve

BootstrapESTR::BootstrapESTR(Date Base, Date Term, string filename, Bootstraps::Interpolationtype interp) {

    file = filename;

    setBase(Base); setTerm(Term);

    base = getBase();

    term = getTerm();



    ReadTxt_bloomcurve(tenor, rate, file, Readindex);



    boost::shared_ptr Index = boost::make_shared<OvernightIndex>(OvernightIndex("ESTR OIS", settlement, EURCurrency(), innercalendar, innerdcc));

    vector <boost::shared_ptr <RateHelper>> helpers;

    for (int i = 0; i < tenor.size(); ++i) {

        helpers.push_back(boost::make_shared<OISRateHelper>(settlement,

            tenor[i],

            Handle <Quote>(boost::make_shared <SimpleQuote>(rate[i])),

            Index,

            Handle<YieldTermStructure> {},

            false,

            pay_delay,

            bcc2,

            Annual,

            innercalendar));

    }



    switch (interp) {

    case linear: {

        curve = boost::shared_ptr <YieldTermStructure>(new PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap>(base, helpers, innerdcc));

        break; }

    case cubic: {

        boost::shared_ptr <YieldTermStructure> curve(new PiecewiseYieldCurve<ZeroYield, Cubic, IterativeBootstrap>(base, helpers, innerdcc));

        break;  }

    default: {

        curve = boost::shared_ptr <YieldTermStructure>(new PiecewiseYieldCurve<ZeroYield, Linear, IterativeBootstrap>(base, helpers, innerdcc));

        break; }

    }

    curve->enableExtrapolation();



    setCurve(curve);

}



void BootstrapESTR::setting(DayCounter dc, Calendar cal, Compounding com, Frequency f, BusinessDayConvention bc) {

    setDayCounter(dc);         setCalendar(cal); setCompounding(com); setFreq(f); setBCC(bc);

}



vector <Date> BootstrapESTR::daily() { return getdaily(); }



vector <Rate> BootstrapESTR::spots() { return getspots(); }



vector <Rate> BootstrapESTR::forwards() { return getforwards(); }



vector <double> BootstrapESTR::tenors() { return gettenors(); }



Rate BootstrapESTR::forward(Date d1, Date d2) { return getforward(d1, d2); }



void BootstrapESTR::saveTxt(string filename, vector <Date> daily, vector <double> tenors, vector <Rate> spot, vector <Rate> fwds) {

    savetoTxt(filename, daily, tenors, spot, fwds);

}



BootstrapESTR::~BootstrapESTR() {

    std::cout << "Using Curve : ESTR OIS" << endl;