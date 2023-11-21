#pragma once
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



class Bootstraps {

private:

    const vector <int> Readindex = { 1,2,3,4,5,6,7 };

    DayCounter dc = Thirty360(Thirty360::BondBasis);

    Calendar calendar = NullCalendar();

    Compounding compound = Continuous;

    Frequency freq = Annual;

    BusinessDayConvention bcc = ModifiedFollowing;



    Date base;

    Date term;

    boost::shared_ptr <YieldTermStructure> curve;



protected:

    void ReadTxt_bloomcurve(vector <Period>&, vector <Rate>&, string, vector <int>);

    void savetoTxt(string, vector <Date>, vector <double>, vector <Rate>, vector <Rate>);

    void setBase(Date);

    void setTerm(Date);

    void setCurve(boost::shared_ptr <YieldTermStructure>);



    Date getBase();

    Date getTerm();



    vector <Date> getdaily();

    vector <Rate> getspots();

    vector <Rate> getforwards();

    vector <double> gettenors();

    Rate getforward(Date, Date);



public:

    enum Interpolationtype { linear, cubic };



    // setting value

    void setDayCounter(DayCounter);

    void setCalendar(Calendar);

    void setCompounding(Compounding);

    void setFreq(Frequency);

    void setBCC(BusinessDayConvention);



    // get value

    DayCounter getDayCounter();

    Calendar getCalendar();

    Compounding getCompounding();

    Frequency getFreq();

    BusinessDayConvention getBCC();

    boost::shared_ptr <YieldTermStructure> getcurve();

};



class BootstrapKRWCD91 : Bootstraps {

    // class를 초기화 할때 curve를 만들고, 맴버 함수를 통해서 spot등을 추출합니다.

private:

    const vector <int> Readindex = { 1,2,3,4,5,6,7 };

    string file;

    vector <Period> tenor;

    vector <Rate> rate;



    // Fixed Inputs for curve genration

    const Integer settlement = 1;

    const DayCounter innerdcc = Actual365Fixed();

    const Period coupon_period = Period(3, Months);

    const BusinessDayConvention bcc1 = ModifiedFollowing; // Fixed

    const BusinessDayConvention bcc2 = ModifiedFollowing; // Float

    const Calendar innercalendar = SouthKorea(); // Curve 생성시 Input 정의에 따른 calendar



    // Inputs for outputs

    DayCounter dc;

    Calendar calendar;

    Compounding compound;

    Frequency freq;

    BusinessDayConvention bcc;



    Date base;

    Date term;

    boost::shared_ptr <YieldTermStructure> curve;



public:

    BootstrapKRWCD91(Date, Date, string, Interpolationtype interp = linear);



    // Input setting for output

    void setting(DayCounter = Actual365Fixed(), Calendar calendar = NullCalendar(), Compounding compound = Continuous, Frequency freq = Annual, BusinessDayConvention bcc = ModifiedFollowing);



    // output method

    vector <Date> daily();

    vector <Rate> spots();

    vector <Rate> forwards();

    vector <double> tenors();

    Rate forward(Date, Date);



    boost::shared_ptr <YieldTermStructure> curveout();



    void saveTxt(string, vector <Date>, vector <double>, vector <Rate>, vector <Rate>);



    ~BootstrapKRWCD91();

};



class BootstrapUSDSofr :Bootstraps {

    // class를 초기화 할때 curve를 만들고, 맴버 함수를 통해서 spot등을 추출합니다.

private:

    const vector <int> Readindex = { 1,2,3,4,5,6,7 };

    string file;

    vector <Period> tenor;

    vector <Rate> rate;



    const Integer settlement = 2;

    const Integer pay_delay = 2;

    const DayCounter innerdcc = Actual360();

    const BusinessDayConvention bcc1 = ModifiedFollowing; // Fixed

    const BusinessDayConvention bcc2 = ModifiedFollowing; // Float

    const Calendar innercalendar = UnitedStates(UnitedStates::FederalReserve); // Curve 생성시 Input 정의에 따른 calendar



    Date base;

    Date term;

    boost::shared_ptr <YieldTermStructure> curve;



public:

    BootstrapUSDSofr(Date, Date, string, Interpolationtype interp = linear);



    // Input setting for output

    void setting(DayCounter = Actual360(), Calendar calendar = NullCalendar(), Compounding compound = Continuous, Frequency freq = Annual, BusinessDayConvention bcc = ModifiedFollowing);



    // output method

    vector <Date> daily();

    vector <Rate> spots();

    vector <Rate> forwards();

    vector <double> tenors();

    Rate forward(Date, Date);



    void saveTxt(string, vector <Date>, vector <double>, vector <Rate>, vector <Rate>);



    ~BootstrapUSDSofr();

};



class BootstrapHKDOIS :Bootstraps {

    // class를 초기화 할때 curve를 만들고, 맴버 함수를 통해서 spot등을 추출합니다.

private:

    const vector <int> Readindex = { 1,2,3,4,5,6,7 };

    string file;

    vector <Period> tenor;

    vector <Rate> rate;



    const Integer settlement = 2;

    const Integer pay_delay = 2;

    const DayCounter innerdcc = Actual365Fixed();

    const BusinessDayConvention bcc1 = ModifiedFollowing; // Fixed

    const BusinessDayConvention bcc2 = ModifiedFollowing; // Float

    const Calendar innercalendar = HongKong(); // Curve 생성시 Input 정의에 따른 calendar



    Date base;

    Date term;

    boost::shared_ptr <YieldTermStructure> curve;



public:

    BootstrapHKDOIS(Date, Date, string, Interpolationtype interp = linear);



    // Input setting for output

    void setting(DayCounter = Actual365Fixed(), Calendar calendar = NullCalendar(), Compounding compound = Continuous, Frequency freq = Quarterly, BusinessDayConvention bcc = ModifiedFollowing);



    // output method

    vector <Date> daily();

    vector <Rate> spots();

    vector <Rate> forwards();

    vector <double> tenors();

    Rate forward(Date, Date);



    void saveTxt(string, vector <Date>, vector <double>, vector <Rate>, vector <Rate>);



    ~BootstrapHKDOIS();

};



class BootstrapESTR :Bootstraps {

    // class를 초기화 할때 curve를 만들고, 맴버 함수를 통해서 spot등을 추출합니다.

private:

    const vector <int> Readindex = { 1,2,3,4,5,6,7 };

    string file;

    vector <Period> tenor;

    vector <Rate> rate;



    const Integer settlement = 2;

    const Integer pay_delay = 1;

    const DayCounter innerdcc = Actual360();

    const BusinessDayConvention bcc1 = ModifiedFollowing; // Fixed

    const BusinessDayConvention bcc2 = ModifiedFollowing; // Float

    const Calendar innercalendar = TARGET(); // Curve 생성시 Input 정의에 따른 calendar



    Date base;

    Date term;

    boost::shared_ptr <YieldTermStructure> curve;



public:

    BootstrapESTR(Date, Date, string, Interpolationtype interp = linear);



    // Input setting for output

    void setting(DayCounter = Actual360(), Calendar calendar = NullCalendar(), Compounding compound = Continuous, Frequency freq = Annual, BusinessDayConvention bcc = ModifiedFollowing);



    // output method

    vector <Date> daily();

    vector <Rate> spots();

    vector <Rate> forwards();

    vector <double> tenors();

    Rate forward(Date, Date);



    void saveTxt(string, vector <Date>, vector <double>, vector <Rate>, vector <Rate>);



    ~BootstrapESTR();

};
