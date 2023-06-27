import numpy as np
import QuantLib as ql

def KTB_Bootstrapping(base_date, day_fraction, Time, r, Interpolation_type = 'LogCubic', frequency = ql.Semiannual):
    
    # Bond rates
    bond_maturities = [ql.Period(int(i*12), ql.Months) for i in Time ]          # n year-> n x 4 + 1
    
    bond_rates = r
    
    #print_curve(bond_maturities, bond_rates)
    
    # some constants and conventions
    # here we just assume for the sake of example
    # that some of the constants are the same for
    # depo rates and bond rates
    
    calc_date = base_date                                                       # ql.Date(31, 12, 2019)
    ql.Settings.instance().evaluationDate = calc_date
    
    calendar =  ql.NullCalendar()
    bussiness_convention = ql.Unadjusted
    #day_count = ql.Actual365Fixed()
    day_count = day_fraction                                                    #ql.ActualActual()
    end_of_month = True
    settlement_days = 0
    face_amount = 100
    coupon_frequency = ql.Period(ql.Semiannual)
    settlement_days = 0
    
    
    # create fixed rate bond helpers from fixed rate bonds
    bond_helpers = []
    for r, m in zip(bond_rates, bond_maturities):
        termination_date = calc_date + m
        if calendar.isEndOfMonth(calc_date)==True:
            termination_date = calendar.endOfMonth(termination_date)
        schedule = ql.Schedule(calc_date,
                       termination_date,
                       coupon_frequency,
                       calendar,
                       bussiness_convention,
                       bussiness_convention,
                       ql.DateGeneration.Backward,
                       end_of_month)
    
        helper = ql.FixedRateBondHelper(ql.QuoteHandle(ql.SimpleQuote(face_amount)),
                                            settlement_days,
                                            face_amount,
                                            schedule,
                                            [r],
                                            ql.Thirty360(ql.Thirty360.BondBasis),
                                            bussiness_convention,
                                            )
        bond_helpers.append(helper)
        
    rate_helpers = bond_helpers
    
    # Choise option for interpolation**********************************************
    if Interpolation_type == 'Linear':
        yieldcurve = ql.PiecewiseLinearZero(calc_date, rate_helpers, day_count)
    elif Interpolation_type == 'LogCubic':
        yieldcurve = ql.PiecewiseLogCubicDiscount(calc_date, rate_helpers, day_count)
    elif Interpolation_type == 'Cubic':   
        yieldcurve = ql.PiecewiseCubicZero(calc_date, rate_helpers, day_count)
    
   
    yieldcurve.enableExtrapolation()
    daily_dates = [ ql.Date(i) for i in range(calc_date.serialNumber(), termination_date.serialNumber()+1) ]
    daily_tenors = []
    tenors = []
    spot_rates = []
    for d in daily_dates:
        yrs = day_count.yearFraction(calc_date, d)
        compounding = ql.Compounded
        freq = frequency
        zero_rate = yieldcurve.zeroRate(yrs, compounding, freq)
        
        tenors.append(yrs)
        eq_rate = zero_rate.equivalentRate(day_count,
                                           compounding,
                                           freq,
                                           calc_date,
                                           d).rate()
        spot_rates.append(eq_rate)
    
    for d in daily_dates:
        dailys = day_count.yearFraction(calc_date,d)
        daily_tenors.append(dailys)
        
#    fwd_rates = [ 0.0 if day_count.yearFraction(d, calendar.advance(d,1,ql.Days)) == 0.0 else yieldcurve.forwardRate(d, calendar.advance(d,1,ql.Days),
#                                                   day_count, ql.Simple).rate()
#                         for d in daily_dates ]
    
    fwd_rates = [ yieldcurve.forwardRate(d, calendar.advance(d,1,ql.Days), day_count, ql.Simple).rate()
                         for d in daily_dates ]
    
    
    del fwd_rates[-1]

    return daily_dates, daily_tenors, fwd_rates, spot_rates

def KB_Bootstrapping(base_date, day_fraction, Time, r,  Interpolation_type = 'LogCubic', frequency = ql.Quarterly):
    
    
    # Bond rates
    bond_maturities = [ql.Period(int(i*12), ql.Months) for i in Time ]          # n year-> n x 4 + 1
    
    bond_rates = r
    
    #print_curve(bond_maturities, bond_rates)
    
    # some constants and conventions
    # here we just assume for the sake of example
    # that some of the constants are the same for
    # depo rates and bond rates
    
    calc_date = base_date                                                       # ql.Date(31, 12, 2019)
    ql.Settings.instance().evaluationDate = calc_date
    
    calendar =  ql.NullCalendar()
    bussiness_convention = ql.Unadjusted
    #day_count = ql.Actual365Fixed()
    day_count = day_fraction                                                    #ql.ActualActual()
    end_of_month = True
    settlement_days = 0
    face_amount = 100
    coupon_frequency = ql.Period(ql.Quarterly)
    settlement_days = 0
    
    
    # create fixed rate bond helpers from fixed rate bonds
    bond_helpers = []
    for r, m in zip(bond_rates, bond_maturities):
        termination_date = calc_date + m
        if calendar.isEndOfMonth(calc_date)==True:
            termination_date = calendar.endOfMonth(termination_date)
        schedule = ql.Schedule(calc_date,
                       termination_date,
                       coupon_frequency,
                       calendar,
                       bussiness_convention,
                       bussiness_convention,
                       ql.DateGeneration.Backward,
                       end_of_month)
    
        helper = ql.FixedRateBondHelper(ql.QuoteHandle(ql.SimpleQuote(face_amount)),
                                            settlement_days,
                                            face_amount,
                                            schedule,
                                            [r],
                                            ql.Thirty360(ql.Thirty360.BondBasis),
                                            bussiness_convention,
                                            )
        bond_helpers.append(helper)
        
    rate_helpers = bond_helpers
    
    # Choise option for interpolation**********************************************
    if Interpolation_type == 'Linear':
        yieldcurve = ql.PiecewiseLinearZero(calc_date, rate_helpers, day_count)
    elif Interpolation_type == 'LogCubic':
        yieldcurve = ql.PiecewiseLogCubicDiscount(calc_date, rate_helpers, day_count)
    elif Interpolation_type == 'Cubic':   
        yieldcurve = ql.PiecewiseCubicZero(calc_date, rate_helpers, day_count)
    
   
    yieldcurve.enableExtrapolation()
    daily_dates = [ ql.Date(i) for i in range(calc_date.serialNumber(), termination_date.serialNumber()+1) ]
    daily_tenors = []
    tenors = []
    spot_rates = []
    for d in daily_dates:
        yrs = day_count.yearFraction(calc_date, d)
        compounding = ql.Compounded
        freq = frequency
        zero_rate = yieldcurve.zeroRate(yrs, compounding, freq)
        
        tenors.append(yrs)
        eq_rate = zero_rate.equivalentRate(day_count,
                                           compounding,
                                           freq,
                                           calc_date,
                                           d).rate()
        spot_rates.append(eq_rate)
    
    for d in daily_dates:
        dailys = day_count.yearFraction(calc_date,d)
        daily_tenors.append(dailys)
        
#    fwd_rates = [ 0.0 if day_count.yearFraction(d, calendar.advance(d,1,ql.Days)) == 0.0 else yieldcurve.forwardRate(d, calendar.advance(d,1,ql.Days),
#                                                   day_count, ql.Simple).rate()
#                         for d in daily_dates ]
    
    fwd_rates = [ yieldcurve.forwardRate(d, calendar.advance(d,1,ql.Days), day_count, ql.Simple).rate()
                         for d in daily_dates ]
    
    
    del fwd_rates[-1]

    return daily_dates, daily_tenors, fwd_rates, spot_rates

def SOFR_Bootstrapping(base_date, day_fraction, Term_dates, r,  Interpolation_type = 'LogCubic', frequency = ql.Semiannual):
    """
    @author: Sunggu Kim
    USD SOFR 금리에 대한 Curve를 생성합니다.
    
    Term_dates에는 SOFR swap 만기일이 들어감에 주의
    settlement day는 ISDA에 따른 일반적인 경우만 고려하여 2일로 고려
    
    """
    
    settlement_day = 2;
    calendar = ql.NullCalendar();
    Index = ql.Sofr();
    DCC = ql.Actual360();
    Term_days = [calendar.businessDaysBetween(base_date, d) for d in Term_dates];
    
    helpers = []
    for d, r in zip(Term_days, r):
        helper = ql.OISRateHelper(settlement_day, ql.Period(d,ql.Days), ql.QuoteHandle(ql.SimpleQuote(r)), Index)
        
        helpers.append(helper)
    

    Interpolations = {'Logcubic': ql.PiecewiseLogCubicDiscount(base_date, helpers, DCC), 
                      'Linear' :  ql.PiecewiseLinearZero(base_date, helpers, DCC)}
    
    curve = Interpolations[Interpolation_type]
    curve.enableExtrapolation()
    
    daily_dates = [ ql.Date(serial) for serial in range(base_date.serialNumber(), Term_dates[-1].serialNumber()+1) ]
    
    daily_tenors = []
    rf_spots = []
    for d in daily_dates:
        yrs = DCC.yearFraction(base_date, d)
        compounding = ql.Compounded
        freq = frequency
        zero_rate = curve.zeroRate(yrs, compounding, freq)
        
        daily_tenors.append(yrs)
        eq_rate = zero_rate.equivalentRate(DCC,
                                           compounding,
                                           freq,
                                           base_date,
                                           d).rate()
        rf_spots.append(eq_rate)
        
    rf_fwds = [ curve.forwardRate(d, calendar.advance(d,1,ql.Days), DCC, ql.Simple).rate()
                         for d in daily_dates ]
    
    
    return daily_dates, daily_tenors, rf_fwds, rf_spots