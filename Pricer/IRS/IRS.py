# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 00:34:29 2023

@author: 82107
"""

import xlwings as xw
import QuantLib as ql
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import date


def to_date(date):
    return ql.Date(date.day, date.month, date.year)


def to_currency(A):
    dic = {"USD"  : ql.USDCurrency()}
    return dic[A]

def to_calendar(A):
    dic = {"USA": ql.UnitedStates()}
    return dic[A]
    
def to_DCC(A):
    dic = {"30U/360": ql.Actual360()}#ql.Thirty360(ql.Thirty360.USA) }
    return(dic[A])

def to_BDC(A):
    dic = {"Modified Following": ql.ModifiedFollowing }
    return(dic[A])

def to_index(A):
    dic = {"LIBOR 3M": ql.EURLibor3M(),
           "Fixed": ql.EURLibor3M() }
    return(dic[A])

def to_Freq(A):
    dic = {"Quarterly": ql.Quarterly}
    return(dic[A])

def to_direction(A):
    dic = {"Leg1": ql.VanillaSwap.Payer,
           "Leg2": ql.VanillaSwap.Receiver}
    return(dic[A])
#%% data imoprt

cell = xw.books("IRS pricer.xlsx")

IRS_S = cell.sheets("IRS")

Position = to_direction(IRS_S.range("C3").value) # 포지션으로 나타난 금액 지급


Valuation_date = to_date(IRS_S.range("C5").value.date())
Effective_date = to_date(IRS_S.range("C8").value.date())
Termination_date = to_date(IRS_S.range("C9").value.date())

# Leg1
Leg1 = IRS_S.range("Leg_1")

Type1 = Leg1(2,2).value
Par1 = Leg1(3,2).value
Currency1 = to_currency(Leg1(4,2).value)

Index1 = to_index(Leg1(6,2).value)
Spread1 = Leg1(7,2).value
Lever1 = Leg1(8,2).value

settle1 = to_date(Leg1(10,2).value.date())
End1 = to_date(Leg1(11,2).value.date())
Calendar1 = to_calendar(Leg1(12,2).value)

AFreq1 = to_Freq(Leg1(14,2).value)
PFreq1 = to_Freq(Leg1(15,2).value)

DCC1 = to_DCC(Leg1(16,2).value)
BDC1 = to_BDC(Leg1(17,2).value)

Tenors1 = Leg1(21,1).expand("down").value
Rate1 = Leg1(21,2).expand("down").value

swaphelpers1 = []
for i in range(len(Tenors1)):
    [tenor, rate] = [Tenors1[i], Rate1[i]]
    temphelper = ql.SwapRateHelper(ql.QuoteHandle(ql.SimpleQuote(rate)),
                                   ql.Period(int(tenor), ql.Days),
                                   Calendar1,
                                   AFreq1,
                                   BDC1,
                                   DCC1,
                                   Index1
                                   )

    swaphelpers1.append(temphelper)

Curve1 = ql.PiecewiseLogLinearDiscount(settle1, swaphelpers1, DCC1)

Schedule1 = ql.Schedule(settle1,
                        End1,
                        ql.Period(PFreq1),
                        Calendar1,
                        BDC1,
                        BDC1,
                        ql.DateGeneration.Forward,
                        False)



# Leg2
Leg2 = IRS_S.range("Leg_2")

Type2 = Leg2(2,2).value
Par2 = Leg2(3,2).value
Currency2 = to_currency(Leg2(4,2).value)

Index2 = to_index(Leg2(6,2).value)
Spread2 = Leg2(7,2).value
Lever2 = Leg2(8,2).value
Fixing_rate2= Leg2(9,2).value

settle2 = to_date(Leg2(10,2).value.date())
End2 = to_date(Leg2(11,2).value.date())
Calendar2 = to_calendar(Leg2(12,2).value)

Fixing_date2 = to_date(Leg2(13,2).value.date())
AFreq2 = to_Freq(Leg2(14,2).value)
PFreq2 = to_Freq(Leg2(15,2).value)

DCC2 = to_DCC(Leg2(16,2).value)
BDC2 = to_BDC(Leg2(17,2).value)

Tenors2 = Leg2(21,1).expand("down").value
Rate2 = Leg2(21,2).expand("down").value



swaphelpers2 = []
for i in range(len(Tenors1)):
    [tenor, rate] = [Tenors2[i], Rate2[i]]
    temphelper = ql.SwapRateHelper(ql.QuoteHandle(ql.SimpleQuote(rate)),
                                   ql.Period(int(tenor), ql.Days),
                                   Calendar2,
                                   AFreq2,
                                   BDC2,
                                   DCC2,
                                   Index2
                                   )

    swaphelpers2.append(temphelper)

Curve2 = ql.PiecewiseLogLinearDiscount(settle2, swaphelpers2, DCC2)

Schedule2 = ql.Schedule(settle2,
                        End2,
                        ql.Period(PFreq2),
                        Calendar2,
                        BDC2,
                        BDC2,
                        ql.DateGeneration.Forward,
                        False)

Index2.addFixing(Fixing_date2, Fixing_rate2)

#%% QuantLib setting
ql.Settings_instance().evaluationDate = Valuation_date

IRS = ql.VanillaSwap(Position,
                     Par1, Schedule1, Spread1, DCC1,
                     Schedule2, Index2, Spread2, DCC2)

discountCurve = ql.RelinkableYieldTermStructureHandle(Curve1)
Engine = ql.DiscountingSwapEngine(discountCurve)
IRS.setPricingEngine(Engine)

IRS.legNPV(1)


# 스케줄 생성 #


