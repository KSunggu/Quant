
#%% Package import
import numpy as np
import QuantLib as ql
import sys
import xlwings as xw
import win32com.client as winc
import glob
import os
import copy
import matplotlib.pyplot as plt
from itertools import combinations
import time
from tqdm import tqdm
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
import inspect
from datetime import date, datetime
import scipy

import Bootstrapping

# from BS_module import European_Analytic

 
def to_date(date):

    # datetime to QuantLib Date

    return ql.Date(date.day, date.month, date.year)

def Basis_Laguerre(data, degree, include_bias = True):
    # The result has the form of data_number X degree

    # Weight by exp(-X/2), alpha = 0;
    if include_bias == True:
        flag = 1;
    else:
        flag = 0;

    R = np.zeros([len(data), degree+flag])
    for i in range(flag, degree):
        R[:,i] = scipy.special.eval_genlaguerre(i,0,data)

    return R

 

#%%
frequency = ql.Continuous
bussiness_convention = ql.Unadjusted
end_of_month = False
settlement_days = 0
date_generation = ql.DateGeneration.Forward
calendar = ql.NullCalendar()
Outset_date =  date(1899,12,31)
Seed = 1;

np.random.seed(Seed)
day_count = ql.Thirty360(ql.Thirty360.BondBasis)

#%% 

Base = ql.Date(31,12,2021);
Term = ql.Date(31,12,2022);
tenors_rf = [0.25, 0.5, 0.75, 1, 1.5, 2];
rf = [0.03, 0.035, 0.04, 0.0425, 0.045, 0.05];
tenors_rj = [0.25, 0.5, 0.75, 1, 1.5, 2];
rj = [0.03, 0.035, 0.04, 0.0425, 0.045, 0.05];
sig = 0.4;
q = 0.01;
#%% Bootstrapping

daily_dates, daily_tenors, rf_fwds, rf_spots =  Bootstrapping.KTB_Bootstrapping(Base,
                                                                                day_count, tenors_rf, rf,
                                                                                Interpolation_type = 'LogCubic')
daily_dates, daily_tenors, rj_fwds, rj_spots =  Bootstrapping.KB_Bootstrapping(Base,
                                                                               day_count, tenors_rj, rj,
                                                                               Interpolation_type = 'LogCubic')

 
print("Bootstrap Clear")

 

#%% Schedule generate
Option_date = [];  Option_index = 0;
Option_start_date = Base; Option_end_date = Term;


Option_dates = ql.Schedule(to_date(Option_start_date),
                        to_date(Option_end_date),
                        ql.Period(1, ql.Days),
                        calendar,
                        bussiness_convention,
                        bussiness_convention,
                        date_generation,
                        end_of_month)

Option_dates = list(Option_dates)

while Option_dates[0] < Base:
    del Option_dates[0]

Option_dates.append(ql.Date(1,1,1901))

         

for d in daily_dates:
    if (d == Option_dates[Option_index]):
        Option_index = Option_index + 1
        Option_date.append(d)

    else:
        Option_date.append(ql.Date(1,1,1901))

Option_dates = Option_date[:]

 

 

V1_date = [];  V1_index = 0;
V1_start_date = Base; V1_end_date = Term;

V1_dates = ql.Schedule(to_date(V1_start_date),
                        to_date(V1_end_date),
                        ql.Period(1, ql.Days),
                        calendar,
                        bussiness_convention,
                        bussiness_convention,
                        date_generation,
                        end_of_month)

V1_dates = list(V1_dates)


while V1_dates[0] < Base:

    del V1_dates[0]

V1_dates.append(ql.Date(1,1,1901))

for d in daily_dates:
    if (d == V1_dates[V1_index]):
        V1_index = V1_index + 1
        V1_date.append(d)

    else:
        V1_date.append(ql.Date(1,1,1901))

V1_dates = V1_date[:]

 

print("Schedule Clear")

 

#%% Path generation

print("Start Create Stock path")

dt = np.diff(np.array(daily_tenors))
index = np.where(np.array([d for d in daily_dates])<=Term)[0][-1]
daily_dates = daily_dates[0:index+1]                                            # Trim daily_dates until option life
daily_tenors = daily_tenors[0:index+1]
Nt = len(daily_dates)

Ns = 1.0e+5;
VarReduction = 'AV';
S0 = 100; K = 100;

if VarReduction == 'AV':
    if (Ns % 2) != 0.0:
        sys.exit("Ns : even")
        
    S_AV1 = np.zeros([int(Ns/2),Nt])
    S_AV1[:,0] = S0
    S_AV2 = np.zeros([int(Ns/2),Nt])
    S_AV2[:,0] = S0

    w = np.random.randn(int(Ns/2),Nt)

    for i in range(1,Nt):

        S_AV1[:,i] = S_AV1[:,i-1]*np.exp((rf_fwds[i-1] -q -0.5*sig**2)*dt[i-1] +sig*np.sqrt(dt[i-1])*w[:,i-1]);
        S_AV2[:,i] = S_AV2[:,i-1]*np.exp((rf_fwds[i-1] -q -0.5*sig**2)*dt[i-1] -sig*np.sqrt(dt[i-1])*w[:,i-1]);

    S = np.vstack([S_AV1,S_AV2])
    Ns = len(S[:,0])
    del S_AV1, S_AV2

else:      

    S = np.zeros([Ns,Nt])
    S[:,0] = S0
    w = np.random.randn(Ns,Nt)

    for i in range(1,Nt):
        S[:,i] = S[:,i-1]*np.exp((rf_fwds[i-1] -q -0.5*sig**2)*dt[i-1] +sig*np.sqrt(dt[i-1])*w[:,i-1]);


 

#%% Strike Price

Strike = []

for i, d in enumerate(daily_dates):
    if d == Option_dates[i]:
        Strike.append(K);

#%% Vesting condition
Vesting_flag = 'Off'; V1_Barriers = [100, 100];
if Vesting_flag == 'On':
    V1_flag = np.zeros([Ns,Nt]);
    Universial_set = np.where(S[:,0] >= -1.0e-32)[0]
    V1_idx = 0;

    for j in tqdm(range(Nt)):
        if daily_dates[j] ==  V1_dates[j]:
            # There is no duplicate element in each array. So for the computation time, we suggest that assume_unique to be True
            Sub_set = np.setdiff1d(Universial_set, np.where(S[:,j] < V1_Barriers[V1_idx])[0], True)
            V1_flag[Sub_set,j:] = 1.0
            Universial_set = np.setdiff1d(Universial_set, np.where(S[:,j] >= V1_Barriers[V1_idx])[0], True) # Find the paths that are not satisfied the condition.
            V1_idx += 1;

    del Universial_set, Sub_set, V1_idx;
    
else:
    V1_flag = np.ones([Ns,Nt]);

 

#%% Backward propagation
O_type = "Call";

if O_type == "Call":
    flag = 1;

elif O_type == "Put":
    flag = -1;


Hv = np.zeros([Ns,Nt])
Cv = np.zeros([Ns,Nt])
decision = np.empty([Ns,Nt], dtype = str )

Basis = "Laguerre";
Degree = 3;

Base_model = LinearRegression() # 회귀 모형

O_idx = len(Strike)-1
Payoff = np.maximum(flag*(S[:,-1] -Strike[O_idx]),0)*V1_flag[:,-1]
Hv[:,-1] = Payoff
O_idx -= 1;

for j in tqdm(range(Nt-2,-1,-1)):
    model= copy.deepcopy(Base_model)
    Hv[:,j] = Hv[:,j+1]/(1+rj_fwds[j]*dt[j])

    if daily_dates[j] == Option_dates[j]:
        Payoff = np.maximum( flag*(S[:,j] -Strike[O_idx]), 0)*V1_flag[:,j]
        ITM_path = np.where( Payoff > 0.0)[0]

        if len(ITM_path) == 0:
            pass;

        else:

            if Basis == "Taylor":
                x_data = PolynomialFeatures(int(Degree), include_bias = True).fit_transform(S[ITM_path,j].reshape(-1,1)); 
                
            elif Basis == "Laguerre":
                x_data = Basis_Laguerre(S[ITM_path,j], int(Degree), include_bias = True)

            model.fit(x_data, Hv[ITM_path,j].reshape(-1,1))
            Cv[ITM_path,j] = model.predict(x_data).reshape(1,-1)

   

        for i in ITM_path: 
            Payoff = max([max(flag*(S[i,j] -K),0)*V1_flag[i,j], 0, Cv[i,j]])

            if Payoff == max(flag*(S[i,j] -K),0)*V1_flag[i,j]:
                Hv[i,j] = max(flag*(S[i,j] -K),0)*V1_flag[i,j]
                decision[i,j] = 'C'

            elif Payoff == Cv[i,j]:
                Hv[i,j] = Cv[i,j]
                decision[i,j] = 'H'

            else:
                Hv[i,j] = 0.0
                decision[i,j] = 'H'

        O_idx -= 1;
    else:

        decision[i,j] = 'H'

 

Price = np.mean(Hv[:,0])
