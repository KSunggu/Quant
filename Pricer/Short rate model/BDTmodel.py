# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import QuantLib as ql
import pandas as pd
import Bootstrapping
from tqdm import tqdm
from scipy import optimize

class BDT:

    def __init__(self, i, B0):

        self.i = i

        self.B0 = B0

        self.ylist = []

       

    

    def BDT_Calibration(self, var):

        # Caution. This function does not calibrate for volatility, but only construct 'the tree' to satisfy martingale.

        B = np.array([np.array([0.0 for l in range(k+1)]) for k in range(self.i+3)])

        B[self.i+2][:] = 1.0

        for j in range(self.i+1,-1,-1):

            for k in range(0,j+1):

                if j == (i+1):

                    B[j][k] = (B[j+1][k]+B[j+1][k+1])*0.5/(1+(r_tree[j][k]+var)*interval_dt[j])

                else:

                    B[j][k] = (B[j+1][k]+B[j+1][k+1])*0.5/(1+(r_tree[j][k])*interval_dt[j])

        return (self.B0 - B[0][0])

   

      

    def BDT_Minimize(self, var):

        # var : [x1, x2] 

        B = np.array([np.array([0.0 for l in range(k+1)]) for k in range(self.i+3)])

        B[self.i+2][:] = 1.0

        for j in range(self.i+1,-1,-1):

            if j == (self.i+1):

                B[j][j] = (B[j+1][j]+B[j+1][j+1])*0.5/(1+(var[0])*interval_dt[j])

                for k in range(j-1,-1,-1):

                    B[j][k] = (B[j+1][k]+B[j+1][k+1])*0.5/(1+(var[0]*np.exp( (j-k)* 2.0*var[1]*np.sqrt(interval_dt[j])))*interval_dt[j])

            else:

                for k in range(0,j+1):

                    B[j][k] = (B[j+1][k]+B[j+1][k+1])*0.5/(1+(r_tree[j][k])*interval_dt[j])

 

 

        return (abs(self.B0 - B[0][0]))

 

    def BDT_Constraint(self, var):

        # var : [x1, x2] 

        B = np.array([np.array([0.0 for l in range(k+1)]) for k in range(self.i+3)])

        B[self.i+2][:] = 1.0

        for j in range(self.i+1,0,-1):

            if j == (self.i+1):

                B[j][j] = (B[j+1][j]+B[j+1][j+1])*0.5/(1+(var[0])*interval_dt[j])

                for k in range(j-1,-1,-1):

                    B[j][k] = (B[j+1][k]+B[j+1][k+1])*0.5/(1+(var[0]*np.exp( (j-k)* 2.0*var[1]*np.sqrt(interval_dt[j])))*interval_dt[j])

            else:

                for k in range(0,j+1):

                    B[j][k] = (B[j+1][k]+B[j+1][k+1])*0.5/(1+(r_tree[j][k])*interval_dt[j])

   

        # !! We suppose that the Bond as zero-coupon bond (Maturity of bond is dt_tenors[i+2])

        y1 = 4* ( (1/B[1][0])**(1/(4*(dt_tenors[self.i+2]-dt_tenors[1])))  -1);

        y2 = 4* ( (1/B[1][1])**(1/(4*(dt_tenors[self.i+2]-dt_tenors[1])))  -1);

       

        

        target = abs(0.5*np.log(y1/y2) -sigma*np.sqrt(interval_dt[self.i+1]));

               

        self.ylist += [y1, y2];

               

        return (target)

   

    def getparams(self, model, cons = None):

        errorfun = model.fun

        if (cons == None):

            errorcons = None;

        else:

            errorcons = cons['fun'](model.x)

 

        y = self.ylist

        approxvol = 0.5*np.log(y[0]/y[1]);

       

        return errorfun, errorcons, y, approxvol;

 


#%% using example
base = ql.Date(31,12,2021);
Term = ql.Date(31,12,2022);
day_fraction = ql.Actual365Fixed();

Time = [0.25, 0.5, 0.75, 1, 1.5, 2];
r = [0.03, 0.035, 0.04, 0.0425, 0.045, 0.05];
sigma = 0.2;





daily_dates, daily_tenors, fwd_rates, spot_rates = Bootstrapping.KB_Bootstrapping(base, day_fraction, Time, r, Interpolation_type = 'LogCubic', frequency = ql.Quarterly)



index = np.where(np.array([d.to_date()  for d in daily_dates]) <= Term.to_date() )[0][-1]

daily_dates = daily_dates[:index+1]
dt_index = []
step = 1;
dt_dates = daily_dates[0:-1:step]


for i in range(len(dt_dates)):

    dt_index.append(np.where(np.array(daily_dates)==dt_dates[i])[0][0])

B0 = 1/(1+np.array(spot_rates[:index+1])/4)**(np.array(daily_tenors[:index+1])*4)

B0 = B0[dt_index]

dt_tenors = np.array(daily_tenors)[dt_index]

interval_dt = np.diff(dt_tenors)

 

up = np.append([1],np.exp(sigma * np.sqrt(interval_dt)))    # len(daily_dates)

down = np.append([1],np.exp(-sigma * np.sqrt(interval_dt))) # len(daily_dates)

prob = 0.5

 

r_tree = np.array([np.array([((1+spot_rates[dt_index[1]]/4)**(4*interval_dt[0])-1)/interval_dt[0] * up[i+1] ** (i - 2*j) for j in range(i+1)]) for i in range(len(interval_dt))])

 
#%% Calibration

factor = np.zeros(len(interval_dt)-1)

 

factor2 = np.zeros([len(interval_dt)-1, 2])

# daily_dates의 수보다 2개 적다. => index-1번 해찾기를 해야함

#for i in range(index):

errorfun = []; errorcons = [];

ylist = []; approxvol = [];

for i in tqdm(range(len(factor))):

    # print(i)

    F = BDT(i,B0[i+2])

   

    ## Method 1 (This method does not calibrate volatility)

    # factor[i] = optimize.newton_krylov(F.BDT_Calibration, 0.1, f_tol=1e-10)

    # r_tree[i+1] = r_tree[i+1] + factor[i]

   

    ## Method 2

    cons = ({'type': 'eq', 'fun': F.BDT_Constraint});

    # BDT model assume that short rate does not have negative value, so we can construct bounds as 0 to inf. This property enalbes L-BFGS-B optimization algorithm. (Note that volatility also does not have negative value.)

    res = optimize.minimize(F.BDT_Minimize, [r_tree[i+1][-1], sigma], method = 'L-BFGS-B', bounds = ((0, None), (0, None)),constraints = cons, tol = 1.0e-9);

    factor2[i] = res.x

    r_tree[i+1][-1] = factor2[i][0]

    for j in range(len(r_tree[i+1])-2,-1,-1):

        r_tree[i+1][j] = r_tree[i+1][j+1]*np.exp(2.0*factor2[i][1]*np.sqrt(interval_dt[i]))

 

    temp1, temp2, temp3, temp4 = F.getparams(res, cons);

   

    errorfun.append(temp1);

    errorcons.append(temp2);   

    ylist += temp3;

    approxvol.append(temp4);

   

    del temp1, temp2, temp3, temp4

 