% 2021.08.03
% kospi 200 미래에셋대우 16201회 조기상환형
clear; clc; close all;
tic;
% simulation time
ns = 10000;

r = 0.012; vol = 0.2044; face = 10000; S0 = 100; ref_S = 100;

strike_price = [0.95 0.95 0.95 0.9 0.9 0.9]*ref_S;
repay_n = length(strike_price);
coupon_rate = [0.025 0.05 0.075 0.1 0.125 0.15];
dummy = 0.15; ki = 0.65*ref_S;

oneyear = 360; total_day = 3*oneyear; dt = 1/oneyear;

S = zeros(total_day+1,1); S(1) = S0;
S_sub = S;
% check_day = ceil(oneyear*0.5):round(oneyear*0.5):total_day;
check_day = ceil(3*oneyear*cumsum(ones(repay_n,1))/repay_n);

payment = face*(1+coupon_rate);
payoff = zeros(repay_n,1);
total_payoff = payoff;
disc_payoff = payoff;
payoff_sub = zeros(repay_n,1);
total_payoff_sub = payoff;
disc_payoff_sub = payoff;

for i = 1:ns
    for j = 1:total_day
        eps = randn;
        S(j+1) = S(j)*exp((r - 0.5*vol^2)*dt + vol*sqrt(dt)*eps);
        S_sub(j+1) = S_sub(j)*exp((r - 0.5*vol^2)*dt + vol*sqrt(dt)*(-eps));
    end
    index = S(check_day+1);
    index_sub = S_sub(check_day +1);
    payoff = zeros(repay_n,1);
    payoff_sub = zeros(repay_n,1);
    repay_event_sub = 0;
    repay_event_sub = 0;
    
    for j = 1:repay_n
        if index(j) >= strike_price(j)
            payoff(j) = payment(j);
            repay_event = 1;
            break;
        end
    end
    for j = 1:repay_n
        if index_sub(j) >= strike_price(j)
            payoff_sub(j) = payment(j);
            repay_event_sub = 1;
            break;
        end
    end
    
    if repay_event == 0
        if min(S) < ki
            payoff(end) = face*(S(end)/ref_S);
        else
            payoff(end) = face*(1 + dummy);
        end
    end
    if repay_event_sub == 0
        if min(S_sub) < ki
            payoff_sub(end) = face*(S_sub(end)/ref_S);
        else
            payoff_sub(end) = face*(1 + dummy);
        end
    end
    total_payoff = total_payoff + 0.5*(payoff + payoff_sub);
end
total_payoff = total_payoff/ns;
for j = 1:repay_n
    disc_payoff(j) = total_payoff(j)*exp(-r*(check_day(j))*dt);
end
ELS_price = sum(disc_payoff)
cputime = toc