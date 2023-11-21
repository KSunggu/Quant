% 2021.08.03
% 미래에셋 2-asset
clear; clc; close all;
tic;
ns = 10000;

face = 10000; r = 0.005; vol = [0.5191 0.4851];
rho = [1 0.7463; 0.7463 1]; q = [0 0];
M = chol(rho); PQ = 100;

oneyear = 360; dt = 1/oneyear; tot_day = 3*oneyear;
S0 = [195.03 97.93]; ref_S1 = S0(1); ref_S2 = S0(2);
S1 = zeros(tot_day+1,1); S2 = S1; S1(1) = S0(1); S2(1) = S0(1);
S1_s = zeros(tot_day+1,1); S2_s = S1_s; S1_s(1) = S0(1); S2_s(1) = S0(1);

strike_rate = [0.85 0.85 0.85 0.85 0.8 0.75];
repay_n = length(strike_rate);
check_day = ceil(tot_day*cumsum(ones(repay_n,1))/repay_n);
coupon_rate = cumsum(0.055*ones(repay_n,1));
dummy = coupon_rate(end); ki_rate = 0.45;

payment = face*(1 + coupon_rate);
disc_payoff = zeros(repay_n,1);
tot_payoff = zeros(repay_n,1);
for i =1:ns
    for j = 1:tot_day
        w = randn(1,2)*M;
        S1(j+1) = S1(j)*exp((r -q(1) - 0.5*vol(1)^2)*dt +vol(1)*sqrt(dt)*w(1));
        S1_s(j+1) = S1_s(j)*exp((r -q(1) - 0.5*vol(1)^2)*dt +vol(1)*sqrt(dt)*(-w(1)));
        S2(j+1) = S2(j)*exp((r -q(2) - 0.5*vol(2)^2)*dt +vol(2)*sqrt(dt)*w(2));
        S2_s(j+1) = S2_s(j)*exp((r -q(2) - 0.5*vol(2)^2)*dt +vol(2)*sqrt(dt)*(-w(2)));
    end
    WP = min(S1/ref_S1, S2/ref_S2);
    WP_check = WP(check_day +1);
    WP_s = min(S1/ref_S1, S2/ref_S2);
    WP_check_s = WP(check_day +1);
%     WP_check = min(S1(check_day +1)/ref_S(1), S2(check_day +1)/ref_S(2));
    repay_event = 0;
    payoff = zeros(repay_n,1);
    repay_event_s = 0;
    payoff_s = zeros(repay_n,1);
    
    for j = 1:repay_n
        if WP_check(j) >= strike_rate(j)
            payoff(j) = payment(j);
            repay_event = 1;
            break
        end
    end
    for j = 1:repay_n
        if WP_check_s(j) >= strike_rate(j)
            payoff_s(j) = payment(j);
            repay_event = 1;
            break
        end
    end
    
    if repay_event == 0
        if min(WP) >= ki_rate
            payoff(end) = face*(1+dummy);
        else
            payoff(end) = face*WP(end);
        end
    end
    if repay_event == 0
        if min(WP_s) >= ki_rate
            payoff_s(end) = face*(1+dummy);
        else
            payoff_s(end) = face*WP(end);
        end
    end
    tot_payoff = tot_payoff +0.5*(payoff +payoff_s);
    
%     if mod(i,ns/100) ==0
%         figure(1); 
%         subplot(2,1,1)
%         plot(1:tot_day+1,S1)
%         hold on
%         subplot(2,1,2)
%         plot(1:tot_day+1,S2)
%         hold on
%     end
end
tot_payoff = tot_payoff/ns;

for j = 1:repay_n
    disc_payoff(j) = tot_payoff(j)*exp(-r*check_day(j)*dt);
end
Price = sum(disc_payoff)
cputime = toc