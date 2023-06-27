clear; clc; close all;
tic;
ns = 10000;

oneyear = 360; tot_day = 3*oneyear; dt = 1/oneyear; 
strike_rate = [0.90 0.90 0.90 0.85 0.85 0.8];
coupon_rate = cumsum(ones(6,1))*0.0375;
ki = 0.45; dummy = coupon_rate(end);

mid_day = cumsum(ones(6,1))*round(tot_day/6);
s_size = length(mid_day);

S0 = [100 100 100];
S1 = zeros(tot_day + 1, 1); S2 = S1; S3 = S2;
S1(1) = S0(1); S2(1) = S0(2); S3(1) = S0(3);
S1_s = S1; S2_s = S2; S3_s = S3;

face = 10000; r = 0.012; vol = [0.2408 0.2936 0.2254];
rho = [0.5495 0.3456 0.3681]; q = [0 0 0];

payment = face*(1 + coupon_rate);

corr_rho = [1 rho(1) rho(3); rho(1) 1 rho(2); rho(3) rho(2) 1];
M = chol(corr_rho);

tot_payoff = zeros(s_size,1);
disc_payoff = zeros(s_size,1);
for i = 1:ns
    w = randn(tot_day, 3)*M;
    for j = 1:tot_day
        S1(j+1) = S1(j)*exp((r -q(1) -0.5*vol(1)^2)*dt +vol(1)*sqrt(dt)*w(j,1));
        S1_s(j+1) = S1_s(j)*exp((r -q(1) -0.5*vol(1)^2)*dt +vol(1)*sqrt(dt)*(-w(j,1)));
        
        S2(j+1) = S2(j)*exp((r -q(2) -0.5*vol(2)^2)*dt +vol(2)*sqrt(dt)*w(j,2));
        S2_s(j+1) = S2_s(j)*exp((r -q(2) -0.5*vol(2)^2)*dt +vol(2)*sqrt(dt)*(-w(j,2)));
        
        S3(j+1) = S3(j)*exp((r -q(3) -0.5*vol(3)^2)*dt +vol(3)*sqrt(dt)*w(j,3));
        S3_s(j+1) = S3_s(j)*exp((r -q(3) -0.5*vol(3)^2)*dt +vol(3)*sqrt(dt)*(-w(j,3)));
    end
    WP = min(min(S1/S0(1), S2/S0(2)), S3/S0(3));
    WP_check = WP(mid_day + 1);
    flag = 0;
    flag_s = 0;    
    payoff = zeros(s_size,1);
    payoff_s = zeros(s_size,1);
    WP_s = min(min(S1_s/S0(1), S2_s/S0(2)), S3_s/S0(3));
    WP_check_s = WP_s(mid_day + 1);
    
    for j = 1:s_size
        if WP_check(j) >= strike_rate(j)
            payoff(j) = payment(j);
            flag = 1;
            break;
        end
    end
    for j = 1:s_size
        if WP_check_s(j) >= strike_rate(j)
            payoff_s(j) = payment(j);
            flag_s = 1;
            break;
        end
    end
    
    if flag == 0
        if min(WP) >= ki
            payoff(end) = face*(1+ dummy);
        else
            payoff(end) = face*WP(end);
        end
    end
    if flag_s == 0
        if min(WP_s) >= ki
            payoff_s(end) = face*(1+ dummy);
        else
            payoff_s(end) = face*WP_s(end);
        end
    end
%     tot_payoff = tot_payoff +payoff;
    tot_payoff = tot_payoff +0.5*(payoff +payoff_s);
end
tot_payoff = tot_payoff/ns;
for j = 1:s_size
    disc_payoff(j) = tot_payoff(j)*exp(-r*mid_day(j)*dt);
end
Price = sum(disc_payoff)

cputime = toc