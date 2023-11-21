clear; clc; close all;

val_date = datenum('11-25-2016');
mid_date = datenum(['05-23-2017'; '11-22-2017'; '05-23-2018'; ...
    '11-23-2018'; '05-22-2019'; '11-21-2019']);
c_rate = [0.025 0.05 0.075 0.1 0.125 0.15];
mid_ch = (mid_date -val_date);
mid_size = length(mid_ch);
strike = [0.95 0.95 0.95 0.9 0.9 0.9];
N = mid_ch(end);
dt = 1/365;
ki = 0.65;
dummy = 0.15;

face = 10000;
r = 0.012; sig = 0.2044;
S0 = 100;
ref_S = S0;

ns = 10000;
payment = zeros(ns,mid_size);
for k = 1:mid_size
    payment(:,k) = face*(1+c_rate(k));
end

S = zeros(ns,N+1);
S(:,1) = S0;
for k = 1:ns
    for m = 1:N
        S(k,m+1) = S(k,m)*exp((r- sig^2/2)*dt +sig*sqrt(dt)*randn(1));
    end
end

% plot(1:N+1,S)

WP = S/ref_S;
WP_check = WP(:,mid_ch +1);

payoff = zeros(ns,mid_size);
for k = 1:ns
    for m = 1:mid_size
        if WP_check(k,m) >= strike(m)
            payoff(k,m) = payment(k,m);
            break;
        end
    end
end

for k =1:ns
   if payoff(k,:) == 0
        ki_event = any(WP(k,:) < ki);
        if ki_event == 1
            payoff(k,end) = face*WP(k,end);
        else
            payoff(k,end) = payment(k,end);
        end
    end
end

disc_payoff = zeros(ns,1);
for k = 1:ns
    for m = 1:mid_size
        if payoff(k,m) ~= 0
            disc_payoff(k,1) = payoff(k,m)*exp(-r*mid_ch(m)*dt);
        end
    end
end

% exp_payoff = mean(payoff);
% disc_payoff = zeros(1, mid_size);
% for m = 1:mid_size
%     disc_payoff(m) = exp_payoff(m)*exp(-r*mid_ch(m)*dt);
% end

ELS_payoff = mean(disc_payoff)