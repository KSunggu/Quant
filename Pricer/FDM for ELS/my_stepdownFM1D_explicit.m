clear; clc; close all;
% Boundary condtion : Left : du = 1; Right : ddu = 0
face = 10000; E = 250; L = 3*E; sig = 0.2044; r = 0.0012; Nx = L;
h = L/Nx; x = linspace(0.5*h,L -0.5*h,Nx)'; h2 = h^2;

% round¸¦ ÀÌ¿ëÇÏ¸é Á¤È®µµ°¡ ¶³¾îÁü?? ±â°£ÀÌ µü ¾È³ª´²Áü..
oneyear = 300; T = 3; dt = 1/oneyear; Nt = round(T/dt);
mid_day = 150*cumsum(ones(6,1));
mid_day(end) = mid_day(end) +2; % tag index¸¦ À§ÇØ
% mid_day = [round(Nt/6) round(2*Nt/6) round(3*Nt/6) round(4*Nt/6) round(5*Nt/6) Nt+2];

coupon = [0.15 0.125 0.10 0.075 0.05 0.025];
R = [0.90 0.90 0.90 0.95 0.95 0.95];
dummy = 0.15; 
kib = 0.65;

% initial condition
u = face*(1+coupon(1))*ones(Nx,1);
ku = u;
for i = 1:Nx
    if x(i) < E*0.65
        u(i) = face*(x(i)/E);
        ku(i) = u(i);
    elseif x(i) < E*0.9
        u(i) = face*(1+dummy);
        ku(i) = face*(x(i)/E);
    end
end

% plot initial condition
figure(1); hold on
plot(x,u,'ko-'); plot(x,ku,'r*-');
tag = 1;
for i = 1:Nt
    if i*dt == mid_day(tag)
        idx = min(find(x >= E*R(tag+1)));
        u(idx:end) = face*(1 + coupon(tag+1));
        ku(idx:end) = face*(1 + coupon(tag+1));
        tag = tag+1; 
    end
    idx = min(find(x >= E*kib));
    u(1:idx-1)= ku(1:idx-1);
    
    u(2:end-1) = u(2:end-1) ...
        + 0.5*dt*sig^2*x(2:end-1).^2.*(u(1:end-2) -2*u(2:end-1) + u(3:end))/h2 ...
        + 0.5*dt*r*x(2:end-1).*(u(3:end) -u(1:end-2))/h - dt*r*u(2:end-1);
    ku(2:end-1) = ku(2:end-1) ...
        + 0.5*dt*sig^2*x(2:end-1).^2.*(ku(1:end-2) -2*ku(2:end-1) + ku(3:end))/h2 ...
        + 0.5*dt*r*x(2:end-1).*(ku(3:end) -ku(1:end-2))/h - dt*r*ku(2:end-1);
    u(end) = 2*u(end-1) -u(end-2); u(1) = -u(2);
    ku(end) = 2*ku(end-1) -u(end-2); ku(1) = -ku(2);
end

% plot final
figure(2); hold on
plot(x,u,'ko-'); 
figure(3)
plot(x,ku,'r*-')

idx = min(find(x >= E));
Price = 0.5*(u(idx) +ku(idx))
