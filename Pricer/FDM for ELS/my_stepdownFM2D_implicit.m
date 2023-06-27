clc; close all; clear;
face = 10000; r = 0.005; sig1 = 0.5191; sig2 = 0.4851;
rho = 0.7463; q1 = 0; q2 = 0; E1 = 195.03; E2 = 97.93;
L1 = E1*3; L2 = E2*3; Nx1 = ceil(L1); Nx2 = ceil(L2);
Nx = max(Nx1, Nx2); h = 1/Nx; h2 = h^2;
x1 = linspace(0.5*h, L1 -0.5*h, Nx);
x2 = linspace(0.5*h, L2 -0.5*h, Nx);
oneyear = 300; T = 3; dt = 1/oneyear; Nt = round(T/dt);

R = [0.75 0.8 0.85 0.85 0.85 0.85];
step = [round(Nt/6) round(2*Nt/6) round(3*Nt/6) round(4*Nt/6) round(5*Nt/6) Nt+2];
coupon = cumsum(0.055*ones(6,1)); coupon = sort(coupon,'descend');
dummy = coupon(end); kib = 0.45;

u = face*(1 +coupon(1))*ones(Nx,Nx); ku = u;
for i = 1:Nx
    for j = 1:Nx
        wp = min(x1(i)/E1,x2(j)/E2);
        if wp < kib
            u(i,j) = face*wp;
            ku(i,j) = u(i);
        elseif wp < R(1)
            u(i,j) = face*(1+ dummy);
            ku(i,j) = face*wp;
        end
    end
end
nu = u; nku = ku;
% plot initial condition
figure(1);
mesh(x1,x2,u);
figure(2);
mesh(x1,x2,ku);

main1 = zeros(1,Nx); sub1 = zeros(1,Nx); sup1 = zeros(1,Nx);
main2 = zeros(1,Nx); sub2 = zeros(1,Nx); sup2 = zeros(1,Nx);
b1 = zeros(1,Nx); b2 = zeros(1,Nx);
for i = 1:Nx
    main1(i) = 1/dt +0.5*r + sig1^2*x1(i)^2/h2;
    sub1(i) = -0.5*sig1^2*x1(i)^2/h2 +0.5*r*x1(i)/h;
    sup1(i) = -0.5*sig1^2*x1(i)^2/h2 -0.5*r*x1(i)/h;
end
main1(1) = main1(1) -sub1(1);
main1(Nx) = main1(Nx) +2*sup1(Nx);
sub1(Nx) = sub1(Nx) -sup1(Nx);

for i =1:Nx
    main2(i) = 1/dt +0.5*r + sig2^2*x2(i)^2/h2;
    sub2(i) = -0.5*sig2^2*x2(i)^2/h2 +0.5*r*x2(i)/h;
    sup2(i) = -0.5*sig2^2*x2(i)^2/h2 -0.5*r*x2(i)/h;
end
main2(1) = main2(1) -sub1(1);
main2(Nx) = main2(Nx) +2*sup1(Nx);
sub2(Nx) = sub2(Nx) -sup2(Nx);

tag = 1;
for n = 1:Nt
    if n == step(tag)
        s1 = min(find(x1 > E1*R(tag+1)));
        s2 = min(find(x2 > E2*R(tag+1)));
        u(s1:end,s2:end) = face*(1+coupon(tag+1));
        ku(s1:end,s2:end) = face*(1+coupon(tag+1));
        tag = tag+1;
    end
% % % %     u calculate
    %%%% step1
    for i =2:Nx-1
        b1(i) = u(i,1)/dt +0.5*rho*x1(i)*x2(1)*sig1*sig2*...
            (u(i+1,1+1) +u(i+1,1) -u(i-1,1+1) -u(i-1,1));
    end
    b1(1) = u(1,1)/dt +0.5*rho*x1(1)*x2(1)*sig1*sig2*...
        (u(1+1,1+1) +u(1+1,1) +u(1,1+1) +u(1,1));
    b1(Nx) = u(Nx,1)/dt +0.5*rho*x1(Nx)*x2(1)*sig1*sig2*...
        (2*u(Nx,1+1) -u(Nx-1,1+1) -(-2*u(Nx,1) +u(Nx-1,1)) -u(Nx-1,1+1) -u(Nx-1,1));
    nu(:,1) = thomas(sub1,main1,sup1,b1);
 
    for j = 2:Nx-1
        for i =2:Nx-1
            b1(i) = u(i,j)/dt +0.5*rho*x1(i)*x2(j)*sig1*sig2*...
                (u(i+1,j+1) -u(i+1,j-1) -u(i-1,j+1) +u(i-1,j-1));
        end
        b1(1) = u(1,j)/dt +0.5*rho*x1(1)*x2(j)*sig1*sig2*...
            (u(1+1,j+1) -u(1+1,j-1) +u(1,j+1) -u(1,j-1));
        b1(Nx) = u(Nx,j)/dt +0.5*rho*x1(Nx)*x2(j)*sig1*sig2*...
            (2*u(Nx,j+1) -u(Nx-1,j+1) -(2*u(Nx,j-1) -u(Nx-1,j-1)) -u(Nx-1,j+1) +u(Nx-1,j-1));
        nu(:,j) = thomas(sub1,main1,sup1,b1);
    end
    
    for i =2:Nx-1
        b1(i) = u(i,Nx)/dt +0.5*rho*x1(i)*x2(Nx)*sig1*sig2*...
            (2*u(i+1,Nx) -u(i+1,Nx-1) -u(i+1,Nx) -(2*u(i-1,Nx) -u(i-1,Nx)) +u(i-1,Nx-1));
    end
    b1(1) = u(1,Nx)/dt +0.5*rho*x1(1)*x2(Nx)*sig1*sig2*...
        (2*u(1+1,Nx) -u(1+1,Nx-1) -u(1+1,Nx) -(-2*u(1,Nx) +u(1,Nx)) -u(1,Nx-1));
    b1(Nx) = u(Nx,Nx)/dt +0.5*rho*x1(Nx)*x2(Nx)*sig1*sig2*...
        ((4*u(Nx,Nx) -2*u(Nx-1,Nx) -2*u(Nx,Nx-1) + u(Nx-1,Nx-1))-(2*u(Nx,Nx-1) -u(Nx-1,Nx-1)) -(2*u(Nx-1,Nx) -u(Nx-1,Nx-1)) +u(Nx-1,Nx-1));
    nu(:,Nx) = thomas(sub1,main1,sup1,b1);
    
    %%%% step2 
    for j =2:Nx-1
        b2(j) = nu(1,j)/dt +0.5*rho*x1(1)*x2(j)*sig1*sig2*...
            (nu(1+1,j+1) -nu(1+1,j-1) +nu(1,j+1) -nu(1,j-1));
    end
    b2(1) = nu(1,1)/dt +0.5*rho*x1(1)*x2(1)*sig1*sig2*...
        (nu(1+1,1+1) +nu(1+1,1) +nu(1,1+1) +nu(1,1));
    b2(Nx) = nu(1,Nx)/dt +0.5*rho*x1(1)*x2(Nx)*sig1*sig2*...
        (2*nu(1+1,Nx) -nu(1+1,Nx-1) -nu(1+1,Nx-1) -(-2*nu(1,Nx) +nu(1,Nx-1)) -nu(1,Nx-1));
    u(1,:) = thomas(sub2,main2,sup2,b2);
 
    for i = 2:Nx-1
        for j =2:Nx-1
            b2(j) = nu(i,j)/dt +0.5*rho*x1(i)*x2(j)*sig1*sig2*...
                (nu(i+1,j+1) -nu(i+1,j-1) -nu(i-1,j+1) +nu(i-1,j-1));
        end
        b2(1) = nu(i,1)/dt +0.5*rho*x1(i)*x2(1)*sig1*sig2*...
            (nu(i+1,1+1) +nu(i+1,1) -nu(i-1,1+1) -nu(i-1,1));
        b2(Nx) = nu(i,Nx)/dt +0.5*rho*x1(i)*x2(Nx)*sig1*sig2*...
            ((2*nu(i+1,Nx) -nu(i+1,Nx)) -nu(i+1,Nx-1) -(2*nu(i-1,Nx) -nu(i-1,Nx-1)) +nu(i-1,Nx-1));
        u(i,:) = thomas(sub2,main2,sup2,b2);
    end
    
    for j =2:Nx-1
        b2(j) = nu(Nx,j)/dt +0.5*rho*x1(Nx)*x2(j)*sig1*sig2*...
            (2*nu(Nx,j+1) -nu(Nx-1,j+1) -(2*nu(Nx,j-1) -nu(Nx,j-1)) -nu(Nx-1,j+1) +nu(Nx-1,j-1));
    end
    b2(1) = nu(Nx,1)/dt +0.5*rho*x1(Nx)*x2(1)*sig1*sig2*...
        (2*nu(Nx,1+1) -nu(Nx-1,1+1) -(-2*nu(Nx,1) +nu(Nx,1)) -nu(Nx-1,1+1) -nu(Nx-1,1));
    b2(Nx) = nu(Nx,Nx)/dt +0.5*rho*x1(Nx)*x2(Nx)*sig1*sig2*...
        (4*nu(Nx,Nx) -2*nu(Nx-1,Nx) -2*nu(Nx,Nx-1) +nu(Nx-1,Nx-1) -(2*nu(Nx,Nx-1) -nu(Nx-1,Nx-1)) -(2*nu(Nx-1,Nx) -nu(Nx-1,Nx-1)) +nu(Nx-1,Nx-1));
    u(Nx,:) = thomas(sub2,main2,sup2,b2);
    
    % % % %     ku calculate
    %%%% step1
    for i =2:Nx-1
        b1(i) = ku(i,1)/dt +0.5*rho*x1(i)*x2(1)*sig1*sig2*...
            (ku(i+1,1+1) +ku(i+1,1) -ku(i-1,1+1) -ku(i-1,1));
    end
    b1(1) = ku(1,1)/dt +0.5*rho*x1(1)*x2(1)*sig1*sig2*...
        (ku(1+1,1+1) +ku(1+1,1) +ku(1,1+1) +ku(1,1));
    b1(Nx) = ku(Nx,1)/dt +0.5*rho*x1(Nx)*x2(1)*sig1*sig2*...
        (2*ku(Nx,1+1) -ku(Nx-1,1+1) -(-2*ku(Nx,1) +ku(Nx-1,1)) -ku(Nx-1,1+1) -ku(Nx-1,1));
    nku(:,1) = thomas(sub1,main1,sup1,b1);
 
    for j = 2:Nx-1
        for i =2:Nx-1
            b1(i) = ku(i,j)/dt +0.5*rho*x1(i)*x2(j)*sig1*sig2*...
                (ku(i+1,j+1) -ku(i+1,j-1) -ku(i-1,j+1) +ku(i-1,j-1));
        end
        b1(1) = ku(1,j)/dt +0.5*rho*x1(1)*x2(j)*sig1*sig2*...
            (ku(1+1,j+1) -ku(1+1,j-1) +ku(1,j+1) -ku(1,j-1));
        b1(Nx) = ku(Nx,j)/dt +0.5*rho*x1(Nx)*x2(j)*sig1*sig2*...
            (2*ku(Nx,j+1) -ku(Nx-1,j+1) -(2*ku(Nx,j-1) -ku(Nx-1,j-1)) -ku(Nx-1,j+1) +ku(Nx-1,j-1));
        nku(:,j) = thomas(sub1,main1,sup1,b1);
    end
    
    for i =2:Nx-1
        b1(i) = ku(i,Nx)/dt +0.5*rho*x1(i)*x2(Nx)*sig1*sig2*...
            (2*ku(i+1,Nx) -ku(i+1,Nx-1) -ku(i+1,Nx) -(2*ku(i-1,Nx) -ku(i-1,Nx)) +ku(i-1,Nx-1));
    end
    b1(1) = ku(1,Nx)/dt +0.5*rho*x1(1)*x2(Nx)*sig1*sig2*...
        (2*ku(1+1,Nx) -ku(1+1,Nx-1) -ku(1+1,Nx) -(-2*ku(1,Nx) +ku(1,Nx)) -ku(1,Nx-1));
    b1(Nx) = ku(Nx,Nx)/dt +0.5*rho*x1(Nx)*x2(Nx)*sig1*sig2*...
        ((4*ku(Nx,Nx) -2*ku(Nx-1,Nx) -2*ku(Nx,Nx-1) + ku(Nx-1,Nx-1))-(2*ku(Nx,Nx-1) -ku(Nx-1,Nx-1)) -(2*ku(Nx-1,Nx) -ku(Nx-1,Nx-1)) +ku(Nx-1,Nx-1));
    nku(:,Nx) = thomas(sub1,main1,sup1,b1);
    
    %%%% step2 
    for j =2:Nx-1
        b2(j) = nku(1,j)/dt +0.5*rho*x1(1)*x2(j)*sig1*sig2*...
            (nku(1+1,j+1) -nku(1+1,j-1) +nku(1,j+1) -nku(1,j-1));
    end
    b2(1) = nku(1,1)/dt +0.5*rho*x1(1)*x2(1)*sig1*sig2*...
        (nku(1+1,1+1) +nku(1+1,1) +nku(1,1+1) +nku(1,1));
    b2(Nx) = nku(1,Nx)/dt +0.5*rho*x1(1)*x2(Nx)*sig1*sig2*...
        (2*nku(1+1,Nx) -nku(1+1,Nx-1) -nku(1+1,Nx-1) -(-2*nku(1,Nx) +nku(1,Nx-1)) -nku(1,Nx-1));
    ku(1,:) = thomas(sub2,main2,sup2,b2);
 
    for i = 2:Nx-1
        for j =2:Nx-1
            b2(j) = nku(i,j)/dt +0.5*rho*x1(i)*x2(j)*sig1*sig2*...
                (nku(i+1,j+1) -nku(i+1,j-1) -nku(i-1,j+1) +nku(i-1,j-1));
        end
        b2(1) = nku(i,1)/dt +0.5*rho*x1(i)*x2(1)*sig1*sig2*...
            (nku(i+1,1+1) +nku(i+1,1) -nku(i-1,1+1) -nku(i-1,1));
        b2(Nx) = nku(i,Nx)/dt +0.5*rho*x1(i)*x2(Nx)*sig1*sig2*...
            ((2*nku(i+1,Nx) -nku(i+1,Nx)) -nku(i+1,Nx-1) -(2*nku(i-1,Nx) -nku(i-1,Nx-1)) +nku(i-1,Nx-1));
        ku(i,:) = thomas(sub2,main2,sup2,b2);
    end
    
    for j =2:Nx-1
        b2(j) = nku(Nx,j)/dt +0.5*rho*x1(Nx)*x2(j)*sig1*sig2*...
            (2*nku(Nx,j+1) -nku(Nx-1,j+1) -(2*nku(Nx,j-1) -nku(Nx,j-1)) -nku(Nx-1,j+1) +nku(Nx-1,j-1));
    end
    b2(1) = nku(Nx,1)/dt +0.5*rho*x1(Nx)*x2(1)*sig1*sig2*...
        (2*nku(Nx,1+1) -nku(Nx-1,1+1) -(-2*nku(Nx,1) +nku(Nx,1)) -nku(Nx-1,1+1) -nku(Nx-1,1));
    b2(Nx) = nku(Nx,Nx)/dt +0.5*rho*x1(Nx)*x2(Nx)*sig1*sig2*...
        (4*nku(Nx,Nx) -2*nku(Nx-1,Nx) -2*nku(Nx,Nx-1) +nku(Nx-1,Nx-1) -(2*nku(Nx,Nx-1) -nku(Nx-1,Nx-1)) -(2*nku(Nx-1,Nx) -nku(Nx-1,Nx-1)) +nku(Nx-1,Nx-1));
    ku(Nx,:) = thomas(sub2,main2,sup2,b2);
end

figure(3)
mesh(x1,x2,u);
figure(4)
mesh(x1,x2,ku);

wp = max(x1/E1,x2/E2);
idx = min(find(wp > 1));
u(idx)