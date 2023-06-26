% BarrierOption.m
close all; clear all; clc; 
E=100; L=200; sigma = 0.3; r = 0.03; 
Nx = 400; dh = L/(Nx-1);

% k = T/Nt; 
T = 0.01; k = 0.0001; Nt = round(T/k); 
x = linspace(0, L, Nx);
h = diff(x); h = [h h(end)];
scale0 = sum(h);
N=Nx-2;

% Boundary condition.
u(1:Nx)=max(x-E,0); u(Nx)=0;

check_h(1) = 0.0;
for n=1:100
%     Implicit method
    % coefficients 
%%%%%%%%%%%%%%%%%%% Sample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    u(2:Nx-1)=thomas(aa, dd, cc, bb);

% % SSP-RK3
% d1(1:Nx) =0.0; d2(1:Nx) = 0.0;
% d1(2:Nx-1) = -h(2:Nx-1).*u(1:Nx-2)./(h(1:Nx-2).*(h(1:Nx-2)+h(2:Nx-1))) ... 
%%%%%%%%%%%%%%%%%%% Sample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADAPTIVE GRID GENERATION 
    if (mod(n,2) == 0)
        figure(1); hold on
        d1(1:Nx) = 0.0;
%%%%%%%%%%%%%%%%%%% Sample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %figure; plot(x(2:Nx-1),max(abs(d1(2:Nx-1)),1),'k*-')
    
    eps = 0.01;
    % monitoring function setting
%%%%%%%%%%%%%%%%%%% Sample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    M = min(M1,M2);
    
%%%%%%%%%%%%%%%%%%% Sample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
%    if (n == 2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EXACT SOLUTION
%%%%%%%%%%%%%%%%%%% Sample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % ERROR ESTIMATION
    clear error tempnorm;
    error = u - exact; 
    tempnorm = (error(2:end).^2 + error(1:end-1).^2).*h(1:end-1);
    l2norm(n) = sqrt(0.5*sum(tempnorm));  
    maxnorm(n) = max(abs(error));

  %  end
    % print out
    %if mod(n,20) == 0
    %    plot(x,u,'k-','linewidth',1)
    %end
    
end
ss = sprintf('DATA_adap_barrier_nx%d.m',Nx);
fid = fopen(ss,'wt');
fprintf(fid,'x = [\n');
fprintf(fid,'%12.10f \n',x);
fprintf(fid,'];\n');
fprintf(fid,'u=[ \n');
fprintf(fid,'%12.10f \t',u);
fprintf(fid,'];\n');
fprintf(fid,'exact=[ \n');
fprintf(fid,'%12.10f \t',exact);
fprintf(fid,'];\n');
fprintf(fid,'A=[ \n');
fprintf(fid,'%12.10f \t',l2norm);
fprintf(fid,'];\n');
fprintf(fid,'B=[ \n');
fprintf(fid,'%12.10f \t',maxnorm);
fprintf(fid,'];\n');
fclose(fid);

