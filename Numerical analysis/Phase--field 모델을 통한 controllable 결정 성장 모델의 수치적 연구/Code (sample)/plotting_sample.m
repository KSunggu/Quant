clear all; clc; close all;
cd data_90_10_K10
remarks; cd ..
fs = 22;
xright = h*nx;
x = linspace(-0.5*xright, 0.5*xright, nx+1);
[xx, yy] = meshgrid(x,x);

fig = figure(1); clf; hold on
for i = [1:3:20]
    ss = sprintf('./data_90_10_K10/cry%d.m',i);
    C = load(ss);
    axis image
    contour(xx, yy, C, [0.5 0.5], 'k-');
end
axis image
axis([-0.5*xright 0.5*xright -0.5*xright 0.5*xright])
set(gca,'fontsize',fs)
box on

print -deps fig_90_10_K10.eps 
ss = sprintf('fig_90_10_K10');
saveas(fig,ss,"jpeg")


% figcounter = 2;
% for i = [1:3:20]
%     figure(figcounter);
%     ss = sprintf('./data/cry%d.m',i);
%     C = load(ss);
%     
%     mesh(xx,yy,C);
%     
%     figcounter = figcounter +1;
% end


