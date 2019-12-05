%对文献中的声压幅值公式进行求导，得到焦点处附近一阶导数等于0的时候得到焦点声压最大值对应的w
clear all;
syms w;
% w=0:2.5e-4:14e-3;
f0=1e6;%定义频率和声功率
medium = set_medium('lossless');%定义介质（单层：水）
lambda = medium.soundspeed/f0;%波长=c/f
u=1;
A=15e-3;
k=2*pi/lambda;%波数
y2=pi/5;
% y1=0.35*y2;
y1=0;
n=sin(y1)/sin(y2);
U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
p=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
dp=diff(p,w);
w0=vpasolve(dp==0,w,14e-3);
w1=subs(w0);
% w=0:2.5e-4:14e-3;

