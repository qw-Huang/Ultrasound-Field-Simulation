clc;
clear all;
close all;
tic 

%对文献中的声压幅值公式进行求导，得到焦点处附近一阶导数等于0的时候得到焦点声压最大值对应的w
syms w;
n=0:0.002:0.5;
a=7.5e-3;
f0=1e6;%定义频率和声功率
medium = set_medium('lossless');%定义介质（单层：水）
lambda = medium.soundspeed/f0;%波长=c/f
P=100;
A=15e-3;
k=2*pi/lambda;%波数

%根据文献中公式，焦点前移距离随开孔比例的变化
for i=1:length(n)
    hole_a=n(i)*a;
        u=normal_velocity(P,A,a,hole_a,medium.density,medium.soundspeed);
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    dp=diff(p_abs,w); %对函数p(w)求导，在ROC附近一阶导数为0的点是最大值点的w坐标
    w0=vpasolve(dp==0,w,14.9e-3);
    w_forward(i)=A-double(w0); %（ROC-w）得到前移距离    
end
figure(1);
scatter(n,w_forward,'.','k');
xlabel('开孔大小');
ylabel('前移距离');
hold on;
toc