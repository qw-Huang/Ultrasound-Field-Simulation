%使用dS1点声源划分方法，在2d平面xz面，归一化对比自定义rayleigh积分和focus自带
%的rayleigh积分，对比声强分布的误差
clc;
clear all;
f0=1e6;u=1;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数

R = 10 * 5 * 2 * lambda;%ROC曲率半径
a = 10 * 5 * lambda;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

%换能器离散化为点声源  注意离散化之后点声源的大小
dr=lambda/6;%球面分割成很多个圆环，对应的半径是r，dr是半径增加的步长
r_back=0:dr:a-dr;%点声源前一段弧长对应的r
r_after=dr:dr:a;%点声源后一段弧长对应的r
r=r_after-dr/2;%第i个环带对应的中心点的r
Sm=lambda/6;%最外层环带（2pi*a）离散化后对应的dS的弧长小于lambda/6
median=length(r)/2+1;
ntheta=round(2*pi*r(median)/Sm);%最外层一个环带的划分点数,取整
dtheta=2*pi./ntheta;%根据取整重新调整每个点声源的对应弧度
theta_after=dtheta:dtheta:2*pi;%每个环带离散成多个点对应的弧度数组
theta_back=0:dtheta:(2*pi-dtheta);%每个环带离散成多个点对应的弧度数组
theta=theta_after-dtheta/2;%第i个环带离散成点声源，dS中间的点对应的弧度数组
X=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
Y=sin(theta)'*r;%点源的y坐标
Z0=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。

%计算划分之后的dS
tic
dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
toc

%比较dS的累加和解析面积
S_discrete=sum(sum(dS));%离散之后的全部面积累加
S_theory=2*pi*R*(R-d);%理论计算球面面积，即解析面积
dSi=sum(dS);%求每一环面积的累加，也就是把每一列（70个数）累加
di_back=sqrt(R^2-r_back.^2);%从0开始作为第i层环的d
di_after=sqrt(R^2-r_after.^2);%从dr开始作为第i层环的d
dSi_theory=2*pi*R*(R-di_after)-2*pi*R*(R-di_back);%第i环的环带面积
error_dSi=abs(dSi-dSi_theory)./dSi;%求解析解环带面积和数值解环带面积的误差

%画图
figure(4);
histogram(dS_ring,10);
title('点源面积dS分布(R=150mm ,a=75mm)');

figure(5);
plot(error_dSi);
xlabel('ith ring ');
ylabel('error ');
title('the error of ith ring (R=150mm, a=75mm)');
