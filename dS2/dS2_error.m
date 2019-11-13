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
phi0=asin(a/R);

%换能器离散化为点声源  注意离散化之后点声源的大小
dR=lambda/6;%每个环带的弧度相等
dphi=dR/R;%xz平面上对应的等角度
nphi=round(phi0/dphi); %环带的个数，
dR=R*phi0/nphi;%nphi取整之后重新调整弧度
dphi=dR/R;%调整之后dphi的值
phi_after=dphi:dphi:phi0;
phi_back=0:dphi:(phi0-dphi);
phi=phi_after-dphi/2;%ri对应的xz平面的角度
ri0=R*sin(phi);%ri是第i环对应半径
Sm=lambda/6;%每个圆环对应的圆弧间隔频率固定
ntheta=round(2*pi*ri0/Sm);  %第i个圆环离散的弧度个数
dtheta=2*pi./ntheta; %根据取整个数，重新调整每一环的dtheta
Sm0=dtheta.*ri0; %调整每一环离散的弧度
Sm=repelem(Sm0,ntheta);
theta=dthetarepet(dtheta,ntheta);%变成一个数组
ri=repelem(ri0,ntheta);
Z0=R-R*cos(phi); %每个半径为ri的圆环对应的z轴坐标
Z=repelem(Z0,ntheta);

%rayleigh积分计算xz平面的声场  
tic
dS=dR.*Sm; % 求点源面积=圆环的采样间隔dR*第i个环的离散的弧度
toc

%比较dS的累加
S_discrete=sum(dS);%离散之后的全部面积累加
S_theory=2*pi*R*(R-d);%理论计算球面面积，即解析面积
dSi=Sm0.*dR.*ntheta;%求每一环面积的累加，也就是把每一列（70个数）累加
ri_back=R*sin(phi_back);%ri是第i环短弧端对应半径
ri_after=R*sin(phi_after);%ri是第i环场=长弧端对应半径
di_back=sqrt(R^2-ri_back.^2);%从0开始作为第i层环的d
di_after=sqrt(R^2-ri_after.^2);%从dr开始作为第i层环的d
dSi_theory=2*pi*R*(R-di_after)-2*pi*R*(R-di_back);%第i环的环带面积
error_dSi=abs(dSi-dSi_theory)./dSi;

%画图
figure(4);
histogram(dS,10);
title('点源面积dS分布(R=150mm ,a=75mm)');

figure(5);
plot(error_dSi);
xlabel('ith ring ');
ylabel('error ');
title('the error of ith ring (R=150mm, a=75mm)');
