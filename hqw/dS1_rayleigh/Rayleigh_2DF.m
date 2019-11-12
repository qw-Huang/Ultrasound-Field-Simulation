%使用内置rayleigh计算xy面声场分布，结果存在问题xy平面
clc;
clear all;clear all;
f0=1e6;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数

R = 5 * 2* lambda;%ROC曲率半径
a = 5 * lambda;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

%划分网格点
xmin=-a;%观察点坐标的范围
xmax=-xmin;
ymax=a;
ymin=-ymax;
z0=R;

nx = 61;%网格点的分割点数 
ny = 61; 

dx = (xmax-xmin)/(nx-1); %网格点的步长
dy = (ymax-ymin)/(ny-1); 
dz = 1;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=z0;

%调用focus内置rayleigh函数rayleigh_cw()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr=get_spherical_shell(a,R);
%xdcr = create_spherical_shell_planar_array(1, 1, a, R, 0.01, 0.01);

delta = [dx dy 0];%网格点步长
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, z0, z0);%网格点的划分
ndiv = 100;%积分的点数
dflag = 0;

tic
prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0);%FOCUS自带的rayleigh计算xy平面
toc
prs_normalized=abs(prs./max(prs(:)));%对prs归一化
% plot the pressures
figure(1);
surf(prs_normalized);   
shading interp;
axis equal;
colorbar;
title('Rayleigh Sommerfeld Result(z=R)');
xlabel('y');
ylabel('x');

