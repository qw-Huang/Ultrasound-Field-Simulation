%使用内置rayleigh计算xy面声场分布，结果存在问题xy平面
clc;
clear all;clear all;
f0=1.4e6;%定义频率和法向阵速
medium = set_medium('muscle');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数
P=100;

R = 75e-3;%ROC曲率半径
a = 30e-3;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

%划分网格点
xmin=-1.5*a;%观察点坐标的范围
xmax=-xmin;
ymax=0;
ymin=0;
zmin=30e-3;
zmax=80e-3;

% nx = 61;%网格点的分割点数 
% ny = 61; 

dx = 2.5e-4; %网格点的步长
dy = 1; 
dz = 2.5e-4;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

%调用focus内置rayleigh函数rayleigh_cw()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr=get_spherical_shell(a,R);
%xdcr = create_spherical_shell_planar_array(1, 1, a, R, 0.01, 0.01);

delta = [dx dy dz];%网格点步长
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%网格点的划分
ndiv = 200;%积分的点数
dflag = 0;

tic
prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0);%FOCUS自带的rayleigh计算xy平面
toc
prs=prs*u;
prs_normalized=abs(prs./max(prs(:)));%对prs归一化
% plot the pressures
figure(1);
surf(z*1000,x*1000,squeeze(prs_normalized));   
shading interp;
axis equal;
colorbar;
title('Rayleigh Sommerfeld Result(z=R)');
xlabel('y');
ylabel('x');

