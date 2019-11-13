%使用内置rayleigh计算3D声场分布，(画出xy平面声场分布，但是结果存在错误) 
clc;
clear all;
f0=1e6;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f

k=2*pi/lambda;%波数

R = 5 * 2 * lambda;%ROC曲率半径
a = 5 * lambda;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

%划分网格点
xmin=-0.7*a;%观察点坐标的范围
xmax=-xmin;
ymin=xmin;
ymax=xmax;
zdiff=0.7*d;
zmin=R-zdiff;
zmax=R+zdiff;

nx = 61;%网格点的分割点数 
ny = 61; 
nz = 101;

dx = (xmax-xmin)/(nx-1); %网格点的步长
dy = (ymax-ymin)/(ny-1); 
dz = (zmax-zmin)/(nz-1);

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

%调用focus内置rayleigh函数rayleigh_cw()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr = get_spherical_shell(a,R); %得到球面换能器

delta = [dx dy dz];%网格点步长
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%网格点的划分
ndiv =100;%积分的点数
dflag = 0;%如果=1，结果有什么不同？

tic
prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS自带的rayleigh计算
I_prs= acousticintensity(prs,medium.density,medium.soundspeed); %声强计算 
toc


%找到感兴趣的区域声强-6dB范围，画出这个范围内的error;找到最大值，返回对应的三维坐标，在最大点对应的xy平面，把感兴趣区域画出来
max_index=find_maxpoint(I_prs);%返回最大点位置坐标
x_index=max_index(1);%最大值位置x下标
y_index=max_index(2);%最大值位置y下标
z_index=max_index(3);%最大值位置z下标
I_prs_max=max(I_prs(:));%自定义rayleigh声强最大值


% plot the pressures
figure(1);
surf(abs(squeeze(I_prs(:,:,z_index)))/abs(I_prs_max));%对prs归一化   为什么对应的坐标是z在前面
axis equal;
shading interp;
colorbar;
title('Rayleigh Sommerfeld Result(z=R）');
xlabel('x ');
ylabel('y ');
zlabel('normalized pressure');
figure(2);
    axis equal;
    surf(abs(squeeze(I_prs(:,y_index,:)))/abs(I_prs_max));%对prs归一化   为什么对应的坐标是z在前面
    colorbar
    title('Rayleigh Sommerfeld Result(y=0)');
    xlabel('z');
    ylabel('x');
    zlabel('normalized pressure');

