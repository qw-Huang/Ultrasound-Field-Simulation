%使用内置FNM（FOCUS）计算3D空间声场分布
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
ymax=0.7*a;
ymin=-ymax;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %网格点的步长
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);
ny=length(y);
nz=length(z);

%调用focus内置FNM函数fnm_call()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr = get_spherical_shell(a,R); %得到球面换能器

if nz > 1,
    dz = 2 * d / (nz - 1);
else
    dz = 0;
end
if nx > 1,
    dx = 2 * a / (nx - 1);
else
    dx = 0;
end
if ny > 1,
    dy = 2 * a / (ny - 1);
else
    dy = 0;
end

delta = [dx dy dz];%网格点步长
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%网格点的划分
ndiv = 60;%积分的点数
dflag = 0;%如果=1，结果有什么不同？

tic
pref=fnm_call(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS自带的rayleigh计算
toc

I=pref.*pref/(medium.soundspeed*medium.density);

z_index=51;
z_index_value=(z_index-1)*dz;

% plot the pressures
if nx > 1 & nz > 1 & ny > 1,
    figure(1);
    axis equal;
    If_normalized=abs(squeeze(I(:,:,z_index)))/abs(max(max(squeeze(I(:,:,z_index)))));
    mesh(If_normalized);%对prs归一化   为什么对应的坐标是z在前面？
    If_xy=squeeze(I(:,:,z_index))/abs(max(max(max(I))));%对内置rayleigh的prs归一化处理，除以三维空间中最大的点，得到复数
    colorbar
    title('Rayleigh Sommerfeld Result');
    xlabel('nx ');
    ylabel('ny ');
    zlabel('normalized intensity');
end
