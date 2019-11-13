%用内置FNM函数，求xz平面声压分布，得到平面轴向和径向声压-3dB（或声强-6dB）大小
clc;
clear all;
lossless = set_medium('lossless');
f = 1e6;
lambda = lossless.soundspeed / f;
omega = 2 * pi * f;
k = 2 * pi / lambda;

R = 6 * 1.1*lambda;
a = 6* lambda;
d = sqrt(R^2 - a^2);
phi0 = asin(a/R);
dBr=1.1*lambda;%验证径向-6dB范围=f-number*波长（是否是声强的-6dB）

xdcr = get_spherical_shell(a,R);% 球形换能器
%划分网格点
xmin = -a;
xmax = a;
ymin = 0;
ymax = 0;
zmin = R-d;
zmax = 1.5*d + R;

nz = 201; % ok to sample the origin
nx = 101;

if nx > 1,
    dx = 2 * a / (nx - 1);
else
    dx = 0;
end

if nz > 1,
    dz = 2.5*d/ (nz - 1);
else
    dz = 0;
end
delta = [dx 0 dz];

x = xmin:dx:xmax;
y = ymin:dx:ymax;
z = zmin:dz:zmax;

ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);

ndiv = 200;
dflag = 0;
tic
pref=fnm_call(xdcr,ps,lossless,ndiv,f,dflag);
toc

figure(1)
mesh(z*1000, x*1000, abs(squeeze(pref)))
xlabel('z (mm)');
ylabel('x (mm)');
zlabel('normalized pressure');

p= abs(squeeze(pref));
[Pm1,Ix]=max(p);%找到声压每一列最大值，返回对应行数
[Pm,Iy]=max(Pm1());%找到有最大值声压的一列，返回列数
Ij=Iy;
Ii=Ix(1,Ij);
Pmax=p(Ii,Ij);%找到声压最大值Pmax

pa=p(Ii,:);%x为0的那一行，也就是轴向声压分布
figure(2);
plot(z*1000,pa);
xlabel('z (mm)');
ylabel('axial pressure');

co=find(pa>0.5*Pmax);%在声压最大点的行，选出声压>0.5*Pmax的列
co1=min(co);co2=max(co);%找到声压>0.5*Pmax的最小的列和最大的列数
range_c=(co2-co1)*dz;%有个问题，得到的范围包括了近场超过0.25倍pmax的数值，但是并不是想要的，希望得到连续的-6dB的范围
pr=p(:,Ij);%x为0的那一行，也就是轴向声压分布

figure(3);
plot(x*1000,pr);
xlabel('x (mm)');
ylabel('radial pressure ');

ro=find(pr>0.5*Pmax);%在声压最大点的列，选出声压>0.5*Pmax的行
ro1=min(ro);ro2=max(ro);%找到声压>0.5*Pmax的最小的列和最大的行数
range_r=(ro2-ro1)*dx;
