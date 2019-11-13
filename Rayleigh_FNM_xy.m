clc;
clear all;
%比较同一个球面换能器FOCUS中使用FNM函数和rayleigh积分函数的声压分布的区别,prs-pref,再取模（在边缘有明显的差别，其他位置差别不大）
%波长=1.5mm，R=15mm f-number=1    x方向分辨率：0.25mm   z方向分辨率：0.2598mm  离散化点声源大小（点声源个数ndiv=200）
%dflag什么含义？一般设为0，改成1试试
lossless = set_medium('lossless');

f0 = 1e6;
lambda = lossless.soundspeed / f0;
omega = 2 * pi * f0;
k = 2 * pi / lambda;

R =  5 * 2 * lambda;
a =  5 * lambda;
d = sqrt(R^2 - a^2);
phi0 = asin(a/R);

xdcr = get_spherical_shell(a,R);
xmin = -a;
xmax = a;
ymin = -a;
ymax = a;
z0=R;
zmin=z0;
zmax=z0;

nx = 121;
ny = 121;

dx= 2 * a / (nx - 1);
dy= 2 * a / (ny - 1);
delta = [dx dy 0];

x = xmin:dx:xmax;
y= ymin:dy:ymax;

ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);

ndiv = 100;
dflag = 0;
tic
pref=fnm_call(xdcr,ps,lossless,ndiv,f0,dflag);
toc

tic
prs=rayleigh_cw(xdcr,ps,lossless,ndiv,f0,dflag);
RStimesphere= toc;

% plot the pressures
figure(1);
pref_normalized=abs(pref)./abs(max(pref(:)));
surf(pref_normalized);
shading interp
colorbar;
title('Fast Nearfield Method Result');
xlabel('y');
ylabel('x ');
zlabel('normalized pressure');

figure(2);
prs_normalized=abs(prs)./abs(max(prs(:)));
surf(prs_normalized);
shading interp
colorbar;
title('Rayleigh Sommerfeld Result');
xlabel('y ');
ylabel('x');
zlabel('normalized pressure');

figure(3);
plot(prs_normalized(30,:));
title('x=0 Rayleigh Sommerfeld Result');
xlabel('y');
ylabel('normalized pressure');

figure(4);
plot(prs_normalized(:,30));
title('y=0 Rayleigh Sommerfeld Result');
xlabel('x');
ylabel('normalized pressure');

error=abs(squeeze(pref - prs))./abs(squeeze(pref));%结果小于1%
figure(5);
surf(error);
shading interp
colorbar;
title('error of between rayleigh and FNM');
xlabel('y ');
ylabel('x');
zlabel('pressure');

figure(6);
error0=abs(pref/abs(max(max(pref)))-prs/abs(max(max(prs))))./abs(pref/abs(max(max(pref)))); %结果小于1%
surf(squeeze(error0));
shading interp
