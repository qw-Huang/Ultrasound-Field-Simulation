%开圆形孔，比较不同开孔的大小，声压焦点前移、-6dB变化情况
clc;
clear all;
close all;
u=1;
R=2*15e-3;
a=2*7.5e-3;
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

% zDiff=d;
xmin=-0.6*a;
xmax=0.6*a;
zDiff=0.99*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = 2.5e-4; %网格点的步长
dz = 2.5e-4;

x=xmin:dx:xmax;%网格点的分布
z=zmin:dz:zmax;
% % n=0.35
%  n=0:0.1:0.5;
n=[0 0.37 0.72];
for i = 1:length(n)
    hole_a=n(i)*a;
    pr=hole_rayleigh2(R,a,u,hole_a,x,z);
    pr_abs=abs(pr);
    pr_max(i)=max(pr_abs(:));
    
    z_range=find(pr_abs(37,:)>=0.5*pr_max(i));
    rangedB_back=z_range(1);
    rangedB_after=z_range(length(z_range));
    
    x_range=find(pr_abs(:,103)>=0.5*pr_abs(37,103));
    xrangedB_back=x_range(1);
    xrangedB_after=x_range(length(x_range));
    
    figure(1);
    plot(z./R,pr_abs(37,:));
    hold on;
    plot(z(rangedB_back)./R,pr_abs(37,rangedB_back),'b-o');
    hold on;
    plot(z(rangedB_after)./R,pr_abs(37,rangedB_after),'b-o');
    hold on;
    figure(2);
    plot(x./R,pr_abs(:,103));
    hold on;
    plot(x(xrangedB_back)./R,pr_abs(xrangedB_back,103),'b-o');
    hold on;
    plot(x(xrangedB_after)./R,pr_abs(xrangedB_after,103),'b-o');
    hold on;
    
end
    xlabel('z/R');
    ylabel('pressure(Pa)');
%  figure(2);
%  plot(n,pr_max);
