clc;
clear all;
P=100;
R=15e-3;
a=7.5e-3;
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

xmin=-a;
xmax=a;
dx =2.5e-4;
x=xmin:dx:xmax;
n=0:2.5:5;
for i = 1:length(n)
    hole_a=n(i)*0.1*a;
    pr_abs=hole_rayleigh2(R,a,P,hole_a);
    plot(x/R,pr_abs(:,37));
    hold on;
end
