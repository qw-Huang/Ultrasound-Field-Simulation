clc;
clear all;
close all;

n=0:0.002:0.5;
f0=1e6;%定义频率和声功率
medium = set_medium('lossless');%定义介质（单层：水）
lambda = medium.soundspeed/f0;%波长=c/f
u=1;
k=2*pi/lambda;%波数

% % %根据公式， 轴向-6dB范围变化范围
w=0:0.000001:25e-3;
A=15e-3;
a=5e-3;
for i=1:length(n)
    hole_a=n(i)*a;
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    p_max(i)=max(p_abs(:));
    index_z=find(p_abs>=max(p_max(i)*0.5));
    dB(i)=w(max(index_z))-w(min(index_z));
end
figure(1);
scatter(n,dB*1000,5,'.','b');
hold on;
figure(2);
scatter(n,p_max/10^6,5,'.','b');
hold on


a=7.5e-3;
for i=1:length(n)
    hole_a=n(i)*a;
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    p_max(i)=max(p_abs(:));
    index_z=find(p_abs>=max(p_max(i)*0.5));
    dB(i)=w(max(index_z))-w(min(index_z));
end
figure(1);
scatter(n,dB*1000,5,'.','r');
hold on;
figure(2);
scatter(n,p_max/10^6,5,'.','r');
hold on;


a=10e-3;
for i=1:length(n)
    hole_a=n(i)*a;
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    p_max(i)=max(p_abs(:));
    index_z=find(p_abs>=max(p_max(i)*0.5));
    dB(i)=w(max(index_z))-w(min(index_z));
end
figure(1);
scatter(n,dB*1000,5,'.','k');
xlabel('hole(a)/a');
ylabel('axial distance -6dB(mm)');
figure(2);
scatter(n,p_max/10^6,5,'.','k');
xlabel('hole(a)/a');
ylabel('pmax(MPa)');
