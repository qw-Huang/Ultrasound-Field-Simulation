
%对比理论推导结果和离散化计算结果关于声学特性的变化
clc;
clear all;
close all;
tic 

%对文献中的声压幅值公式进行求导，得到焦点处附近一阶导数等于0的时候得到焦点声压最大值对应的w
syms w;
n=0:0.002:0.5;
a=7.5e-3;
f0=1e6;%定义频率和声功率
medium = set_medium('lossless');%定义介质（单层：水）
lambda = medium.soundspeed/f0;%波长=c/f
u=1;
A=15e-3;
ROC=15e-3;
k=2*pi/lambda;%波数

%根据文献中公式，焦点前移距离随开孔比例的变化
for i=1:length(n)
    hole_a=n(i)*a;
    y2=asin(a/A);
    y1=asin(hole_a/A);
    U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
    U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
    p_abs=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
    dp=diff(p_abs,w); %对函数p(w)求导，在ROC附近一阶导数为0的点是最大值点的w坐标
    w0=vpasolve(dp==0,w,14.9e-3);
    w_forward(i)=A-double(w0); %（ROC-w）得到前移距离    
end
% figure(1);
% scatter(n,w_forward,'.','k');
% xlabel('开孔大小');
% ylabel('前移距离');
% hold on;
toc

x=0;
z=0:0.000001:22e-3;
%调用rayleigh积分计算前移、焦点声压、轴向-6dB
for i=1:length(n)
   hole_a(i)=n(i)*a;
   pr_axial=hole_rayleigh2(ROC,a,u,hole_a(i),x,z);
   pr_axial=abs(pr_axial);
   pr_max(i)=max(pr_axial(:));
   index_focus_z(i)=find(pr_axial==pr_max(i));
   focus_forward(i)=ROC-z(index_focus_z(i));
    
   %轴向-6dB和开孔大小的关系
   index_dB_axial=find(pr_axial>=0.5*pr_max(i));
   index_axial_forward(i)=min(index_dB_axial);
   index_axial_backward(i)=max(index_dB_axial);
   axial_dB(i)=z(index_axial_backward(i))-z(index_axial_forward(i));
end
% figure(1);
% scatter(n,focus_forward,'.','b');
% figure(2);
% scatter(n,pr_max,'.','b');
% hold on;
% figure(3);
% scatter(n,axial_dB,'.','b');
% 

% % %根据公式， 轴向-6dB范围变化范围
w=0:0.000001:22e-3;
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
% figure(3);
% scatter(n,dB,'.','k');
% xlabel('开孔大小');
% ylabel('轴向-6dB范围变化');
% figure(2);
% scatter(n,p_max,'.','k');
% xlabel('开孔大小');
% ylabel('pmax');

%画图
figure(1);
scatter(n,w_forward*1000,5,'.','k');
xlabel('hole(a)/a');
ylabel('focus forward distance(mm)');
hold on;
scatter(n,focus_forward*1000,5,'.','b');

figure(2);
scatter(n,dB*1000,5,'.','k');
xlabel('hole(a)/a');
ylabel('axial distance -6dB(mm)');
hold on;
scatter(n,axial_dB*1000,5,'.','b');

figure(3);
scatter(n,p_max/10^6,5,'.','k');
xlabel('hole(a)/a');
ylabel('pmax(MPa)');
hold on;
scatter(n,pr_max/10^6,5,'.','b');


% 
% f_number=0.51:0.05:5;
% hole_a=0;
% %改变f-number 换能器焦点前移情况
% for i=1:length(f_number)
%     a=A/(2*f_number(i));
%     y2=asin(a/A);
%     y1=asin(hole_a/A);
%     U2=sqrt(w.^2+2*A.*(A-w).*(1-cos(y2)));
%     U1=sqrt(w.^2+2*A.*(A-w).*(1-cos(y1)));
%     p=2*medium.density*medium.soundspeed*u*A./abs(A-w).*abs(sin(k.*(U2-U1)./2));
%     index=find(p==max(p(:)));
%     focus_forward(i)=A-w(index);
% end
% figure(3);
% scatter(f_number,focus_forward,'.','k');
% xlabel('f-number');
% ylabel('前移距离');
% toc