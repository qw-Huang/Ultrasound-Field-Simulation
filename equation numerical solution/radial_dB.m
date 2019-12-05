%用rayleigh积分计算轴向声压分布
clc
close all;
clear all;

n=0:0.002:0.5;
% n=0;
a=7.5e-3;
f0=1e6;%定义频率和声功率
medium = set_medium('lossless');%定义介质（单层：水）
lambda = medium.soundspeed/f0;%波长=c/f
u=1;
A=15e-3;
R=-2e-3:0.000001:2e-3;
k=2*pi/lambda;%波数
x=0;
% tic
% for i=1:length(n)
%     hole_a=n(i)*a;
%     y2=asin(a/A);
%     y1=asin(hole_a/A);
%     r=sqrt(w^2-2*A*R*sin(theta)*cos(phi)+R.^2+2*A*(A-w)*(1-cos(phi)));
%     f1=exp(-j*k.*r)./r;
%     f2=sin(theta).*int(f1,0,pi);
%     p=j*medium.density*medium.soundspeed*u*k*A^2/pi.*int(f2,theta,y1,y2);
%     p_abs=abs(p);
%     p_abs=double(p_abs);
%     p_max(i)=max(p_abs(:));
%     index_R=find(p_abs>=max(p_max(i)*0.5));
%     dB(i)=R(max(index_R))-R(min(index_R));
% end
% toc
z=10e-3:0.000001:15e-3;
%调用rayleigh积分计算前移、焦点声压、轴向-6dB
tic
for i=1:length(n)
   hole_a(i)=n(i)*a;
   pr_axial=hole_rayleigh2(A,a,u,hole_a(i),x,z);
   pr_axial=abs(pr_axial);
   pr_max(i)=max(pr_axial(:));
   index_focus_z(i)=find(pr_axial==pr_max(i));
   pr_radial(i,:)=hole_rayleigh2(A,a,u,hole_a(i),R,z(index_focus_z(i)));
   pr_radial(i,:)=abs(pr_radial(i,:));
   index_radial=find(pr_radial(i,:)>=0.5*pr_max(i));
   dB_radial(i)=R(max(index_radial))-R(min(index_radial));
end
toc
scatter(n,dB_radial*1000,5,'.','b');



