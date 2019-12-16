%用rayleigh积分
% clc
close all;
clear all;

n=0:0.005:0.5;
% n=0;
f0=1e6;%定义频率和声功率
medium = set_medium('lossless');%定义介质（单层：水）
lambda = medium.soundspeed/f0;%波长=c/f
P=100;
k=2*pi/lambda;%波数
x=0;

A=75e-3;
% a=30e-3;
R=-2e-3:0.000001:2e-3;
z=70e-3:0.000005:75e-3;
%调用rayleigh积分计算前移、焦点声压、轴向-6dB
% tic
% for i=1:length(n)
%    hole_a(i)=n(i)*a;
%            u(i)=normal_velocity(P,A,a,hole_a(i),medium.density,medium.soundspeed);
%    pr_axial=hole_rayleigh2(A,a,u(i),hole_a(i),x,z);
%    pr_axial=abs(pr_axial);
%    pr_max(i)=max(pr_axial(:));
%    index_focus_z(i)=find(pr_axial==pr_max(i));
%    pr_radial(i,:)=hole_rayleigh2(A,a,u(i),hole_a(i),R,z(index_focus_z(i)));
%    pr_radial(i,:)=abs(pr_radial(i,:));
%    index_radial=find(pr_radial(i,:)>=0.5*pr_max(i));
%    dB_radial(i)=R(max(index_radial))-R(min(index_radial));
% end
% toc
% scatter(n,dB_radial*1000,5,'.','b');
% hold on;
 
% a=37.5e-3;
% R=-2e-3:0.000001:2e-3;
% z=70e-3:0.000005:75e-3;
% %调用rayleigh积分计算前移、焦点声压、轴向-6dB
% tic
% for i=1:length(n)
%    hole_a(i)=n(i)*a;
%            u(i)=normal_velocity(P,A,a,hole_a(i),medium.density,medium.soundspeed);
%    pr_axial=hole_rayleigh2(A,a,u(i),hole_a(i),x,z);
%    pr_axial=abs(pr_axial);
%    pr_max(i)=max(pr_axial(:));
%    index_focus_z(i)=find(pr_axial==pr_max(i));
%    pr_radial(i,:)=hole_rayleigh2(A,a,u(i),hole_a(i),R,z(index_focus_z(i)));
%    pr_radial(i,:)=abs(pr_radial(i,:));
%    index_radial=find(pr_radial(i,:)>=0.5*pr_max(i));
%    dB_radial(i)=R(max(index_radial))-R(min(index_radial));
% end
% toc
% scatter(n,dB_radial*1000,5,'.','r');
% hold on;

a=45e-3;
R=-2e-3:0.000001:2e-3;
z=70e-3:0.000005:75e-3;
%调用rayleigh积分计算前移、焦点声压、轴向-6dB
tic
for i=1:length(n)
   hole_a(i)=n(i)*a;
           u(i)=normal_velocity(P,A,a,hole_a(i),medium.density,medium.soundspeed);
   pr_axial=hole_rayleigh2(A,a,u(i),hole_a(i),x,z);
   pr_axial=abs(pr_axial);
   pr_max(i)=max(pr_axial(:));
   index_focus_z(i)=find(pr_axial==pr_max(i));
   pr_radial(i,:)=hole_rayleigh2(A,a,u(i),hole_a(i),R,z(index_focus_z(i)));
   pr_radial(i,:)=abs(pr_radial(i,:));
   index_radial=find(pr_radial(i,:)>=0.5*pr_max(i));
   dB_radial(i)=R(max(index_radial))-R(min(index_radial));
end
toc
scatter(n,dB_radial*1000,5,'.','k');
xlabel('hole(a)/a');
ylabel('radial distance -6dB(mm)');