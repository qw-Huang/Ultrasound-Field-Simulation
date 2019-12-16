%用rayleigh积分计算轴向声压分布
clc
% close all;
clear all;

n=0:0.002:0.5;
% n=0;
a=7.5e-3;
f0=1e6;%定义频率和声功率
medium = set_medium('lossless');%定义介质（单层：水）
lambda = medium.soundspeed/f0;%波长=c/f
P=100;
A=15e-3;
R=-2e-3:0.000001:2e-3;
k=2*pi/lambda;%波数
x=0;

z=10e-3:0.000001:15e-3;
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
figure(4);
scatter(n,dB_radial*1000,5,'.','k');
xlabel('hole(a)/a');
ylabel('radial distance -6dB(mm)');

save('radial_-6dB_simulation_R15_a7.5.mat','dB_radial');

