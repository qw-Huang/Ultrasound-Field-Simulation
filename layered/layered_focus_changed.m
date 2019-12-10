%两层组织 用改进的rayleigh积分
clc
close all;
clear all;

n=0:0.05:0.5;
% n=0;
f0=1e6;%定义频率和声功率
medium = set_layered_medium([0,5e-3],[set_medium('water'),set_medium('muscle')]);
P=100;
x=0;

A=15e-3;
a=7.5e-3;
R=-2e-3:0.00001:2e-3;
z=8e-3:0.00001:21e-3;
dB_axial=zeros(1,length(n));
%调用rayleigh积分计算前移、焦点声压、轴向-6dB
tic
for i=1:length(n)
   hole_a(i)=n(i)*a;
   pr_axial=layered_hole_rayleigh_2D(A,a,hole_a(i),P,medium,0,z);
   pr_axial=abs(pr_axial);
   pr_max(i)=max(pr_axial(:));
   index_axial=find(pr_axial>=0.5*pr_max(i));
    dB_axial(i)=z(max(index_axial))-z(min(index_axial));
    
   index_focus_z(i)=find(pr_axial==pr_max(i));
   pr_radial(i,:)=layered_hole_rayleigh_2D(A,a,hole_a(i),P,medium,R,z(index_focus_z(i)));
   pr_radial(i,:)=abs(pr_radial(i,:));
   index_radial=find(pr_radial(i,:)>=0.5*pr_max(i));
   dB_radial(i)=R(max(index_radial))-R(min(index_radial));
end
toc
figure(1);
scatter(n,dB_axial*1000,5,'.','b');
figure(2);
scatter(n,dB_radial*1000,5,'.','b');