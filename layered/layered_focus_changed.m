%两层组织 用改进的rayleigh积分  有问题程序会崩MATLAB自动关闭
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
% R=-1.5e-3:0.00001:1.5e-3;
z=8e-3:0.00001:22e-3;
dB_axial=zeros(1,length(n));
hole_a=zeros(1,length(n));
pr_max=zeros(1,length(n));
index_focus_z=zeros(1,length(n));
focus_forward=zeros(1,length(n));
%调用rayleigh积分计算前移、焦点声压、轴向-6dB
tic
for i=1:length(n)
   hole_a(i)=n(i)*a;
   pr_axial=layered_hole_rayleigh_2D(A,a,hole_a(i),P,medium,0,z);
   pr_axial=abs(pr_axial);
   pr_max(i)=max(pr_axial(:));
   index_axial=find(pr_axial>=0.5*pr_max(i));
   index_focus_z(i)=find(pr_axial==pr_max(i));
    dB_axial(i)=z(max(index_axial))-z(min(index_axial));
    focus_forward(i)=A-z(index_focus_z(i));
end
toc
figure(1);
scatter(n,dB_axial*1000,25,'.','k');
figure(2);
scatter(n,focus_forward*1000,25,'.','k');
figure(3)
scatter(n,pr_max/10^6,25,'.','k');
% for i=1:length(n)
% 
%    pr_radial(i,:)=layered_hole_rayleigh_2D(A,a,hole_a(i),P,medium,R,z(index_focus_z(i)));
%    pr_radial(i,:)=abs(pr_radial(i,:));
%    index_radial=find(pr_radial(i,:)>=0.5*pr_max(i));
%    dB_radial(i)=R(max(index_radial))-R(min(index_radial));
% end
% figure(2);
% scatter(n,dB_radial*1000,5,'.','b');