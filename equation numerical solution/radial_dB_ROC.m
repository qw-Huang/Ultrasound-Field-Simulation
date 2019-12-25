%��rayleigh����
clc
close all;
clear all;

n=0:0.002:0.5;
% n=0;
f0=1e6;%����Ƶ�ʺ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ��
lambda = medium.soundspeed/f0;%����=c/f
u=1;
k=2*pi/lambda;%����
x=0;

A=15e-3;
a=7.5e-3;
R=-2e-3:0.000001:2e-3;
z=10e-3:0.000001:15e-3;
%����rayleigh���ּ���ǰ�ơ�������ѹ������-6dB
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

A=75e-3;
a=37.5e-3;
R=-1.5e-3:0.000001:1.5e-3;
z=73e-3:0.000001:75e-3;
%����rayleigh���ּ���ǰ�ơ�������ѹ������-6dB
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
scatter(n,dB_radial*1000,5,'.','r');

A=150e-3;
a=75e-3;
R=-1.5e-3:0.000001:1.5e-3;
z=143e-3:0.000001:150e-3;
%����rayleigh���ּ���ǰ�ơ�������ѹ������-6dB
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
scatter(n,dB_radial*1000,5,'.','k');