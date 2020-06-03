%ʹ������rayleigh����xy�������ֲ��������������xyƽ��
clc;
clear all;clear all;
f0=1.4e6;%����Ƶ�ʺͷ�������
medium = set_medium('muscle');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����
P=100;

R = 75e-3;%ROC���ʰ뾶
a = 30e-3;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

%���������
xmin=-1.5*a;%�۲������ķ�Χ
xmax=-xmin;
ymax=0;
ymin=0;
zmin=30e-3;
zmax=80e-3;

% nx = 61;%�����ķָ���� 
% ny = 61; 

dx = 2.5e-4; %�����Ĳ���
dy = 1; 
dz = 2.5e-4;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

%����focus����rayleigh����rayleigh_cw()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr=get_spherical_shell(a,R);
%xdcr = create_spherical_shell_planar_array(1, 1, a, R, 0.01, 0.01);

delta = [dx dy dz];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%�����Ļ���
ndiv = 200;%���ֵĵ���
dflag = 0;

tic
prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0);%FOCUS�Դ���rayleigh����xyƽ��
toc
prs=prs*u;
prs_normalized=abs(prs./max(prs(:)));%��prs��һ��
% plot the pressures
figure(1);
surf(z*1000,x*1000,squeeze(prs_normalized));   
shading interp;
axis equal;
colorbar;
title('Rayleigh Sommerfeld Result(z=R)');
xlabel('y');
ylabel('x');

