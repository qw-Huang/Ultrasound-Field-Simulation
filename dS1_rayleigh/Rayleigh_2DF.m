%ʹ������rayleigh����xy�������ֲ��������������xyƽ��
clc;
clear all;clear all;
f0=1e6;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����

R = 5 * 2* lambda;%ROC���ʰ뾶
a = 5 * lambda;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���

%���������
xmin=-a;%�۲������ķ�Χ
xmax=-xmin;
ymax=a;
ymin=-ymax;
z0=R;

nx = 61;%�����ķָ���� 
ny = 61; 

dx = (xmax-xmin)/(nx-1); %�����Ĳ���
dy = (ymax-ymin)/(ny-1); 
dz = 1;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=z0;

%����focus����rayleigh����rayleigh_cw()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr=get_spherical_shell(a,R);
%xdcr = create_spherical_shell_planar_array(1, 1, a, R, 0.01, 0.01);

delta = [dx dy 0];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, z0, z0);%�����Ļ���
ndiv = 100;%���ֵĵ���
dflag = 0;

tic
prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0);%FOCUS�Դ���rayleigh����xyƽ��
toc
prs_normalized=abs(prs./max(prs(:)));%��prs��һ��
% plot the pressures
figure(1);
surf(prs_normalized);   
shading interp;
axis equal;
colorbar;
title('Rayleigh Sommerfeld Result(z=R)');
xlabel('y');
ylabel('x');

