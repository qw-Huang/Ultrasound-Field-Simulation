% 使用自定义rayleigh+ASA(单层)计算3D声强分布，针对多层组织  
% 该方法得到的计算结果误差较大，焦点声压幅值误差8%   因此不考虑使用
clc;
clear all;

% Set up the array
R=75e-3;
a=30e-3;

% Use a layered medium
medium1 = set_medium('water');
medium2 = set_medium('water');
% P=100;
% u=normal_velocity(P,R,a,0,medium1.density,medium1.soundspeed);
u=1;
% Center frequency and wavelength
f0 = 1e6;
lambda = medium1.soundspeed/f0;
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

% Set up the coordinate grid
xmin = -1.5*a;
xmax = -xmin;
ymin = -1.5*a;
ymax = -ymin;
zmin = 30e-3;
zmax = 90e-3;

dx = lambda/6;
dy = lambda/6;
dz = lambda/6;

x = xmin:dx:xmax;
y = ymin:dy:ymax;
z = zmin:dz:zmax;

% Determine where the source pressure will be calculated
z0=29e-3;
y_median = floor((ymax-ymin)/2/dy)+1;
x_median=floor(length(x)/2)+1;

zmin1=z0;
zmax1=30e-3;
z1=zmin1:dz:zmax1;

cg_3d1 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin1,zmax1);
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);

tic
pr=rayleigh_2D_xy(R,a,f0,u,medium1,x,y,z0);
[p_asa1,p_interface]=layer_cw_angular_spectrum(pr,cg_3d1,medium1,f0,1024,'Pa');
p_asa=cw_angular_spectrum(p_interface,cg_3d,medium2,f0,1024,'Pa');
toc

p_asa_abs=abs(p_asa);
p_asa_max=max(p_asa_abs(:));
focus_index_1=find(p_asa_abs==p_asa_max);
s=size(p_asa_abs);
[x_index,y_index,z_index]=ind2sub(s,focus_index_1);%将最大值单下标转为三维多下标
focus_index=z_index;
focus_forward=R-z(focus_index);
%轴向-6dB
axial_dB_index=find(p_asa_abs(x_median,y_median,:)>=0.5*p_asa_max);
axial_dB=z(max(axial_dB_index))-z(min(axial_dB_index));
radial_dB_index=find(p_asa_abs(:,y_median,z_index)>=0.5*p_asa_max);
radial_dB=x(max(radial_dB_index))-x(min(radial_dB_index));

%两部分拼接起来
p1=squeeze(p_asa1(:,y_index,:));
p2=p_asa(:,y_index,:);
p2(:,1)=[];
p=[p1 p2];
p_abs=abs(p);
z(1)=[];
z=[z1 z];

figure(3);
surf(z*1000, x*1000, p_abs,C);
xlabel('z (mm)');
ylabel('x (mm)');
axis equal;
shading flat;
title('rayleigh + ASA ');
