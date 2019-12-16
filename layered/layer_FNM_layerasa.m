%使用内置FNM和ASA，计算两层组织声强分布，没考虑界面的反射和折射

clear all;
clear all;
%set up the tranducer

R = 75e-3;
a = 30e-3;
fnumber=R/(2*a);
d = sqrt(R^2 - a^2);

% Use a layered medium
medium = set_layered_medium([0,30e-3],[set_medium('water'),set_medium('skin')]);


% Center frequency and wavelength
f0=1e6;
lambda = medium(1).soundspeed/f0;


%划分网格点
xmin=-1.5*a;
xmax=1.5*a;
ymin=-1.5*a;
ymax=1.5*a;
zmin=20e-3;
zmax=80e-3;
dx=lambda/6;
dy=lambda/6;
dz=lambda/6;
x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;
nx=length(x);
ny=length(y);
nz=length(z);

% Create planar array of rectangular transducers
xdcr_array = get_spherical_shell(a,R);% 球形换能器

% Determine where the source pressure will be calculated
z0 = 30e-3;
y_median = floor((ymax-ymin)/2/dy)+1;
x_median=floor(length(x)/2)+1;

% pressure (x-z plane)
cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
% cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);
% Calculate the pressure
ndiv = 200;
tic();
p0 = fnm_call(xdcr_array,cg_p0,medium,ndiv,f0,0);%'fnm sse'

tic();
p_asa = layerasa(p0,z,medium,1024,dz,f0);
% p_asa=layerasa(p0,z,medium,f0,1024,'Pa');
p_asa_abs=abs(p_asa);
p_asa_max=max(p_asa_abs(:));

focus_index_1=find(p_asa_abs==p_asa_max);
s=size(p_asa_abs);
[x_index,y_index,z_index]=ind2sub(s,focus_index_1);%将最大值单下标转为三维多下标
%轴向-6dB
axial_dB_index=find(p_asa_abs(x_median,y_median,:)>=0.5*p_asa_max);
axial_dB=z(max(axial_dB_index))-z(min(axial_dB_index));

radial_dB_index=find(p_asa_abs(:,y_median,z_index)>=0.5*p_asa_max);
radial_dB=x(max(radial_dB_index))-x(min(radial_dB_index));

figure(2);
p=abs(squeeze(p_asa(:,y_index,:)));
surf(z*1000, x*1000, p);
xlabel('z (mm)');
ylabel('x (mm)');
shading flat;
title('normalized ASA Pressure (y=0)');
