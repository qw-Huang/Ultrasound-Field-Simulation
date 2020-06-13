%使用内置FNM和ASA，计算两层组织声强分布，考虑界面的反射和折射

clear all;
clear all;
%set up the tranducer
R = 120e-3;
fnumber = 0.8;
a=R/(2*fnumber);
d = sqrt(R^2 - a^2);
P=100;
hole_a=0;

% Use a layered medium
medium1 = set_medium('water');%
medium2 = set_medium('skin');
medium3 = set_medium('fat');
medium4 = set_medium('muscle');
medium5 = set_medium('blood');
medium6 = set_medium('muscle');

% Center frequency and wavelength
f0=1e6;
lambda = medium1.soundspeed/f0;
u=normal_velocity(P,R,a,hole_a,medium1.density,medium1.soundspeed);
depth=0.055; %组织深度
z0=R-d+lambda;
z_interface1=R-depth;
z_interface2=z_interface1+2e-3;
z_interface3=z_interface2+15e-3;
z_interface4=z_interface3+38e-3;
z_interface5=z_interface4+5e-3;
%划分网格点
dx=lambda/6;
dy=lambda/6;
dz=lambda/6;
xmin=-1.1*a;
xmax=-xmin;
ymin=xmin;
ymax=-ymin;
x=xmin:dx:xmax;
y=ymin:dy:ymax;

zmin1=z0;
zmax1=z_interface1;
z1=zmin1:dz:zmax1;

zmin2=z_interface1;
zmax2=z_interface2;
z2=zmin2:dz:zmax2;

zmin3=z_interface2;
zmax3=z_interface3;
z3=zmin3:dz:zmax3;

zmin4=z_interface3;
zmax4=z_interface4;
z4=zmin4:dz:zmax4;
zmin5=z_interface4;
zmax5=z_interface5;
z5=zmin5:dz:zmax5;
zmin6=z_interface5;
zmax6=R+20e-3;
z6=zmin6:dz:zmax6;

nx=length(x);
ny=length(y);
nz1=length(z1);
nz2=length(z2);
nz3=length(z3);
nz4=length(z4);
nz5=length(z5);
nz6=length(z6);

% Create planar array of rectangular transducers
xdcr_array = get_spherical_shell(a,R);% 球形换能器
% xdcr_array_hole = get_spherical_shell(hole_a,R);

% Determine where the source pressure will be calculated
y_median = floor(length(y)/2)+1; %floor 的使用有点问题
x_median=floor(length(x)/2)+1;

% pressure (x-z plane)
cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
cg_3d1 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin1,zmax1);
cg_3d2 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin2,zmax2);
cg_3d3 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin3,zmax3);
cg_3d4 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin4,zmax4);
cg_3d5 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin5,zmax5);
cg_3d6 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin6,zmax6);
% Calculate the pressure
ndiv = 200;
    tic();
    p0_all = fnm_call(xdcr_array,cg_p0,medium1,ndiv,f0,0);
%     p0_hole = fnm_call(xdcr_array_hole,cg_p0,medium1,ndiv,f0,0);
%     p0_res=p0_all-p0_hole;
%     p0=p0_res.*u;
    p0=p0_all.*u;
[p_asa1,p_interface1]=layer_cw_angular_spectrum(p0,cg_3d1,medium1,medium2,f0,1024,'Pa');
[p_asa2,p_interface2]=layer_cw_angular_spectrum(p_interface1,cg_3d2,medium2,medium3,f0,1024,'Pa');
[p_asa3,p_interface3]=layer_cw_angular_spectrum(p_interface2,cg_3d3,medium3,medium4,f0,1024,'Pa');
[p_asa4,p_interface4]=layer_cw_angular_spectrum(p_interface3,cg_3d4,medium4,medium5,f0,1024,'Pa');
[p_asa5,p_interface5]=layer_cw_angular_spectrum(p_interface4,cg_3d5,medium5,medium6,f0,1024,'Pa');
p_asa6=cw_angular_spectrum(p_interface5,cg_3d6,medium6,f0,1024,'Pa');
toc

%组织部分拼接起来，三维拼接
p_asa3(:,:,1)=[];
p_asa4(:,:,1)=[];
p_asa5(:,:,1)=[];
p_asa6(:,:,1)=[];
p=cat(3,p_asa2, p_asa3, p_asa4, p_asa5,p_asa6);%三维数组串联
p_abs=abs(p);
z3(1)=[];z4(1)=[];z5(1)=[];z6(1)=[];
z=[z2 z3 z4 z5 z6];

p_asa_max=max(p_abs(:));
focus_index_1=find(p_abs==p_asa_max);
s=size(p_abs);
[x_index,y_index,z_index]=ind2sub(s,focus_index_1);%将最大值单下标转为三维多下标
%轴向-6dB
axial_dB_index=find(p_abs(x_index,y_index,:)>=0.5*p_asa_max);
axial_dB=z(max(axial_dB_index))-z(min(axial_dB_index));
radial_dB_index=find(p_abs(:,y_index,z_index)>=0.5*p_asa_max);
radial_dB=x(max(radial_dB_index))-x(min(radial_dB_index));

figure();
surf(z*1000, x*1000, squeeze(p_abs(:,y_index,:)));
xlabel('z (mm)');
ylabel('x (mm)');
axis equal;
shading interp;
title('vascular');
