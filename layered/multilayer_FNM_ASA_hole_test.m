%使用内置FNM和ASA，计算两层组织声强分布，考虑界面的反射和折射

clear all;
clear all;
%set up the tranducer

R = 120e-3;
fnumber = 0.8;
a=R/(2*fnumber);
d = sqrt(R^2 - a^2);
P=300;
hole_a=2e-2;
% Use a layered medium
medium1 = set_medium('water');%
medium2 = set_medium('skin');
medium3 = set_medium('fat');
medium4 = set_medium('muscle');
medium5 = set_medium('muscle');
% Center frequency and wavelength
f0=1e6;
lambda = medium1.soundspeed/f0;
u=normal_velocity(P,R,a,hole_a,medium1.density,medium1.soundspeed);
depth=0.08;
z0=R-d+lambda;
z_interface1=R-depth;
z_interface2=z_interface1+2e-3;
z_interface3=z_interface2+18e-3;
z_interface4=z_interface3+40e-3;
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
zmax5=R+20e-3;
z5=zmin5:dz:zmax5;

nx=length(x);
ny=length(y);
nz1=length(z1);
nz2=length(z2);
nz3=length(z3);
nz4=length(z4);
nz5=length(z5);

% Create planar array of rectangular transducers
xdcr_array = get_spherical_shell(a,R);% 球形换能器
xdcr_array_hole = get_spherical_shell(hole_a,R);

% Determine where the source pressure will be calculated
y_median = floor((ymax-ymin)/2/dy)+1; %floor 的使用有点问题
x_median=floor(length(x)/2)+1;

% pressure (x-z plane)
cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
cg_3d1 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin1,zmax1);
cg_3d2 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin2,zmax2);
cg_3d3 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin3,zmax3);
cg_3d4 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin4,zmax4);
cg_3d5 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin5,zmax5);
% Calculate the pressure
ndiv = 200;
    tic();
    p0_all = fnm_call(xdcr_array,cg_p0,medium1,ndiv,f0,0);
    p0_hole = fnm_call(xdcr_array_hole,cg_p0,medium1,ndiv,f0,0);
    p0_res=p0_all-p0_hole;
    p0=p0_res.*u;
[p_asa1,p_interface1]=layer_cw_angular_spectrum(p0,cg_3d1,medium1,medium2,f0,1024,'Pa');
[p_asa2,p_interface2]=layer_cw_angular_spectrum(p_interface1,cg_3d2,medium2,medium3,f0,1024,'Pa');
[p_asa3,p_interface3]=layer_cw_angular_spectrum(p_interface2,cg_3d3,medium3,medium4,f0,1024,'Pa');
[p_asa4,p_interface4]=layer_cw_angular_spectrum(p_interface3,cg_3d4,medium4,medium5,f0,1024,'Pa');
p_asa5=cw_angular_spectrum(p_interface4,cg_3d5,medium5,f0,1024,'Pa');
toc
p_asa_abs=abs(p_asa5);
p_asa_max=max(p_asa_abs(:));
focus_index_1=find(p_asa_abs==p_asa_max);
s=size(p_asa_abs);
[x_index,y_index,z_index]=ind2sub(s,focus_index_1);%将最大值单下标转为三维多下标
%轴向-6dB
% axial_dB_index=find(p_asa_abs(x_median,y_median,:)>=0.5*p_asa_max);
% axial_dB=z2(max(axial_dB_index))-z2(min(axial_dB_index));
% radial_dB_index=find(p_asa_abs(:,y_median,z_index)>=0.5*p_asa_max);
% radial_dB=x(max(radial_dB_index))-x(min(radial_dB_index));
%两部分拼接起来
p1=squeeze(p_asa1(:,y_index,:));
p2=squeeze(p_asa2(:,y_index,:));
p3=squeeze(p_asa3(:,y_index,:));
p4=squeeze(p_asa4(:,y_index,:));
p5=squeeze(p_asa5(:,y_index,:));
p2(:,1)=[];
p3(:,1)=[];
p4(:,1)=[];p5(:,1)=[];
p=[p1 p2 p3 p4 p5];
p_abs=abs(p);
z2(1)=[];
z3(1)=[];z4(1)=[];z5(1)=[];
z=[z1 z2 z3 z4 z5];

figure(3);
surf(z*1000, x*1000, p_abs);
xlabel('z (mm)');
ylabel('x (mm)');
axis equal;
shading flat;
title('FNM + ASA ');
