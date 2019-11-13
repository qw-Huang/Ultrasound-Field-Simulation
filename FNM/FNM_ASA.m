%内置FNM(FOCUS)计算一个平面，然后用ASA计算三维空间声场分布
clc;clear all;
f0=1e6;u=1;%定义频率和法向阵速
a = 5 * 1.5e-3;%注意这里a是孔径的一半
R = 2*a;%ROC曲率半径
medium = set_medium('lossless');%定义介质：单层->水
lambda = medium.soundspeed/f0;%波长=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%波数
fnumber=1;
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
%划分网格点
xmin=-0.7*a;%观察点坐标的范围
xmax=-xmin;
ymax=0.7*a;
ymin=-ymax;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %网格点的步长
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);
ny=length(y);
nz=length(z);

% Create planar array of rectangular transducers
xdcr_array = get_spherical_shell(a,R);% 球形换能器

% Determine where the source pressure will be calculated
z0 = zmin;
y_index = floor((ymax-ymin)/2/dy);
% Coordinate grids to calclate the initial pressure (x-y plane) and final
% pressure (x-z plane)
cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);

% Calculate the pressure
ndiv = 100;
dflag=0;
fprintf('Calculating initial pressure plane with FNM... ');
tic();
p0 = fnm_call(xdcr_array,cg_p0,medium,ndiv,f0,dflag);
fprintf('done in %f s.\n', toc());

fprintf('Calculating 3D pressure (%i points) with ASA... ', (length(x) * length(y) * length(z)));
tic();
p_asa = cw_angular_spectrum(p0,cg_3d,medium,f0,1024,'Pa');
fprintf('done in %f s.\n', toc());

% Show the initial pressure
figure(1);
mesh(x*1000, y*1000, rot90(abs(squeeze(p0(:,:,1)))));
xlabel('x (mm)');
ylabel('y (mm)');
shading flat;
title(['p0 (Calculated with FNM at z = ', num2str(z0*1000), ' mm)']);
% Show the 3D field calculated with ASA
figure(2);
p_normalized=abs(squeeze(p_asa(:,y_index,:)))/max(max(abs(squeeze(p_asa(:,y_index,:)))));
mesh(z*1000, x*1000, p_normalized);
xlabel('z (mm)');
ylabel('x (mm)');
shading flat;
title('normalized ASA Pressure (y=0)');
