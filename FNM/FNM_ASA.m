%内置FNM(FOCUS)计算一个平面，然后用ASA计算三维空间声场分布
clc;clear all;
f0=1e6;u=1;%定义频率和法向阵速
a = 30e-3;%注意这里a是孔径的一半
R = 75e-3;%ROC曲率半径
medium = set_medium('water');%定义介质：单层->水
lambda = medium.soundspeed/f0;%波长=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%波数
fnumber=1;
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
%划分网格点
xmin=-1.5*a;%观察点坐标的范围
xmax=-xmin;
ymax=1.5*a;
ymin=-ymax;

zmin=30e-3;
zmax=30e-3;

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
% p_asa = cw_angular_spectrum(p0,cg_3d,medium,f0,1024,'Pa');
fprintf('done in %f s.\n', toc());

% Show the initial pressure
pr=rayleigh_2D_xy(R,a,f0,u,medium,x,y,zmax);
error=abs(pr-p_asa(:,:,length(z)))./abs(pr);
