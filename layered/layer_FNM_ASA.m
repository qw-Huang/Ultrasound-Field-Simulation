%使用内置FNM和ASA，计算多层组织声强分布，代码有问题，需要重新看
clc;
clear all;
%set up the tranducer
ntheta=70;
R = 5 * 3 * lambda;
a = 5 * 1.5 * lambda;
fnumber=R/(2*a);
d = sqrt(R^2 - a^2);

% Use a layered medium
medium = set_layered_medium([0,2e-2,22e-3],[set_medium('water'),set_medium('skin'),set_medium('fat')]);%水->皮肤->脂肪

% Center frequency and wavelength
f0=1e6;
lambda = medium(1).soundspeed/f0;
k=2*pi/lambda;

%划分网格点
xmin=-a;
xmax=a;
ymin=-a;
ymax=a;
zmin=0;
zmax=30*lambda;
dx=lambda/6;
dy=lambda/6;
dz=lambda/6;
x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;
nx=(xmax-xmin)/dx+1;
ny=(ymax-ymin)/dy+1;
nz=(zmax-zmin)/dz+1;

% Create planar array of rectangular transducers
xdcr_array = get_spherical_shell(a,R);% 球形换能器

% Determine where the source pressure will be calculated
z0 = 2e-2;
y_index = floor((ymax-ymin)/2/dy);
% Coordinate grids to calclate the initial pressure (x-y plane) and final
% pressure (x-z plane)
cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);

% Calculate the pressure
ndiv = 100;
fprintf('Calculating initial pressure plane with FNM... ');
tic();
p0 = cw_pressure(xdcr_array,cg_p0,medium,ndiv,f0);%'fnm sse'
fprintf('done in %f s.\n', toc());

fprintf('Calculating 3D pressure (%i points) with ASA... ', (length(x) * length(y) * length(z)));
tic();
p_asa = layerasa(p0,z,medium,1024,dz,f0);
fprintf('done in %f s.\n', toc());

% Show the initial pressure
figure(1);
pcolor(x*1000, y*1000, rot90(abs(squeeze(p0(:,:,1)))));
xlabel('x (mm)');
ylabel('y (mm)');
shading flat;
title(['p0 (Calculated with FNM at z = ', num2str(z0*1000), ' mm)']);
% Show the 3D field calculated with ASA
figure(2);
p_normalized=abs(squeeze(p_asa(:,y_index,:)))/max(max(abs(squeeze(p_asa(:,y_index,:)))));
pcolor(z*1000, x*1000, p_normalized);
xlabel('z (mm)');
ylabel('x (mm)');
shading flat;
title('normalized ASA Pressure (y=0)');
