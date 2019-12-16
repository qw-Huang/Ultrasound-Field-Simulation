% 使用自定义rayleigh+ASA计算3D声强分布，针对多层组织（或者单层有衰减的组织）
clc;
clear all;

% Set up the array

R=15e-3;
a=7.5e-3;
u=1;
% Use a layered medium
medium = set_layered_medium([0,5e-3],[set_medium('water'),set_medium('muscle')]);

% Center frequency and wavelength
f0 = 1e6;
lambda = medium(1).soundspeed/f0;
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离


% Set up the coordinate grid
xmin = -7.5e-3;
xmax = -xmin;
ymin = -7.5e-3;
ymax = -ymin;
zmin = R-0.9*d;
zmax = R+10e-3;



dx = lambda/6;
dy = lambda/6;
dz = lambda/6;

x = xmin:dx:xmax;
y = ymin:dy:ymax;
z = zmin:dz:zmax;

% Determine where the source pressure will be calculated
z0 = R-0.9*d;
y_index = floor((ymax-ymin)/2/dy);



pr=rayleigh_2D_xy(R,a,f0,u,medium(1),x,y,zmin);

% % Focus the array
% xdcr_array = find_single_focus_phase(xdcr_array,focus_x,focus_y,focus_z,medium,f0,200);

% % Calculate the pressure
% ndiv = 10;
% fprintf('Calculating p0 with FNM... ');
% tic();
% p0 = cw_pressure(xdcr_array,cg_p0,medium,ndiv,f0);
% fprintf('done in %f s.\n', toc());

tic();
p_asa = layerasa(pr,z,medium,1024,dz,f0);
fprintf('done in %f s.\n', toc());

% figure(1);
% pcolor(x*1000, y*1000, rot90(abs(squeeze(p0(:,:,1)))));
% xlabel('x (mm)');
% ylabel('y (mm)');
% shading flat;
% title(['p0 (Calculated with FNM at z = ', num2str(z0*1000), ' mm)']);
% 
% figure(2);
% pcolor(z*1000, x*1000, abs(squeeze(p_asa(:,y_index,:))));
% xlabel('z (mm)');
% ylabel('x (mm)');
% shading flat;
% title('ASA Pressure (y=0)');
