%rayleigh积分+ASA 改变肌肉的层厚，即改变z0的位置
clc;
close all;
clear all;

% Set up the array
R=75e-3;
a=30e-3;

% Use a layered medium
medium1 = set_medium('water');
medium2 = set_medium('muscle');
P=100;

d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,0,medium1.density,medium1.soundspeed);

z_interface=30e-3;

% Set up the coordinate grid
xmin = -1.5*a;
xmax = -xmin;
ymin = -1.5*a;
ymax = -ymin;
zmin = z_interface;
zmax = 90e-3;

dx = 2.5e-4;
dy = 2.5e-4;
dz = 2.5e-4;

x = xmin:dx:xmax;
y = ymin:dy:ymax;
z = zmin:dz:zmax;

z0=7.75e-3;
zmin1=z0;
zmax1=z_interface;
z1=zmin1:dz:zmax1;

% Determine where the source pressure will be calculated
y_median = floor((ymax-ymin)/2/dy)+1;
x_median=floor(length(x)/2)+1;

xdcr_array = get_spherical_shell(a,R);% 球形换能器

cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
cg_3d1 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin1,zmax1);
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);
   
f=1e6:0.1e6:2.5e6;
% f=[1e6 1.3e6 1.6e6 1.9e6 2.2e6 2.5e6];
for i=1:length(f)
    f0 = f(i);
    lambda = medium1.soundspeed/f0;
    tic
    ndiv = 200;
    p0 = fnm_call(xdcr_array,cg_p0,medium1,ndiv,f0,0);
    p0_u=p0.*u;
    [p_asa1,p_interface]=layer_cw_angular_spectrum(p0_u,cg_3d1,medium1,medium2,f0,1024,'Pa');
    p_asa=cw_angular_spectrum(p_interface,cg_3d,medium2,f0,1024,'Pa');%衰减系数要考虑频率的影响
   toc
    p_asa_abs=abs(p_asa);
    p_asa_max(i)=max(p_asa_abs(:));
    focus_index_1(i)=find(p_asa_abs==p_asa_max(i));
    s=size(p_asa_abs);
    [x_index,y_index,z_index]=ind2sub(s,focus_index_1(i));%将最大值单下标转为三维多下标
    focus_index(i)=z_index;
    focus_forward(i)=R-z(focus_index(i));
%     %轴向-6dB
%     axial_dB_index=find(p_asa_abs(x_median,y_median,:)>=0.5*p_asa_max(i));
%     axial_dB(i)=z(max(axial_dB_index))-z(min(axial_dB_index));
%     radial_dB_index=find(p_asa_abs(:,y_median,z_index)>=0.5*p_asa_max(i));
%     radial_dB(i)=x(max(radial_dB_index))-x(min(radial_dB_index));

end
toc
figure(1);
scatter(f,focus_forward,25,'.','k');
xlabel('f0(Hz)');
ylabel('focus forward distance(m)');
figure(2);
scatter(f,p_asa_max,25,'.','k');
xlabel('f0(Hz)');
ylabel('pmax(Pa)');
% figure(3);
% scatter(z_interface,axial_dB,25,'.','k');
% xlabel('z-interface(m)');
% ylabel('axial -6dB distance(m)');
% figure(4);
% scatter(z_interface,radial_dB,25,'.','k');
% xlabel('z-interface(m)');
% ylabel('radial -6dB distance(m)');
