%使用内置FNM和ASA，计算多层组织  加开孔 但是结果不能用，内置FNM只能保证u相等，无法保证每一次的声功率相等
clc;
clear all;
close all;
%set up the tranducer

R = 75e-3;
a = 30e-3;
fnumber=R/(2*a);
d = sqrt(R^2 - a^2);
P=100;
n=0.05:0.05:0.5;
% 
% Use a layered medium
medium = set_layered_medium([0,30e-3],[set_medium('water'),set_medium('muscle')]);

% Center frequency and wavelength
f0=1e6;
lambda = medium(1).soundspeed/f0;

%划分网格点
xmin=-1.5*a;
xmax=1.5*a;
ymin=-1.5*a;
ymax=1.5*a;
zmin=29.75e-3;
zmax=90e-3;
dx = lambda/6;
dy = lambda/6;
dz = lambda/6;
x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;
nx=length(x);
ny=length(y);
nz=length(z);
% Determine where the source pressure will be calculated
z0 = 30e-3;
x_median = floor((xmax-xmin)/2/dx)+1;
y_median = floor((ymax-ymin)/2/dy)+1;

% pressure (x-z plane)
cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);

xdcr_array_all = get_spherical_shell(a,R);
% Calculate the pressure
ndiv = 200;
p0_all = cw_pressure(xdcr_array_all,cg_p0,medium,ndiv,f0);%'fnm sse'

tic
for i=1:length(n)
% Create planar array of rectangular transducers
hole_a(i)=n(i)*a;
    u(i)=normal_velocity(P,R,a,hole_a(i),medium1.density,medium1.soundspeed);
xdcr_array_hole = get_spherical_shell(hole_a(i),R);
% Calculate the pressure
p0_hole = cw_pressure(xdcr_array_hole,cg_p0,medium,ndiv,f0);%'fnm sse'
p0=p0_all-p0_hole;
p0_u=p0.*u(i);
p_asa = layerasa(p0_u,z,medium,1024,dz,f0);
p_asa_abs=abs(p_asa);
p_asa_max(i)=max(p_asa_abs(:));
focus_index_1(i)=find(p_asa_abs==p_asa_max(i));

s=size(p_asa_abs);
[x_index,y_index,z_index]=ind2sub(s,focus_index_1(i));%将最大值单下标转为三维多下标
focus_index(i)=z_index;
focus_forward(i)=R-z(focus_index(i));
%轴向-6dB
axial_dB_index=find(p_asa_abs(x_median,y_median,:)>=0.5*p_asa_max(i));
axial_dB(i)=z(max(axial_dB_index))-z(min(axial_dB_index));

radial_dB_index=find(p_asa_abs(:,y_median,z_index)>=0.5*p_asa_max(i));
radial_dB(i)=x(max(radial_dB_index))-x(min(radial_dB_index));
end
toc
figure(1);
scatter(n,focus_forward,25,'.','k');
hold on;
figure(2);
scatter(n,p_asa_max,25,'.','k');
hold on;
figure(3);
scatter(n,axial_dB,25,'.','k');
hold on;
figure(4);
scatter(n,radial_dB,25,'.','k');
hold on;
% 
% 
% %介质全部是水
medium_w = set_layered_medium([0,30e-3],[set_medium('water'),set_medium('water')]);
ndiv = 200;
p0_all_w = cw_pressure(xdcr_array_all,cg_p0,medium_w,ndiv,f0);%'fnm sse'
p_asa_all_w = layerasa(p0_all_w,z,medium_w,1024,dz,f0);
tic
for i=1:length(n)
% Create planar array of rectangular transducers
hole_a(i)=n(i)*a;
xdcr_array_hole = get_spherical_shell(hole_a(i),R);
% Calculate the pressure
p0_hole_w = cw_pressure(xdcr_array_hole,cg_p0,medium_w,ndiv,f0);%'fnm sse'
p_asa_hole_w = layerasa(p0_hole_w,z,medium_w,1024,dz,f0);

p_asa_res_w=p_asa_all_w-p_asa_hole_w;
p_res_abs_w=abs(p_asa_res_w);
p_res_max_w(i)=max(p_res_abs_w(:));
focus_index_1w(i)=find(p_res_abs_w==p_res_max_w(i));

s=size(p_res_abs_w);
[x_index,y_index,z_index_w]=ind2sub(s,focus_index_1w(i));%将最大值单下标转为三维多下标
focus_index_w(i)=z_index_w;
focus_forward_w(i)=R-z(focus_index_w(i));

%轴向-6dB
axial_dB_index_w=find(p_res_abs_w(x_median,y_median,:)>=0.5*p_res_max_w(i));
axial_dBw(i)=z(max(axial_dB_index_w))-z(min(axial_dB_index_w));

radial_dB_indexw=find(p_res_abs_w(:,y_median,z_index_w)>=0.5*p_res_max_w(i));
radial_dBw(i)=x(max(radial_dB_indexw))-x(min(radial_dB_indexw));

end
toc
figure(1);
scatter(n,focus_forward_w,25,'.','b');
xlabel('hole/a');
ylabel('focus forward distance(m)');
figure(2);
scatter(n,p_res_max_w,25,'.','b');
xlabel('hole/a');
ylabel('pmax(Pa)');
figure(3);
scatter(n,axial_dBw,25,'.','b');
xlabel('hole/a');
ylabel('axial -6dB distance(m)');
figure(4);
scatter(n,radial_dBw,25,'.','b');
xlabel('hole/a');
ylabel('radial -6dB distance(m)');
