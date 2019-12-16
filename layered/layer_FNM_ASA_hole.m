%ʹ������FNM��ASA����������֯
clc;
clear all;
close all;
%set up the tranducer

R = 75e-3;
a = 30e-3;
fnumber=R/(2*a);
d = sqrt(R^2 - a^2);
n=0.05:0.01:0.5;
% 
% Use a layered medium
medium = set_layered_medium([0,30e-3],[set_medium('water'),set_medium('muscle')]);

% Center frequency and wavelength
f0=1e6;
lambda = medium(1).soundspeed/f0;
% k=2*pi/lambda;

%���������
xmin=-1.5*a;
xmax=1.5*a;
ymin=-1.5*a;
ymax=1.5*a;
zmin=20e-3;
zmax=90e-3;
dx = lambda/15;
dy = lambda/15;
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
p_asa_all = layerasa(p0_all,z,medium,1024,dz,f0);
tic
for i=1:length(n)
% Create planar array of rectangular transducers
hole_a(i)=n(i)*a;
xdcr_array_hole = get_spherical_shell(hole_a(i),R);
% Calculate the pressure
p0_hole = cw_pressure(xdcr_array_hole,cg_p0,medium,ndiv,f0);%'fnm sse'
p_asa_hole = layerasa(p0_hole,z,medium,1024,dz,f0);

p_asa_res=p_asa_all-p_asa_hole;
p_res_abs=abs(p_asa_res);
p_res_max(i)=max(p_res_abs(:));
focus_index_1(i)=find(p_res_abs==p_res_max(i));

s=size(p_res_abs);
[x_index,y_index,z_index]=ind2sub(s,focus_index_1(i));%�����ֵ���±�תΪ��ά���±�
focus_index(i)=z_index;
focus_forward(i)=R-z(focus_index(i));
%����-6dB
axial_dB_index=find(p_res_abs(x_median,y_median,:)>=0.5*p_res_max(i));
axial_dB(i)=z(max(axial_dB_index))-z(min(axial_dB_index));

radial_dB_index=find(p_res_abs(:,y_median,z_index)>=0.5*p_res_max(i));
radial_dB(i)=x(max(radial_dB_index))-x(min(radial_dB_index));
end
toc
figure(1);
scatter(n,focus_forward,25,'.','k');
hold on;
figure(2);
scatter(n,p_res_max,25,'.','k');
hold on;
figure(3);
scatter(n,axial_dB,25,'.','k');
hold on;
figure(4);
scatter(n,radial_dB,25,'.','k');
hold on;
% 
% 
% %����ȫ����ˮ
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
[x_index,y_index,z_index_w]=ind2sub(s,focus_index_1w(i));%�����ֵ���±�תΪ��ά���±�
focus_index_w(i)=z_index_w;
focus_forward_w(i)=R-z(focus_index_w(i));

%����-6dB
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
