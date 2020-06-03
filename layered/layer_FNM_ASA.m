%ʹ������FNM��ASA������������֯��ǿ�ֲ������ǽ���ķ��������

clear all;
clear all;
%set up the tranducer

R = 150e-3;
a = 90e-3;
fnumber=R/(2*a);
d = sqrt(R^2 - a^2);
P=100;

% Use a layered medium
medium1 = set_medium('muscle');%
medium2 = set_medium('muscle');

% Center frequency and wavelength
f0=1.2e6;
lambda = medium1.soundspeed/f0;
u=normal_velocity(P,R,a,0,medium1.density,medium1.soundspeed);
depth=0.06;
z_interface=R-depth;

%���������
dx=lambda/6;
dy=lambda/6;
dz=lambda/6;
xmin=-1.5*a;
xmax=1.5*a;
ymin=-1.5*a;
ymax=1.5*a;
zmin2=z_interface;
zmax2=151e-3;
x=xmin:dx:xmax;
y=ymin:dy:ymax;
z2=zmin2:dz:zmax2;
nx=length(x);
ny=length(y);
nz=length(z2);

% Create planar array of rectangular transducers
xdcr_array = get_spherical_shell(a,R);% ���λ�����

% Determine where the source pressure will be calculated
% z0 = R-d+lambda/4;
z0=R-d+lambda;
y_median = floor((ymax-ymin)/2/dy)+1; %floor ��ʹ���е�����
x_median=floor(length(x)/2)+1;

zmin1=z0;
zmax1=z_interface;
z1=zmin1:dz:zmax1;

% pressure (x-z plane)
cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
cg_3d1 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin1,zmax1);
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin2,zmax2);

% Calculate the pressure
ndiv = 200;
tic();
p0 = fnm_call(xdcr_array,cg_p0,medium1,ndiv,f0,0);%'fnm sse'
p0_u=p0.*u;
%  p0_u=p0;
[p_asa1,p_interface]=layer_cw_angular_spectrum(p0_u,cg_3d1,medium1,medium2,f0,2048,'Pa');
p_asa=cw_angular_spectrum(p_interface,cg_3d,medium2,f0,2048,'Pa');
toc
p_asa_abs=abs(p_asa);
p_asa_max=max(p_asa_abs(:));
focus_index_1=find(p_asa_abs==p_asa_max);
s=size(p_asa_abs);
[x_index,y_index,z_index]=ind2sub(s,focus_index_1);%�����ֵ���±�תΪ��ά���±�
%����-6dB
axial_dB_index=find(p_asa_abs(x_median,y_median,:)>=0.5*p_asa_max);
axial_dB=z2(max(axial_dB_index))-z2(min(axial_dB_index));

radial_dB_index=find(p_asa_abs(:,y_median,z_index)>=0.5*p_asa_max);
radial_dB=x(max(radial_dB_index))-x(min(radial_dB_index));
%������ƴ������
p1=squeeze(p_asa1(:,y_index,:));
p2=p_asa(:,y_index,:);
p2(:,1)=[];
p=[p1 p2];
p_abs=abs(p);
z2(1)=[];
z=[z1 z2];


figure(3);
surf(z*1000, x*1000, p_abs);
xlabel('z (mm)');
ylabel('x (mm)');
axis equal;
shading flat;
title('FNM + ASA ');
