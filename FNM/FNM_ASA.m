%����FNM(FOCUS)����һ��ƽ�棬Ȼ����ASA������ά�ռ������ֲ�
clc;clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������
a = 30e-3;%ע������a�ǿ׾���һ��
R = 75e-3;%ROC���ʰ뾶
medium = set_medium('water');%������ʣ�����->ˮ
lambda = medium.soundspeed/f0;%����=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%����
fnumber=1;
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
%���������
xmin=-1.5*a;%�۲������ķ�Χ
xmax=-xmin;
ymax=1.5*a;
ymin=-ymax;

zmin=30e-3;
zmax=30e-3;

dx = lambda/6; %�����Ĳ���
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);
ny=length(y);
nz=length(z);

% Create planar array of rectangular transducers
xdcr_array = get_spherical_shell(a,R);% ���λ�����

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
