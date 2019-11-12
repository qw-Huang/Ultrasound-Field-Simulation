%ʹ������FNM��FOCUS������3D�ռ������ֲ�
clc;
clear all;
f0=1e6;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f

k=2*pi/lambda;%����

R = 5 * 2 * lambda;%ROC���ʰ뾶
a = 5 * lambda;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���

%���������
xmin=-0.7*a;%�۲������ķ�Χ
xmax=-xmin;
ymax=0.7*a;
ymin=-ymax;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %�����Ĳ���
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);
ny=length(y);
nz=length(z);

%����focus����FNM����fnm_call()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr = get_spherical_shell(a,R); %�õ����滻����

if nz > 1,
    dz = 2 * d / (nz - 1);
else
    dz = 0;
end
if nx > 1,
    dx = 2 * a / (nx - 1);
else
    dx = 0;
end
if ny > 1,
    dy = 2 * a / (ny - 1);
else
    dy = 0;
end

delta = [dx dy dz];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%�����Ļ���
ndiv = 60;%���ֵĵ���
dflag = 0;%���=1�������ʲô��ͬ��

tic
pref=fnm_call(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS�Դ���rayleigh����
toc

I=pref.*pref/(medium.soundspeed*medium.density);

z_index=51;
z_index_value=(z_index-1)*dz;

% plot the pressures
if nx > 1 & nz > 1 & ny > 1,
    figure(1);
    axis equal;
    If_normalized=abs(squeeze(I(:,:,z_index)))/abs(max(max(squeeze(I(:,:,z_index)))));
    mesh(If_normalized);%��prs��һ��   Ϊʲô��Ӧ��������z��ǰ�棿
    If_xy=squeeze(I(:,:,z_index))/abs(max(max(max(I))));%������rayleigh��prs��һ������������ά�ռ������ĵ㣬�õ�����
    colorbar
    title('Rayleigh Sommerfeld Result');
    xlabel('nx ');
    ylabel('ny ');
    zlabel('normalized intensity');
end
