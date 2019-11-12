%ʹ������rayleigh����3D�����ֲ���(����xyƽ�������ֲ������ǽ�����ڴ���) 
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
ymin=xmin;
ymax=xmax;
zdiff=0.7*d;
zmin=R-zdiff;
zmax=R+zdiff;

nx = 61;%�����ķָ���� 
ny = 61; 
nz = 101;

dx = (xmax-xmin)/(nx-1); %�����Ĳ���
dy = (ymax-ymin)/(ny-1); 
dz = (zmax-zmin)/(nz-1);

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

%����focus����rayleigh����rayleigh_cw()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr = get_spherical_shell(a,R); %�õ����滻����

delta = [dx dy dz];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%�����Ļ���
ndiv =100;%���ֵĵ���
dflag = 0;%���=1�������ʲô��ͬ��

tic
prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS�Դ���rayleigh����
I_prs= acousticintensity(prs,medium.density,medium.soundspeed); %��ǿ���� 
toc


%�ҵ�����Ȥ��������ǿ-6dB��Χ�����������Χ�ڵ�error;�ҵ����ֵ�����ض�Ӧ����ά���꣬�������Ӧ��xyƽ�棬�Ѹ���Ȥ���򻭳���
max_index=find_maxpoint(I_prs);%��������λ������
x_index=max_index(1);%���ֵλ��x�±�
y_index=max_index(2);%���ֵλ��y�±�
z_index=max_index(3);%���ֵλ��z�±�
I_prs_max=max(I_prs(:));%�Զ���rayleigh��ǿ���ֵ


% plot the pressures
figure(1);
surf(abs(squeeze(I_prs(:,:,z_index)))/abs(I_prs_max));%��prs��һ��   Ϊʲô��Ӧ��������z��ǰ��
axis equal;
shading interp;
colorbar;
title('Rayleigh Sommerfeld Result(z=R��');
xlabel('x ');
ylabel('y ');
zlabel('normalized pressure');
figure(2);
    axis equal;
    surf(abs(squeeze(I_prs(:,y_index,:)))/abs(I_prs_max));%��prs��һ��   Ϊʲô��Ӧ��������z��ǰ��
    colorbar
    title('Rayleigh Sommerfeld Result(y=0)');
    xlabel('z');
    ylabel('x');
    zlabel('normalized pressure');

