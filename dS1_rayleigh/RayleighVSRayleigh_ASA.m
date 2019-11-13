%ʹ������Rayleigh��������ά�ռ��е���ǿ�ֲ�,�Ա�����rayleigh����3D��focus����rayleigh+ASA�ĸ���Ȥ�����  
clc;
clear all;
clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ�����->ˮ���ɸĳɶ�㣬�ú���set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%����

R = 5 * 2 * lambda;%ROC���ʰ뾶
a = 5 * lambda;%ע������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���

%���������
xmin=-a;%�۲������ķ�Χ
xmax=-xmin;
ymax=xmax;
ymin=-ymax;
zDiff=d;
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
error_zeros=zeros(nx,ny,nz);

% rayleigh��FOCUS�����ASA����
tic
xdcr = get_spherical_shell(a,R); %�õ����滻����
delta = [dx dy 0];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmin);%�����Ļ���
ndiv = 100;%���ֵĵ���
dflag = 0;%���=1�������ʲô��ͬ��
pre=rayleigh_cw(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS�Դ���fnm����
z0=zmin;  %��z=z0ƽ��������ֲ�
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);%ASA����ά�ռ�����������㻮��
pre_asa= cw_angular_spectrum(pre,cg_3d,medium,f0,1024,'Pa');%��������ASA��������ά�����ֲ�
% I_pref_asa=acousticintensity(pref_asa,medium.density,medium.soundspeed); %��ǿ����
toc

%rayleigh��FOCUS��
tic
pre=rayleigh_cw(xdcr,cg_3d,medium,ndiv,f0,dflag);%FOCUS�Դ���rayleigh����
toc

%�ҵ�����Ȥ��������ǿ-6dB��Χ�����������Χ�ڵ�error;�ҵ����ֵ�����ض�Ӧ����ά���꣬�������Ӧ��xyƽ�棬�Ѹ���Ȥ���򻭳���
max_index=find_maxpoint(abs(pre));%��������λ������
x_index=max_index(1);%���ֵλ��x�±�
y_index=max_index(2);%���ֵλ��y�±�
z_index=max_index(3);%���ֵλ��z�±�
pre_asa_max=max(abs(pre_asa(:)));%�Զ���rayleigh 3D��ǿ���ֵ
pre_max=max(abs(pre(:)));

error=abs((pre_asa-pre)./pre_asa);%������ά�ռ������ַ��������

%ѡ����ά�ռ����Ȥ�������
[index]=find(abs(pre)>=0.25*pre_max);
error_zeros(index)=error(index);

%��������Ȥ�������
figure(1);
axis equal;
surf(error_zeros(:,:,z_index));
shading interp
colorbar;
xlabel('x');
ylabel('y');
title('Rayleigh(FPCUS)3D VS Rayleigh(FOCUS)+ASA ���㴦xyƽ��������');

figure(2);
axis equal;
error_xz=squeeze(error_zeros(:,y_index,:));
surf(error_xz);
shading interp
colorbar;
xlabel('z');
ylabel('x');
title('Rayleigh(FPCUS)3D VS Rayleigh(FOCUS)+ASA ���㴦xzƽ��������');

figure(3);
plot(max(error_xz));
xlabel('z');
ylabel('error');
title('rayleigh 3D VS Rayleigh(FOCUS)+ASA ����Ȥ��xy����������ı仯���');
% figure(4);
% plot(zz,error_average);
% xlabel('z');
% ylabel('error');
% title('rayleigh 3D VS Rayleigh(FOCUS)+ASA ����Ȥ��xy����ƽ�����ı仯���');


