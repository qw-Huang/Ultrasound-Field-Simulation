%ʹ��focus����FNM����ά�ռ��е���ǿ�ֲ�,�ԱȺ�focus����fnm+ASA�ĸ���Ȥ�����  
clc;
clear all;
clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ�����->ˮ���ɸĳɶ�㣬�ú���set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����

R = 75e-3;%ROC���ʰ뾶
a = 30e-3;%ע������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
%���������
xmin=-1.5*a;%�۲������ķ�Χ
xmax=-xmin;
ymax=xmax;
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
error_zeros=zeros(nx,ny,nz);

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
ntheta=100;%
dtheta=2*pi/ntheta;dr=lambda/6;
theta=dtheta:dtheta:2*pi;
r0=0:dr:a-dr;%
r=dr:dr:a;
X=cos(theta)'*r;%��Դ��x����  ������ʽ ntheta*nr
Y=sin(theta)'*r;%��Դ��y����
Z0=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z=repmat(Z0,ntheta,1); %����Դ��ά�ռ��е�z���� repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�

z0=zmin;  %��z=z0ƽ��������ֲ�
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);%ASA����ά�ռ�����������㻮��

%FNM����ά�ռ�
xdcr = get_spherical_shell(a,R); %�õ����滻����
delta = [dx dy dz];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%�����Ļ���
ndiv = 80;%���ֵĵ���
dflag = 0;%���=1�������ʲô��ͬ��
tic
pref=fnm_call(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS�Դ���fnm����
I_pref= acousticintensity(pref,medium.density,medium.soundspeed); %��ǿ���� 
toc

% FNM���ASA����
tic
delta2 = [dx dy 1];%����㲽��
ps2 = set_coordinate_grid(delta2, xmin, xmax, ymin, ymax, zmin, zmin);%�����Ļ���
pref=fnm_call(xdcr,ps2,medium,ndiv,f0,dflag);%FOCUS�Դ���fnm����
pref_asa= cw_angular_spectrum(pref,cg_3d,medium,f0,1024,'Pa');%��������ASA��������ά�����ֲ�
I_pref_asa=acousticintensity(pref_asa,medium.density,medium.soundspeed); %��ǿ����
toc


%�ҵ�����Ȥ��������ǿ-6dB��Χ�����������Χ�ڵ�error;�ҵ����ֵ�����ض�Ӧ����ά���꣬�������Ӧ��xyƽ�棬�Ѹ���Ȥ���򻭳���
max_index=find_maxpoint(I_pref);%��������λ������
x_index=max_index(1);%���ֵλ��x�±�
y_index=max_index(2);%���ֵλ��y�±�
z_index=max_index(3);%���ֵλ��z�±�
I_pref_asa_max=max(I_pref_asa(:));%�Զ���rayleigh 3D��ǿ���ֵ
I_pref_max=max(I_pref(:));


error_I=abs(I_pref_asa-I_pref)./I_pref;%������ά�ռ������ַ��������
%ѡ����ά�ռ����Ȥ�������
[index]=find(I_pref>=0.25*I_pref_max);
error_zeros(index)=error_I(index);

%��������Ȥ�������
figure(1);
surf(error_zeros(:,:,z_index));
shading interp;
axis equal;
xlabel('x');
ylabel('y');
title('FNM 3D VS FNM+ASA ���㴦xyƽ��������');

figure(2);
error_xz=squeeze(error_zeros(:,y_index,:));
surf(error_xz);
shading interp;
axis equal;
xlabel('z');
ylabel('x');
title('FNM 3D VS FNM+ASA ���㴦xzƽ��������');


figure(3);
I_pref_xz=finddB(I_pref(:,y_index,:),nx,nz);
pcolor(I_pref_xz);
axis equal
shading interp; 
colorbar;
xlabel('z');
ylabel('x');
title('FNM ���㴦xzƽ����ǿ�ֲ����');

figure(4);
I_pref_asa_xz=finddB(I_pref_asa(:,y_index,:),nx,nz);
pcolor(I_pref_asa_xz);
axis equal
shading interp; 
colorbar;
xlabel('z');
ylabel('x');
title('FNM+ASA ���㴦xzƽ����ǿ�ֲ����');
figure(6);
plot(max(error_xz));
xlabel('z');
ylabel('error');
title('FNM VS FNM+ASA ����Ȥ��xz����������ı仯���');

% zz=17:1:61;
% for i=1:length(zz)
%     error_max(i)=max(max(error_zeros(:,:,zz(i))));
%     error_average(i)=mean(mean(error_zeros(:,:,zz(i))));
%     %error_median(i)=median(error_zeros(:,:,zz(i)),'all');
% end
% 
% figure(3);
% plot(zz,error_max);
% xlabel('z');
% ylabel('error');
% title('FNM 3D VS FNM+ASA ����Ȥ��xy����������ı仯���');
% figure(4);
% plot(zz,error_average);
% xlabel('z');
% ylabel('error');
% title('����Ȥ��xy����ƽ�����ı仯���');


