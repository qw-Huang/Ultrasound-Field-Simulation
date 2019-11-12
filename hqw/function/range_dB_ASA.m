function [output]=range_dB_ASA(R,a,f0)
%rayleigh+ASA����3D��ǿ�ֲ����󽹵㾶������-6dB��С��ֵ�ͽ���ǰ�Ƶľ���
%rayleigh+ASA
%clc;clear;close all;
f0=1e6;
u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ�����->ˮ
lambda = medium.soundspeed/f0;%����=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%����

%fnumber=1.15;
%a = 5 * 1.5e-3;%ע������a�ǿ׾���һ��
%R = 2*a*fnumber;%ROC���ʰ뾶
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���

%���������
xmin=-a;%�۲������ķ�Χ
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

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
r_back=0:dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=dr:dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda/6;%����㻷����2pi*a����ɢ�����Ӧ��dS�Ļ���С��lambda/6
ntheta=round(2*pi*(a-dr/2)/Sm);%�����һ�������Ļ��ֵ���,ȡ��
dtheta=2*pi./ntheta;%����ȡ�����µ���ÿ������Դ�Ķ�Ӧ����
theta_after=dtheta:dtheta:2*pi;%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta_back=0:dtheta:(2*pi-dtheta);%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta=theta_after-dtheta/2;%��i��������ɢ�ɵ���Դ��dS�м�ĵ��Ӧ�Ļ�������
X=cos(theta)'*r;%��Դ��x����  ������ʽ ntheta*nr
Y=sin(theta)'*r;%��Դ��y����
Z_ring=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z=repmat(Z_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�

% z0=zmin;  %��z=z0ƽ��������ֲ�
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);%ASA����ά�ռ�����������㻮��

%�Զ���rayleigh���ּ���xy�������  
tic
for ix=1:nx  %�۲������x��������껮��
    for iy=1:ny %�۲������y��������껮��
            rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-zmin).^2);%�۲�㵽��Դ�ľ���
            dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
            dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
            A=dS.*exp(-1i.*k.*rn)./rn;
            B=sum(sum(A));%��������õ�ֵ�ۼ�
            pr(ix,iy)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp        
    end
end

%����3D�����ֲ�
p_asa = cw_angular_spectrum(pr,cg_3d,medium,f0,1024,'Pa');%��������ASA��������ά�����ֲ�
I_asa=acousticintensity(p_asa,medium.density,medium.soundspeed); %��ǿ����
toc

%�ҵ������Ӧ����������
max_index=find_maxpoint(I_asa);%��������λ������
x_index=max_index(1);%���ֵλ��x�±�
y_index=max_index(2);%���ֵλ��y�±�
z_index=max_index(3);%���ֵλ��z�±�

%����ǰ�ƾ���
z_center=fix((nz-1)/2+1);%R��Ӧ��z���±�
range_value_z=(z_center-z_index)*dz;%ǰ�ƾ���

%���ҵ������������Ӧ�����棬�ѷ���-6dB��Χ����ǿ�ŵ�һ���µ�������
I_asa_xy=squeeze(I_asa(:,:,z_index));%��ǿ�����Ӧ��xy��
I_asa_xz=squeeze(I_asa(:,y_index,:));%��ǿ�����Ӧ��xz��
I_range_xy=finddB(I_asa_xy,nx,ny);%xyƽ���ҵ�����-6dB��Χ
I_range_xz=finddB(I_asa_xz,nx,nz);%xzƽ���ҵ�����-6dB��Χ

%xy��������ֵ����Զ��������������Ĳ������Բ�������-6dB��Χ�ľ���
range_point_xy=find_range_point(I_range_xy,nx,ny);  %���ܴ�������Խ�磬����ڱ߽紦���ݲ�Ϊ0�������ͻ��������
range_point_xz=find_range_point(I_range_xz,nx,nz);

range_value_xy=range_point_xy*dy;%��ֵ=����*����
range_value_xz=range_point_xz*dz;

output=[range_value_xy,range_value_xz,range_value_z];
end