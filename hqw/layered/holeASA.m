%function [output]=rayleigh_ASA(R,a,u,f0)
%����R,a,u,f0,����м俪�׵�3D��ǿ�ֲ���ʹ���Զ���rayleigh����+ASA������ǿ
clc;clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������

a = 5 * 1.5e-3;%ע������a�ǿ׾���һ��
R = 2*a;%ROC���ʰ뾶
rhole=0.004;
medium = set_medium('lossless');%������ʣ�����->ˮ
lambda = medium.soundspeed/f0;%����=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%����
%fnumber=1;
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
%���������
xmin=-0.7*a;%�۲������ķ�Χ
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
ntheta=70;%
dtheta=2*pi/ntheta;dr=lambda/6;
theta=dtheta:dtheta:2*pi;
r0=rhole:dr:a-dr;%
r=(rhole+dr):dr:a;
X=cos(theta)'*r;%��Դ��x����  ������ʽ ntheta*nr
Y=sin(theta)'*r;%��Դ��y����
Z0=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z=repmat(Z0,ntheta,1); %����Դ��ά�ռ��е�z���� repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�

z0=zmin;  %��z=z0ƽ��������ֲ�
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);%ASA����ά�ռ�����������㻮��

%�Զ���rayleigh���ּ���xy�������  
tic
for ix=1:nx  %�۲������x��������껮��
    for iy=1:ny %�۲������y��������껮��
           rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z0).^2);%�۲�㵽��Դ�ľ���
           dS0=r.*dtheta.*R.*(asin(r./R)-asin(r0./R));%��Դ���������Ϊ����,���=����Բ��������Բ����
           dS=repmat(dS0,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
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
I_xz=squeeze(I_asa(:,y_index,:));
mesh(I_xz);
%end