%��������xzƽ��ƫתһ���Ƕ�
clear all;
f0=1e6;%����Ƶ�ʺͷ�������
P=100;
angle_rotation=pi/36;
medium = set_medium('muscle');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
dBperNeper = 20 * log10(exp(1));
attenuationNeperspermeter=medium.attenuationdBcmMHz/dBperNeper*100*f0/1e6;
k=2*pi/lambda;%����

R = 75e-3;%ROC���ʰ뾶
a=30e-3;
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

%���������
xmin=-1.5*a-10e-3;%�۲������ķ�Χ
xmax=1.5*a-10e-3;
ymax=0;
ymin=0;

zmin=10e-3;
% zmin=R-zDiff;
zmax=90e-3;

dx = lambda/6; %�����Ĳ���
dz = lambda/6;

x=xmin:dx:xmax;%�����ķֲ�
z=zmin:dz:zmax;

nx=length(x);%�����ĵ���
nz=length(z);

pr_zeros_xz=zeros(nx,nz);
Xp=zeros(nx,nz);
Zp=zeros(nx,nz);

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
r_back=0:dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=dr:dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda/6;%�м价����ɢ�����Ӧ��dS�Ļ�������lambda/6
median=round(length(r)/2);%ȡ�м价��������
ntheta=round(2*pi*r(median)/Sm);%�м��һ�������Ļ��ֵ���,ȡ��
dtheta=2*pi./ntheta;%����ȡ�����µ���ÿ������Դ�Ķ�Ӧ����
theta_after=dtheta:dtheta:2*pi;%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta_back=0:dtheta:(2*pi-dtheta);%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta=theta_after-dtheta/2;%��i��������ɢ�ɵ���Դ��dS�м�ĵ��Ӧ�Ļ�������
% X=cos(theta)'*r;%��Դ��x����  ������ʽ ntheta*nr
X=cos(theta)'*r;
Y=sin(theta)'*r;%��Դ��y����
% Z0=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z0=(R-sqrt(R*R-r.*r));
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
Rx=0;Rz=R; %�Ƽ��ν�����ת����xzƽ����
Z_angle=cos(angle_rotation)*(Z-Rz)-sin(angle_rotation)*(X-Rx)+Rz; %��ά�ռ���������ת��ʽ
X_angle=sin(angle_rotation)*(Z-Rz)+cos(angle_rotation)*(X-Rx)+Rx;

% ��y=0ʱxzƽ��
y=0;

%rayleigh���ּ���xzƽ�������  
tic
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������
        rn=sqrt((X_angle-x(ix)).^2+(Y-y).^2+(Z_angle-z(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*k.*rn)./rn.*exp(-attenuationNeperspermeter.*rn);
        B=sum(sum(A));%��������õ�ֵ�ۼ� 
        pr(ix,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc
pr_max=max(abs(pr(:)));
pr_abs=abs(pr);
%ѡ�����Ȥ����
pr_nor=abs(pr)./pr_max;
[index]=find(pr_nor>=0.5);
pr_zeros_xz(index)=pr_nor(index);
X_all=repmat(x',1,nz);
Xp(index)=X_all(index);
Z_all=repmat(z,nx,1);
Zp(index)=Z_all(index);
%������Ȥ������ת�غ�z��ƽ�е�λ�ã��ٽ��бȽ�
Zp_angle=cos(-angle_rotation)*(Zp-Rz)-sin(-angle_rotation)*(Xp-Rx)+Rz; 
Xp_angle=sin(-angle_rotation)*(Zp-Rz)+cos(-angle_rotation)*(Xp-Rx)+Rx;
Xp(index)=Xp_angle(index);%����Ϊ������ת�س����z���غ�
Zp(index)=Zp_angle(index);
axial=max(Zp(:))-min(Zp(:));%������ ��С��0
radial=max(Xp(:))-min(Xp(:));

figure(3)
surf(z*1000,x*1000,abs(pr)./pr_max);
shading interp;
colorbar;
axis equal;
