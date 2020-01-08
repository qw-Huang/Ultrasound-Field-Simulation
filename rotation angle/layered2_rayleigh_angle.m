%�Ľ���rayleigh���֣���������������еķֲ� ˮ30mm+����30mm-90mm
% ��ƫת�ǣ��������������ѹ�Ͳ���ƫת�ǵĲ�һ��

% clc;
clear all;
j=1i;
f0=1e6;
P=100;%����Ƶ�ʺͷ�������
angle_rotation= pi/18;
medium = set_layered_medium([0,30e-3],[set_medium('water'),set_medium('muscle')]);
lambda1 = medium(1).soundspeed/f0;%����=c/f
lambda2 = medium(2).soundspeed/f0;%����=c/f

dBperNeper = 20 * log10(exp(1));
attenuationNeperspermeter1=medium(1).attenuationdBcmMHz/dBperNeper*100*f0/1e6;
k1=2*pi/lambda1-j*attenuationNeperspermeter1;%����
attenuationNeperspermeter2=medium(2).attenuationdBcmMHz/dBperNeper*100*f0/1e6;
k2=2*pi/lambda2-j*attenuationNeperspermeter2;%����

R = 75e-3;%ROC���ʰ뾶
a = 30e-3;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,0,medium(1).density,medium(1).soundspeed);
%���������
% dx = lambda1/6; %�����Ĳ���
dx=2.5e-4;
dy = 2.5e-4;
dz = 2.5e-4;

xmin=-1.5*a-10e-3;%�۲������ķ�Χ
xmax=1.5*a-10e-3;
ymax=1.5*a;
ymin=-ymax;

zmin1=7.75e-3;
zmax1=30e-3;
zmin2=30e-3+dz;
zmax2=90e-3;

x=xmin:dx:xmax;%�����ķֲ�
y=ymin:dy:ymax;
z2=zmin2:dz:zmax2;
z1=zmin1:dz:zmax1;

nx=length(x);%�����ĵ���
ny=length(y);
nz1=length(z1);
nz2=length(z2);

pr_zeros_xz=zeros(nx,nz2);
Xp=zeros(nx,nz2);
Zp=zeros(nx,nz2);

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda1/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
r_back=0:dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=dr:dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda1/6;%�м价����ɢ�����Ӧ��dS�Ļ�������lambda/6
median=length(r)/2+1;%ȡ�м价��������
ntheta=round(2*pi*r(median)/Sm);%�м��һ�������Ļ��ֵ���,ȡ��
dtheta=2*pi./ntheta;%����ȡ�����µ���ÿ������Դ�Ķ�Ӧ����
theta_after=dtheta:dtheta:2*pi;%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta_back=0:dtheta:(2*pi-dtheta);%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta=theta_after-dtheta/2;%��i��������ɢ�ɵ���Դ��dS�м�ĵ��Ӧ�Ļ�������
X=cos(theta)'*r;%��Դ��x����  ������ʽ ntheta*nr
Y=sin(theta)'*r;%��Դ��y����
Z0=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
Rx=0;Rz=R; %�Ƽ��ν�����ת����xzƽ����
Z_angle=cos(angle_rotation)*(Z-Rz)-sin(angle_rotation)*(X-Rx)+Rz; %��ά�ռ���������ת��ʽ
X_angle=sin(angle_rotation)*(Z-Rz)+cos(angle_rotation)*(X-Rx)+Rx;

% rayleigh���ּ���ˮ��xzƽ�������  
tic
for ix=1:nx  %�۲������x���������
    for iz=1:nz1  %�۲������z���������
        rn=sqrt((X_angle-x(ix)).^2+(Y-0).^2+(Z_angle-z1(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*k1.*rn)./rn;
        B=sum(sum(A));%��������õ�ֵ�ۼ� 
        pr1(ix,iz)=1i*medium(1).density*u*medium(1).soundspeed*k1/(2*pi)*B; %������ز����õ���ѹp
    end
end

% interface
z_interface1=30e-3;

%rayleigh���ּ��㾭������1������1�ͽ���2֮�䣩���ķ�������
tic
for ix=1:nx  %�۲������x���������
    for iy=1:ny  %�۲������z���������      
        rm1=sqrt((X_angle-x(ix)).^2+(Y-y(iy)).^2+(Z_angle-z_interface1).^2);%interface1����Դ�ľ���
        theta_m1=acos((z_interface1-Z_angle)./rm1);
        theta_m2=asin(medium(2).soundspeed/medium(1).soundspeed.*sin(theta_m1));
        Tm1=2*medium(1).soundspeed*medium(1).density.*cos(theta_m2)./(medium(2).soundspeed*medium(2).density.*cos(theta_m1)+medium(1).soundspeed*medium(1).density.*cos(theta_m2));
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS_m1=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=exp(-j.*k1.*rm1)./rm1.*(1-j./(k1.*rm1)).*abs(Tm1).*(cos(theta_m1)).*dS_m1;
        B=sum(sum(A));%��������õ�ֵ�ۼ� 
        u_m2(ix,iy)=j*u*k1/(2*pi)*B; 
    end
end
toc
% x1=(xmin+dx/2):dx:(xmax-dx/2)
X1=repmat(x',1,ny);
Y1=repmat(y,nx,1);
Z1=repmat(z_interface1,nx,ny);
%�������2�е���ѹ�ֲ�
for ix=1:nx  %�۲������x���������
    for iz=1:nz2  %�۲������z���������      
        rm2=sqrt((X1-x(ix)).^2+(Y1-0).^2+(Z1-z2(iz)).^2);%interface1����Դ�ľ���
        dS_m2=dx*dy; % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        Am2=u_m2.*exp(-j.*k2.*rm2)./rm2*dS_m2;
        Bm2=sum(sum(Am2));%��������õ�ֵ�ۼ� 
        pr2(ix,iz)=j*k2*medium(2).soundspeed*medium(2).density/(2*pi)*Bm2; %������ز����õ���ѹp
    end
end

pr2_abs=abs(pr2);
pr_max=max(abs(pr2(:)));

%ѡ�����Ȥ����
pr_nor=abs(pr2)./pr_max;
[index]=find(pr_nor>=0.5);
pr_zeros_xz(index)=pr_nor(index);
X_all=repmat(x',1,nz2);
Xp(index)=X_all(index);
Z_all=repmat(z2,nx,1);
Zp(index)=Z_all(index);
%������Ȥ������ת�غ�z��ƽ�е�λ�ã��ٽ��бȽ�
Zp_angle=cos(-angle_rotation)*(Zp-Rz)-sin(-angle_rotation)*(Xp-Rx)+Rz; 
Xp_angle=sin(-angle_rotation)*(Zp-Rz)+cos(-angle_rotation)*(Xp-Rx)+Rx;
Xp(index)=Xp_angle(index);%����Ϊ������ת�س����z���غ�
Zp(index)=Zp_angle(index);
axial=max(Zp(:))-min(Zp(:));%������ ��С��0
radial=max(Xp(:))-min(Xp(:));



%�����־���ƴ��
pr=[pr1 pr2];
pr_abs=abs(pr);
z=[z1 z2];

figure(1);
pr_nor=pr_abs./abs(max(pr2(:)));
surf(z*1000,x*1000,pr_nor);
shading interp;
xlabel('z (mm)');
ylabel('x (mm)');
axis equal;
shading flat;
title('Rayleigh ');


