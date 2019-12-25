%�Ľ���rayleigh���֣���������������еķֲ� ˮ5mm+����10mm-20mm
%ROC =15mm

% clc;
clear all;
j=1i;
f0=1e6;
% P=100;%����Ƶ�ʺͷ�������
% medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
medium = set_layered_medium([0,30e-3],[set_medium('fat'),set_medium('fat')]);
lambda1 = medium(1).soundspeed/f0;%����=c/f
lambda2 = medium(2).soundspeed/f0;%����=c/f
% k1=2*pi*f0/medium(1).soundspeed-j*medium(1).attenuationdBcmMHz;%����
% k2=2*pi*f0/medium(2).soundspeed-j*medium(2).attenuationdBcmMHz;%����

dBperNeper = 20 * log10(exp(1));
attenuationNeperspermeter1=medium(1).attenuationdBcmMHz/dBperNeper*100*f0/1e6;
k1=2*pi/lambda1-j*attenuationNeperspermeter1;%����
attenuationNeperspermeter2=medium(2).attenuationdBcmMHz/dBperNeper*100*f0/1e6;
k2=2*pi/lambda2-j*attenuationNeperspermeter2;%����

R = 75e-3;%ROC���ʰ뾶
a = 30e-3;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
% u=normal_velocity(P,R,a,0,medium(1).density,medium(1).soundspeed);
u=1;
%���������
xmin=-1.5*a;%�۲������ķ�Χ
xmax=-xmin;
ymax=1.5*a;
ymin=-ymax;

zmin=70e-3;
zmax=75e-3;

dx = lambda1/6; %�����Ĳ���
dy = lambda1/6;
dz = lambda1/6;


x=xmin:dx:xmax;%�����ķֲ�
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);%�����ĵ���
ny=length(y);
nz=length(z);


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

% interface
z_interface1=30e-3;

%rayleigh���ּ��㾭������1������1�ͽ���2֮�䣩���ķ�������
tic
for ix=1:nx  %�۲������x���������
    for iy=1:ny  %�۲������z���������      
        rm1=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z_interface1).^2);%interface1����Դ�ľ���
        theta_m1=acos((z_interface1-Z)./rm1);
        theta_m2=asin(medium(2).soundspeed/medium(1).soundspeed.*sin(theta_m1));
%         theta_m22=pi-theta_m21;
%         if abs(theta_m22-theta_m1)<=abs(theta_m21-theta_m1)
%             theta_m2=theta_m22;
%         else theta_m2=theta_m21;
%         end
        Tm1=2*medium(2).soundspeed*medium(2).density.*cos(theta_m1)./(medium(2).soundspeed*medium(2).density.*cos(theta_m1)+medium(1).soundspeed*medium(1).density.*cos(theta_m2));
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS_m1=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=exp(-j.*k1.*rm1)./rm1.*(1-j./(k1.*rm1)).*Tm1.*(cos(theta_m1)).*dS_m1;
        B=sum(sum(A));%��������õ�ֵ�ۼ� 
        u_m2(ix,iy)=j*u*k1/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc

X1=repmat(x',1,ny);
Y1=repmat(y,nx,1);
Z1=repmat(z_interface1,nx,ny);
%�������2�е���ѹ�ֲ�
% for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������      
        rm2=sqrt((X1-0).^2+(Y1-0).^2+(Z1-z(iz)).^2);%interface1����Դ�ľ���
        dS_m2=dx*dy; % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        Am2=u_m2.*exp(-j.*k2.*rm2)./rm2*dS_m2;
        Bm2=sum(sum(Am2));%��������õ�ֵ�ۼ� 
        pr(ix,iz)=j*k2*medium(2).soundspeed*medium(2).density/(2*pi)*Bm2; %������ز����õ���ѹp
    end
% end
figure(1);
pr_nor=abs(pr)./abs(max(pr(:)));
surf(z*1000,x*1000,pr_nor);
shading interp;
