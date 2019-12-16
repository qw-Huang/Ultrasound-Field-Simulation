function [output] = layered_hole_rayleigh_2D(R,a,hole_a,P,medium,x,z)
%LAYERED_RAYLEIGH_2D ����������ʵ���ѹ�ֲ���ˮ+���⣩
%   �˴���ʾ��ϸ˵��
j=sqrt(-1);
f0=1e6;
lambda1 = medium(1).soundspeed/f0;%����=c/f
lambda2 = medium(2).soundspeed/f0;%����=c/f
k1=2*pi/lambda1-j*medium(1).attenuationdBcmMHz;%����
k2=2*pi/lambda2-j*medium(2).attenuationdBcmMHz;%����

fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,hole_a,medium(1).density,medium(1).soundspeed);

dx1=lambda1/6;
dy1=lambda1/6;
x1=-a:dx1:a;
y1=x1;

nx1=length(x1);%�����ĵ���
ny1=length(y1);
nx=length(x);
nz=length(z);

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda1/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
ndr=round((a-hole_a)/dr);
dr=(a-hole_a)/ndr;
r_back=(0+hole_a):dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=(dr+hole_a):dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda1/6;%�м价����ɢ�����Ӧ��dS�Ļ�������lambda/6
median=round(length(r)/2);%ȡ�м价��������
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
z_interface1=5e-3;

%rayleigh���ּ��㾭������1������1�ͽ���2֮�䣩���ķ�������
tic
for ix=1:nx1  %�۲������x���������
    for iy=1:ny1  %�۲������z���������      
        rm1=sqrt((X-x1(ix)).^2+(Y-y1(iy)).^2+(Z-z_interface1).^2);%interface1����Դ�ľ���
        theta_m1=acos((Z-z_interface1)./rm1);
        theta_m2=asin(medium(2).soundspeed/medium(1).soundspeed*sin(theta_m1));
        Tm1=2*medium(2).soundspeed*medium(2).density.*cos(theta_m1)./(medium(2).soundspeed*medium(2).density.*cos(theta_m1)+medium(1).soundspeed*medium(1).density.*cos(theta_m2));
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS_m1=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS_m1.*exp(-j.*k1.*rm1)./rm1.*(1-j./(k1.*rm1)).*Tm1.*cos(theta_m1);
        B=sum(sum(A));%��������õ�ֵ�ۼ� 
        u_m2(ix,iy)=j*u*k1/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc

X1=repmat(x1',1,ny1);
Y1=repmat(y1,nx1,1);
Z1=repmat(z_interface1,nx1,ny1);
%�������2�е���ѹ�ֲ�
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������      
        rm2=sqrt((X1-x(ix)).^2+(Y1-0).^2+(Z1-z(iz)).^2);%interface1����Դ�ľ���
        dS_m2=dx1*dy1; % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        Am2=u_m2.*exp(-j.*k2.*rm2)./rm2*dS_m2;
        Bm2=sum(sum(Am2));%��������õ�ֵ�ۼ� 
        pr(ix,iz)=j*k2*medium(2).soundspeed*medium(2).density/(2*pi)*Bm2; %������ز����õ���ѹp
    end
end

output = pr;

end

