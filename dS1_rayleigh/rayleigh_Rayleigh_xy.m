%ʹ������rayleigh����xy�������ֲ��������������xyƽ��
clc;
clear all;
clear all;
f0=1e6;%����Ƶ�ʺͷ�������
u=1;
medium = set_medium('water');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f

k=2*pi/lambda;%����

R = 75e-3;%ROC���ʰ뾶
a = 30e-3;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���

%���������
xmin=-1.5*a;%�۲������ķ�Χ
xmax=-xmin;
ymin=xmin;
ymax=-ymin;
z0=30e-3;

dx = lambda/6; %�����Ĳ���
dy = lambda/6; 
dz = 0;

x=xmin:dx:xmax;%�����ķֲ�
y=ymin:dy:ymax;
z=z0;

nx=length(x);%�����ĵ���
ny=length(y);
% error_zeros_xz=zeros(nx,nz); %����ȫ��������

% %����focus����rayleigh����rayleigh_cw()
% omega = 2 * pi * f0;
% phi0 = asin(a/R);
% xdcr=get_spherical_shell(a,R);
% 
% delta = [dx dy 0];%����㲽��
% ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, z0, z0);%�����Ļ���
% ndiv = 100;%���ֵĵ���
% dflag = 0;
% 
% tic
% prs=fnm_call(xdcr,ps,medium,ndiv,f0);%FOCUS�Դ���rayleigh����
% toc
% I_prs=acousticintensity(prs,medium.density,medium.soundspeed);
% I_prs_nor=I_prs./max(I_prs(:));

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
Z0=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�

%rayleigh���ּ���xzƽ�������  
tic
for ix=1:nx  %�۲������x���������
    for iy=1:ny  %�۲������z���������
        rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z0).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        pr(ix,iy)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc
% 
% I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
% I_pr_nor=I_pr./max(I_pr(:));
% error=abs(I_pr_nor-I_prs_nor)./I_pr_nor;

% % plot the pressures
% figure(1);
% surf(y*1000,x*1000,I_prs_nor);%��prs��һ��   
% shading interp
% colorbar;
% axis equal;
% xlabel('y(mm)');
% ylabel('x(mm)');
% 
% 
% figure(2);
% surf(y*1000,x*1000,I_pr_nor);%��pr��һ��   
% shading interp;
% axis equal;
% colorbar;
% % title('z=R,Rayleigh Sommerfeld Result');
% xlabel('y(mm)');
% ylabel('x(mm)');
% 
% figure(3);
% surf(y*1000,x*1000,error);%��pr��һ��   
% shading interp;
% axis equal;
% colorbar;
% xlabel('y(mm)');
% ylabel('x(mm)');