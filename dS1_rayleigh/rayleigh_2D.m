%ʹ��dS1����Դ���ַ�������2dƽ��xz�棬��һ���Ա��Զ���rayleigh���ֺ�focus�Դ�
%��rayleigh���֣��Ա���ǿ�ֲ������
clc;
clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����

R = 1 * 5 * 2 * lambda;%ROC���ʰ뾶
a = 1 * 5 * lambda;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���

%���������
xmin=-a;%�۲������ķ�Χ
xmax=-xmin;
ymax=0;
ymin=0;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %�����Ĳ���
dz = lambda/6;

x=xmin:dx:xmax;%�����ķֲ�
z=zmin:dz:zmax;

nx=length(x);%�����ĵ���
nz=length(z);

error_zeros_xz=zeros(nx,nz); %����ȫ��������

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
r_back=0:dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=dr:dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda/6;%�м价����ɢ�����Ӧ��dS�Ļ�������lambda/6
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

% ��y=0ʱxzƽ��
y=0;

%rayleigh���ּ���xzƽ�������  
tic
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        pr(ix,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc

% 
% %�Ƚ�dS���ۼӺͽ������
% S_discrete=sum(sum(dS));%��ɢ֮���ȫ������ۼ�
% S_theory=2*pi*R*(R-d);%���ۼ���������������������
% dSi=sum(dS);%��ÿһ��������ۼӣ�Ҳ���ǰ�ÿһ�У�70�������ۼ�
% di_back=sqrt(R^2-r_back.^2);%��0��ʼ��Ϊ��i�㻷��d
% di_after=sqrt(R^2-r_after.^2);%��dr��ʼ��Ϊ��i�㻷��d
% dSi_theory=2*pi*R*(R-di_after)-2*pi*R*(R-di_back);%��i���Ļ������
% error_dSi=abs(dSi-dSi_theory)./dSi;%������⻷���������ֵ�⻷����������

%����focus����rayleigh����rayleigh_cw()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr = get_spherical_shell(a,R); %�õ����滻����
delta = [dx 0 dz];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%�����Ļ���
ndiv = 100;%���ֵĵ���
dflag = 0;%���=1�������ʲô��ͬ��

tic
prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0);%FOCUS�Դ���rayleigh����
toc

%��ѹת��Ϊ��������
I_prs=acousticintensity(prs,medium.density,medium.soundspeed);
I_prs_nor=squeeze(I_prs./max(I_prs(:)));
I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
I_pr_nor=I_pr./max(I_pr(:));
error=abs(I_pr_nor-I_prs_nor)./I_pr_nor;

%ѡ�����Ȥ�������
[index]=find(I_prs_nor>=0.25);
error_zeros_xz(index)=error(index);

%��ͼ
figure(1);
surf(z*1000,x*1000,I_pr_nor); %������Ϊz��������Ϊx����ɫ����p1�Ĵ�С  
shading interp
title('Rayleigh(dS1)  ');
axis equal;%��������ı�����ͬ
colorbar
xlabel('z��mm�� ');
ylabel('x (mm) ');
zlabel('normalized pressure');

figure(2);
surf(z*1000, x*1000, I_prs_nor);%��prs��һ��   
shading interp;
colorbar;
axis equal;
title('Rayleigh(FOCUS)');
xlabel('z (mm) ');
ylabel('x (mm) ');
zlabel('normalized pressure');

%�Ա����ַ����������Ľ��
figure(3);
surf(z*1000, x*1000, error_zeros_xz); 
shading interp %ȥ������ƽ������
axis equal;
colorbar  %����ɫ��
xlabel('z (mm) ');
ylabel('x (mm) ');
title('error between Rayleigh and Rayleigh (FOCUS)');

figure(4);
histogram(dS_ring,10);
title('��Դ���dS�ֲ�');

figure(5);
plot(error_dSi);
xlabel('ith ring ');
ylabel('error ');
title('error between the sum of dS and theory calculation');