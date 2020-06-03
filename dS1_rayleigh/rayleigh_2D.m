%ʹ��dS1����Դ���ַ�������2dƽ��xz�棬��һ���Ա��Զ���rayleigh���ֺ�focus�Դ�
%��rayleigh���֣��Ա���ǿ�ֲ������
% clc;
clear all;
f0=0.8e6;%����Ƶ�ʺͷ�������
P=100;
medium = set_medium('muscle');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
dBperNeper = 20 * log10(exp(1));
% dBperNeper =8.6886;
attenuationNeperspermeter=medium.attenuationdBcmMHz/dBperNeper*100*(f0/1e6)^1.1;
kk=2*pi/lambda-j*attenuationNeperspermeter;%����
k=2*pi/lambda;%����

% k=2*pi/lambda-j*medium.attenuationdBcmMHz;%����

R = 150e-3;%ROC���ʰ뾶
a=60e-3;
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);
% u=1;

%���������
xmin=-1.5*a;%�۲������ķ�Χ
xmax=-xmin;
ymax=0;
ymin=0;

zmin=140e-3;
% zmin=R-zDiff;
zmax=150e-3;

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

% ��y=0ʱxzƽ��
y=0;

%rayleigh���ּ���xzƽ�������  
tic
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*kk.*rn)./rn;
        B=sum(sum(A));%��������õ�ֵ�ۼ� 
        pr(ix,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc

pr_max=max(abs(pr(:)));
pr_abs=abs(pr);
focus_index_1=find(pr_abs==pr_max);
s=size(pr_abs);
[x_index,z_index]=ind2sub(s,focus_index_1);%�����ֵ���±�תΪ��ά���±�
focus_index=z_index;
% 
% %�Ƚ�dS���ۼӺͽ������
% S_discrete=sum(sum(dS));%��ɢ֮���ȫ������ۼ�
% S_theory=2*pi*R*(R-d);%���ۼ���������������������
% dSi=sum(dS);%��ÿһ��������ۼӣ�Ҳ���ǰ�ÿһ�У�70�������ۼ�
% di_back=sqrt(R^2-r_back.^2);%��0��ʼ��Ϊ��i�㻷��d
% di_after=sqrt(R^2-r_after.^2);%��dr��ʼ��Ϊ��i�㻷��d
% dSi_theory=2*pi*R*(R-di_after)-2*pi*R*(R-di_back);%��i���Ļ������
% error_dSi=abs(dSi-dSi_theory)./dSi;%������⻷���������ֵ�⻷����������
% 
% %����focus����rayleigh����rayleigh_cw()
% omega = 2 * pi * f0;
% phi0 = asin(a/R);
% xdcr = get_spherical_shell(a,R); %�õ����滻����
% delta = [dx 0 dz];%����㲽��
% ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%�����Ļ���
% ndiv = 100;%���ֵĵ���
% dflag = 0;%���=1�������ʲô��ͬ��
% 
% tic
% % prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0);%FOCUS�Դ���rayleigh����
% prs=cw_pressure(xdcr,ps,medium,ndiv,f0);%FOCUS�Դ���fnm����
% 
% toc
% 
% %��ѹת��Ϊ��������
% I_prs=acousticintensity(prs,medium.density,medium.soundspeed);
% I_prs_nor=squeeze(I_prs./max(I_prs(:)));
% I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
% I_pr_nor=I_pr./max(I_pr(:));
% error=abs(I_pr_nor-I_prs_nor)./I_pr_nor;
% 
% %ѡ�����Ȥ�������
% [index]=find(I_prs_nor>=0.25);
% error_zeros_xz(index)=error(index);

% figure(3)
% surf(z*1000,x*1000,abs(pr));
% shading interp;
% colorbar;
% axis equal;
