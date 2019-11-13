%���Լ���ʱ��
clc;
clear all;
profile on;

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
xmax=0;
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

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
nr=round(a/dr);
dr=a/nr;
r_back=0:dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=dr:dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda/6;%����㻷����2pi*a����ɢ�����Ӧ��dS�Ļ���С��lambda/6
median=length(r)/2+1;
ntheta=round(2*pi*r(median)/Sm);%�м��һ�������Ļ��ֵ���,ȡ��
dtheta=2*pi./ntheta;%����ȡ�����µ���ÿ������Դ�Ķ�Ӧ����
theta_after=dtheta:dtheta:2*pi;%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta_back=0:dtheta:(2*pi-dtheta);%ÿ��������ɢ�ɶ�����Ӧ�Ļ�������
theta=theta_after-dtheta/2;%��i��������ɢ�ɵ���Դ��dS�м�ĵ��Ӧ�Ļ�������
X_source=cos(theta)'*r;%��Դ��x����  ������ʽ ntheta*nr
X_row=X_source(:)';%����Դ����ά�ռ�ĺ����꣬ÿ�д�������ÿ�д���һ����ɢ�ĵ���������չ��չ��һ�У�һ���ǰ����ڵ���Ļ����ſ�
Y_source=sin(theta)'*r;%��Դ��y����
Y_row=Y_source(:)';%ÿ�д�������ÿ�д���һ����ɢ�ĵ���������չ��չ��һ�У�һ���ǰ����ڵ���Ļ����ſ�
Z0=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z_source=repmat(Z0,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
Z_row=Z_source(:)';
% ��y=0ʱxzƽ��
y=0;

%��X Y Z�ֱ���չ
X=repmat(X_row,nx,1);
Y=repmat(Y_row,nx,1);
Z=repmat(Z_row,nx,1);

x_colnum=x';
x_repeat=repmat(x_colnum,1,ntheta*nr);
%rayleigh���ּ���xzƽ�������  
tic

    for iz=1:nz  %�۲������z���������
        rn=sqrt((X-x_repeat).^2+(Y-y).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS_source=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        dS_row=dS_source(:)';
        dS=repmat(dS_row,nx,1);
        A=dS.*exp(-1i*k.*rn)./rn;
        
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        pr_up(:,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end

toc
pr_down=flipud(pr_up);%����ѹ������з�ת���൱�ڹ���x=0��Գ�
pr_all=[pr_up;pr_down];%��������������ƴ�������������м�һ�лḴ������
pr=pr_all;
[row,column]=size(pr);
index_median=row/2;%ȡ�м���
pr(index_median,:)=[]; %ɾ���м��ظ���һ��

%��ѹת��Ϊ��������
I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
I_pr_nor=I_pr./max(I_pr(:));

profile viewer
