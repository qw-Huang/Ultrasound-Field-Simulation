%ʹ��dS1����Դ���ַ�������2dƽ��xz�棬�ӿ��׵Ļ�������ǿ����
clear all;
f0=1e6;%����Ƶ�ʺ�������
P=100;
n=8;
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����
R = n * 5 * 2 * lambda;%ROC���ʰ뾶
a = n * 5 * lambda;%ע�������a�ǿ׾���һ��
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

% ��y=0ʱxzƽ��
y=0;

error_zeros_xz=zeros(nx,nz); %����ȫ��������
I_pr_zeros=zeros(nx,nz);
I_pr_hole_zeros=zeros(nx,nz);
%�ӿ���
hole_a=n*7.5e-4;
u_hole=normal_velocity(P,R,a,hole_a,medium.density,medium.soundspeed);

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
ndr=round((a-dr-hole_a)/dr);
dr=(a-dr-hole_a)/ndr;
r_back=(0+hole_a):dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=(dr+hole_a):dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda/6;%����㻷����2pi*a����ɢ�����Ӧ��dS�Ļ���С��lambda/6
median=fix(length(r)/2)+1;%fix��λȡ��
ntheta=round(2*pi*r(median)/Sm);%�����һ�������Ļ��ֵ���,ȡ��
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
    for iz=1:nz  %�۲������z���������
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        pr_hole(ix,iz)=1i*medium.density*u_hole*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc
%��ѹת��Ϊ��ǿ����
I_pr_hole=acousticintensity(pr_hole,medium.density,medium.soundspeed); 

%���ӿ���
hole_a=0;
u=normal_velocity(P,R,a,hole_a,medium.density,medium.soundspeed);

%��������ɢ��Ϊ����Դ  ע����ɢ��֮�����Դ�Ĵ�С
dr=lambda/6;%����ָ�ɺܶ��Բ������Ӧ�İ뾶��r��dr�ǰ뾶���ӵĲ���
ndr=round((a-dr-hole_a)/dr);
dr=(a-dr-hole_a)/ndr;
r_back=(0+hole_a):dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=(dr+hole_a):dr:a;%����Դ��һ�λ�����Ӧ��r
r=r_after-dr/2;%��i��������Ӧ�����ĵ��r
Sm=lambda/6;%����㻷����2pi*a����ɢ�����Ӧ��dS�Ļ���С��lambda/6
median=fix(length(r)/2)+1;%fix��λȡ��
ntheta=round(2*pi*r(median)/Sm);%�����һ�������Ļ��ֵ���,ȡ��
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
%��ѹת��Ϊ��ǿ����
I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 

%��ǿ���� ��Ӧ������λ�� -6dB��С ǰ�����
%�ҵ�����Ȥ��������ǿ-6dB��Χ�����������Χ�ڵ�error;�ҵ����ֵ�����ض�Ӧ����ά���꣬�������Ӧ��xyƽ�棬�Ѹ���Ȥ���򻭳���
max_index_hole=find_maxpoint(I_pr_hole);%��������λ������
max_index=find_maxpoint(I_pr);
z_index_hole=max_index_hole(3);%���ֵλ��z�±�
z_index=max_index(3);%���ֵλ��z�±�
focusz_change=(z_index_hole-z_index)*dz;%�Ա����޿��׵Ľ���ǰ����� �����Ŵ���ǰ�ƣ�
I_pr_max=max(I_pr(:));%�޿�����ǿ���ֵ
I_pr_hole_max=max(I_pr_hole(:));%�п�����ǿ����
focus_Ichange=(I_pr_max-I_pr_hole_max)/I_pr_max; %������ǿ�仯

%ѡ�����Ȥ�������
[index]=find(I_pr>=0.25*I_pr_max);
error_I=abs(I_pr-I_pr_hole)./I_pr;
error_zeros_xz(index)=error_I(index);
[index_hole]=find(I_pr_hole>=0.25*I_pr_hole_max);
I_pr_zeros(index)=I_pr(index); %�޿��׵�-6dB��ǿ�ֲ�
I_pr_hole_zeros(index_hole)=I_pr_hole(index_hole);%�п��׵�-6dB��ǿ�ֲ�

%xy��������ֵ����Զ��������������Ĳ������Բ�������-6dB��Χ�ľ���  
range_point_xz_hole=find_range_point(I_pr_hole_zeros,nx,nz);%���ܴ�������Խ�磬����ڱ߽紦���ݲ�Ϊ0�������ͻ��������
range_point_xz=find_range_point(I_pr_zeros,nx,nz);
range_value_xz=range_point_xz*dz;%��ֵ=����*����
range_value_xz_hole=range_point_xz_hole*dz;%��ֵ=����*����


% %��ͼ
% figure(1);
% surf(z*1000,x*1000,zeros_I_pr_hole); %������Ϊz��������Ϊx����ɫ����p1�Ĵ�С  
% shading interp
% title('Rayleigh(dS1 hole_a=7.5e-4)  ');
% axis equal;%��������ı�����ͬ
% colorbar
% xlabel('z��mm�� ');
% ylabel('x (mm) ');
% zlabel('normalized pressure');
% 
% figure(2);
% surf(z*1000, x*1000, zeros_I_pr);%��prs��һ��   
% shading interp;
% colorbar;
% axis equal;
% title('Rayleigh');
% xlabel('z (mm) ');
% ylabel('x (mm) ');
% zlabel('normalized pressure');

% �Ա����ַ����������Ľ��
figure(3);
surf(z*1000, x*1000, error_zeros_xz); 
shading interp %ȥ������ƽ������
axis equal;
colorbar  %����ɫ��
xlabel('z (mm) ');
ylabel('x (mm) ');
title('error between hole and no hole(R=120mm,a=60mm,10%a)');

% figure(4);
% histogram(dS_ring,10);
% title('��Դ���dS�ֲ�');
% 
% figure(5);
% plot(error_dSi);
% xlabel('ith ring ');
% ylabel('error ');
% title('error between the sum of dS and theory calculation');
