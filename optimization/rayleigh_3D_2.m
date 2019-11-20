%����ʱ�� ����ά�Ĵ�����иĽ���ֻ�����ķ�֮һ������ ��δ�ģ�
clc;
clear all;
clear all;
f0=1e6;%����Ƶ�ʺͷ�������
P=100;
medium = set_medium('lossless');%������ʣ�����->ˮ���ɸĳɶ�㣬�ú���set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����

R = 5 * 2 * lambda;%ROC���ʰ뾶
a = 5 * lambda;%ע������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,medium.density,medium.soundspeed);

%���������
xmin=-a;%�۲������ķ�Χ
xmax=-xmin;
ymin=xmin;
ymax=-ymin;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %�����Ĳ���
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;%�����ķֲ�
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);%�����ĵ���
ny=length(y);
nz=length(z);

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
Z_ring=R-sqrt(R*R-r.*r);%����Դ��ά�ռ��е�z����
Z=repmat(Z_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�

%rayleigh���ּ�����ά�ռ������  
tic
for ix=1:nx  %�۲������x��������껮��
    for iy=1:ny %�۲������y��������껮��
       for iz=1:nz  %�۲������z��������껮��
            rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
            dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
            dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
            A=dS.*exp(-1i.*k.*rn)./rn;
            B=sum(sum(A));%��������õ�ֵ�ۼ�
            pr(ix,iy,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp   
        end
    end
end
toc
pr_max=max(pr(:));
I_pr=abs(pr).^2/(medium.density*medium.soundspeed);

%�ҵ�����Ȥ��������ǿ-6dB��Χ;�ҵ����ֵ�����ض�Ӧ����ά���꣬�������Ӧ��xyƽ�棬�Ѹ���Ȥ���򻭳���
max_index=find_maxpoint(I_pr);%��������λ������
x_index=max_index(1);%���ֵλ��x�±�
y_index=max_index(2);%���ֵλ��y�±�
z_index=max_index(3);%���ֵλ��z�±�
I_pr_max=max(I_pr(:));%�Զ���rayleigh��ǿ���ֵ
I_pr_xz=I_pr(:,y_index,:)./I_pr_max;
I_pr_xy=I_pr(:,:,z_index)./I_pr_max;

figure(1);
surf(squeeze(I_pr_xz));%�����xzƽ�棬(z,x,p)
axis equal;
shading interp;
colorbar;
title('�Զ���Rayleigh����xzƽ����ǿ�ֲ� ');
xlabel('z ');
ylabel('x ');

figure(2);
surf(squeeze(I_pr_xy));%�����xzƽ�棬(z,x,p)
axis equal;
shading interp;
colorbar;
title('�Զ���Rayleigh����xyƽ����ǿ�ֲ� ');
xlabel('y ');
ylabel('x ');

