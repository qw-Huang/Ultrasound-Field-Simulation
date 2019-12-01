%ʹ��dS1����Դ���ַ�������2dƽ��xz�棬�ӿ����ο׵Ļ�������ѹ���㣬�Ƚϻ������������֡����ײ��֡�ʣ�ಿ�ֵ���ѹ��ֵ����λ
clc;
close all;
clear all;
f0=1e6;%����Ƶ�ʺ�������
P=100;
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����
hole_ring_back=1e-3;
% hole_ring_back=0;
hole_ring_after=3e-3;
n=1;
R = n * 5 * 2 * lambda;%ROC���ʰ뾶
a = n * 5 * lambda;%ע�������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

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
ndr=round((a-0)/dr);
dr=(a-0)/ndr;
r_back=0:dr:a-dr;%����Դǰһ�λ�����Ӧ��r
r_after=dr:dr:a;%����Դ��һ�λ�����Ӧ��r
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

%���廻������������
u_ring_res=repmat(u,ntheta,ndr);
back_index=find(r_back<=hole_ring_back,1,'last');%�жϻ����ھ������ĸ���Χ��
after_index=find(r_after>=hole_ring_after,1,'first');%�жϻ����⾶�����ĸ���Χ��
if back_index==after_index
    index_ring=back_index;
else index_ring=back_index:after_index;
end
u_ring_res(:,index_ring)=[0]; %�����ײ��ַ���������Ϊ0

% ��y=0ʱxzƽ��
y=0;

%�󿪿�Ϊ����ʣ�ಿ�ֵ���ѹ
tic
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*k.*rn)./rn.*u_ring_res;
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        pr_ring_res(ix,iz)=1i*medium.density*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc

%�󿪿ײ��ֵ���ѹ
%���忪�ײ��ַ�������
u_ring_hole=zeros(ntheta,ndr);
u_ring_hole(:,index_ring)=[u]; %�����ײ��ַ���������Ϊ0
tic
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*k.*rn)./rn.*u_ring_hole;
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        pr_ring_hole(ix,iz)=1i*medium.density*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc

%����������������ѹ
tic
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������
        rn=sqrt((X-x(ix)).^2+(Y-y).^2+(Z-z(iz)).^2);%�۲�㵽��Դ�ľ���
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%ÿһ������ɢ��dS�Ĵ�С����i����ɢ�ĵ���Դ���dS=��i������ɢ����*dr��Ӧ�Ķ̻�
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )���������������ڴ�ֱ������m�Σ���ˮƽ������n�Ρ�
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        pr_all(ix,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
    end
end
toc

%�Ƚ���ѹ�ķ�ֵ����λ
pr_all_amplitude=abs(pr_all);
pr_all_phase=angle(pr_all);
pr_sum=pr_ring_hole+pr_ring_res;
pr_sum_amplitude=abs(pr_sum);
pr_sum_phase=angle(pr_sum);
error_phase=abs(pr_sum_phase-pr_all_phase)./pr_all_phase;
error_amplitude=abs(pr_sum_amplitude-pr_all_amplitude)./pr_all_amplitude;
I_all=abs(pr_all).^2./(2*medium.soundspeed*medium.density);
I_hole=abs(pr_ring_hole).^2./(2*medium.soundspeed*medium.density);
I_res=abs(pr_ring_res).^2./(2*medium.soundspeed*medium.density);
max_hole=find_maxpoint(I_hole);
max_res=find_maxpoint(I_res);
max_all=find_maxpoint(I_all);


% %��ͼ
figure(1);
surf(z*1000,x*1000,abs(pr_all));
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the pressure amplitude of all(R=15mm,a=7.5mm,hole_a=1mm-3mm)');
% hold on;
figure(2);
surf(z*1000,x*1000,abs(pr_ring_res));
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the pressure amplitude of residue(R=15mm,a=7.5mm,hole_a=1mm-3mm)');
figure(3);
surf(z*1000,x*1000,abs(pr_ring_hole));
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the pressure amplitude of ring(R=15mm,a=7.5mm,hole_a=1mm-3mm)');
 figure(4);
surf(z*1000, x*1000,error_amplitude); 
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the error of amplitude(R=15mm,a=7.5mm,hole_a=1mm-3mm)');
figure(5);
surf(z*1000, x*1000,error_phase);
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the error of phase(R=15mm,a=7.5mm,hole_a=1mm-3mm))');
