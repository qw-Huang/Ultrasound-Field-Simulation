%ʹ���Զ���rayleigh��������ά�ռ��е���ǿ�ֲ�,�Ա��Զ���rayleigh���ֺ�focus����rayleigh����+ASA�ĸ���Ȥ�����  
clc;
clear all;
clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ�����->ˮ���ɸĳɶ�㣬�ú���set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%����

R = 5 * 2 * lambda;%ROC���ʰ뾶
a = 5 * lambda;%ע������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���

%���������
xmin=-a;%�۲������ķ�Χ
xmax=-xmin;
ymax=xmax;
ymin=-ymax;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %�����Ĳ���
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);
ny=length(y);
nz=length(z);

error_zeros=zeros(nx,ny,nz);
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

%�Զ���rayleigh���ּ�����ά�ռ������  
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
I_pr= acousticintensity(pr,medium.density,medium.soundspeed); %��ǿ����
toc

% FNM���ASA����
tic
xdcr = get_spherical_shell(a,R); %�õ����滻����
delta = [dx dy 1];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmin);%�����Ļ���
ndiv = 100;%���ֵĵ���
dflag = 0;%���=1�������ʲô��ͬ��
pre=rayleigh_cw(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS�Դ���fnm����
z0=zmin;  %��z=z0ƽ��������ֲ�
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);%ASA����ά�ռ�����������㻮��
pre_asa= cw_angular_spectrum(pre,cg_3d,medium,f0,1024,'Pa');%��������ASA��������ά�����ֲ�
I_pre_asa=acousticintensity(pre_asa,medium.density,medium.soundspeed); %��ǿ����
toc


%�ҵ�����Ȥ��������ǿ-6dB��Χ�����������Χ�ڵ�error;�ҵ����ֵ�����ض�Ӧ����ά���꣬�������Ӧ��xyƽ�棬�Ѹ���Ȥ���򻭳���
max_index=find_maxpoint(I_pr);%��������λ������
x_index=max_index(1);%���ֵλ��x�±�
y_index=max_index(2);%���ֵλ��y�±�
z_index=max_index(3);%���ֵλ��z�±�
I_pre_max=max(I_pre_asa(:));%�Զ���rayleigh 3D��ǿ���ֵ
I_pr_max=max(max(I_pr(:)));

I_pr_nor=I_pr./I_pr_max;
I_pre_nor=I_pre_asa./I_pre_max;
error_I=abs(I_pre_nor-I_pr_nor)./I_pre_nor;%������ά�ռ������ַ��������

%ѡ����ά�ռ����Ȥ�������
[index]=find(I_pre_nor>=0.25);
error_zeros(index)=error_I(index);

%��������Ȥ�������
figure(1);
surf(error_zeros(:,:,z_index));
axis equal;
shading interp;
colorbar
xlabel('y');
ylabel('x');
title('rayleigh 3D VS Rayleigh(FOCUS)+ASA ���㴦xyƽ��������');

figure(2);
error_xz=squeeze(error_zeros(:,y_index,:));
axis equal;
surf(error_xz);
shading interp;
colorbar
xlabel('z');
ylabel('x');
title('rayleigh 3D VS Rayleigh(FOCUS)+ASA ���㴦xzƽ��������');

figure(3);
surf(I_pre_nor(:,:,z_index));
axis equal;
shading interp;
colorbar;

figure(4);
surf(I_pr_nor(:,:,z_index));
axis equal;
shading interp;
colorbar;

figure(5);
surf(squeeze(I_pre_nor(:,y_index,:)));
axis equal;
shading interp 
colorbar;
figure(6);
surf(squeeze(I_pr_nor(:,y_index,:)));
axis equal;
shading interp 
colorbar;
% zz=17:1:61;
% for i=1:length(zz)
%     error_max(i)=max(max(error_zeros(:,:,zz(i))));
%     error_average(i)=mean(mean(error_zeros(:,:,zz(i))));
% end
% 
% figure(6);
% plot(zz,error_max);
% xlabel('z');
% ylabel('error');
% title('rayleigh 3D VS FNM+ASA ����Ȥ��xy����������ı仯���');
% figure(7);
% plot(zz,error_average);
% xlabel('z');
% ylabel('error');
% title('rayleigh 3D VS FNM+ASA ����Ȥ��xy����ƽ�����ı仯���');


