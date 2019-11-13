%rayleigh����ά�ռ��е���ǿ��FOCUS����rayleigh��ά�ռ��е���������ѹת��Ϊ��ǿ���бȽϼ����������������ͼֻѡȡ�˸���Ȥ����
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

%����focus����rayleigh����
xdcr = get_spherical_shell(a,R); %�õ����滻����
delta = [dx dy dz];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%�����Ļ���
ndiv = 80;%���ֵĵ���
dflag = 0;%���=1�������ʲô��ͬ��

tic
pre=fnm_call(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS�Դ���fnm����
I_pre= acousticintensity(pre,medium.density,medium.soundspeed); %��ǿ���� 
toc

%�ҵ�����Ȥ��������ǿ-6dB��Χ�����������Χ�ڵ�error;�ҵ����ֵ�����ض�Ӧ����ά���꣬�������Ӧ��xyƽ�棬�Ѹ���Ȥ���򻭳���
max_index=find_maxpoint(I_pre);%��������λ������
x_index=max_index(1);%���ֵλ��x�±�
y_index=max_index(2);%���ֵλ��y�±�
z_index=max_index(3);%���ֵλ��z�±�
I_pr_max=max(I_pr(:));%�Զ���rayleigh��ǿ���ֵ
I_pre_max=max(I_pre(:));

%������ά�ռ������ַ��������
I_pre_normalized=I_pre./I_pre_max;%����ǿ��һ��
I_pr_normalized=I_pr./I_pr_max;%����ǿ��һ��
error_I=abs(I_pr_normalized-I_pre_normalized)./I_pr_normalized;

%ѡ����ά�ռ����Ȥ�������
[index]=find(I_pr_normalized>=0.25);
error_zeros(index)=error_I(index);

%��������Ȥ�������
figure(1);
surf(error_zeros(:,:,z_index));
axis equal;
shading interp;
colorbar;
xlabel('x');
ylabel('y');
title('rayleigh����3D VS FNM(FOCUS) 3D�ڽ��㴦��xyƽ��������');

figure(2);
error_xz=squeeze(error_zeros(:,y_index,:));
surf(error_xz);
axis equal;
shading interp;
colorbar;
xlabel('z');
ylabel('x');
title('rayleigh����3D VS FNM(FOCUS)3D�ڽ��㴦��xzƽ��������');

figure(3);
surf(I_pre_normalized(:,:,z_index));
axis equal;
shading interp;
colorbar;
xlabel('y');
ylabel('x');
title('rayleigh(FOCUS) 3D�ڽ��㴦��xyƽ��');

figure(4);
surf(I_pr_normalized(:,:,z_index));
axis equal;
shading interp;
colorbar;
xlabel('y');
ylabel('x');
title('rayleigh����3D �ڽ��㴦��xyƽ��');

figure(5);
plot(max(error_xz));
xlabel('z');
ylabel('error');
title('rayleigh VS rayleigh+ASA ����Ȥ��xz����������ı仯���');

figure(6);
histogram(dS0,10);%����dS��ֱ��ͼ�ֲ�

% 
% figure(5);
% plot(max(error_zeros(:,y_index,:)));
% xlabel('z');
% ylabel('error');
% title('rayleigh����3D VS FNM 3D ����Ȥ����1mm��ÿ��������������');

% 
% %����-6dB �����17-62  ÿ��4������1mm
% zz=17:1:61;
% for i=1:length(zz)
%     error_max(i)=max(max(error_zeros(:,:,zz(i))));
%     error_average(i)=mean(mean(error_zeros(:,:,zz(i))));
%     error_median(i)=median(error_zeros(:,:,zz(i)),'all');
% end

% figure(4);
% plot(zz,error_average);
% xlabel('z');
% ylabel('error');
% title('����Ȥ����1mm��ÿ�����ƽ��������');
% figure(5);
% plot(zz,error_median);


