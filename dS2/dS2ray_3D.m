%rayleigh����ά�ռ��е���ǿ��FOCUS����rayleigh��ά�ռ��е���������ѹת��Ϊ��ǿ���бȽϼ����������������ͼֻѡȡ�˸���Ȥ����
clc;
clear all;
clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ�����->ˮ���ɸĳɶ�㣬�ú���set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f 1.5mm   lambda/6=0.25mm
k=2*pi/lambda;%����


Roc = 5 * 2 * lambda;%ROC���ʰ뾶
a = 5 * lambda;%ע������a�ǿ׾���һ��
fnumber=Roc/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(Roc^2 - a^2);%���۽��㵽�׾����ĵľ���
phi0=asin(a/Roc);
%���������
xmin=-a;%�۲������ķ�Χ
xmax=-xmin;
ymax=xmax;
ymin=-ymax;
zDiff=0.7*d;
zmin=Roc-zDiff;
zmax=Roc+zDiff;

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
dR=lambda/6;%ÿ�������Ļ������
dphi=dR/Roc;%xzƽ���϶�Ӧ�ĵȽǶ�
nphi=round(phi0/dphi); %�����ĸ�����
phi1=dphi:dphi:phi0;
phi=phi1-dphi/2;%ri��Ӧ��xzƽ��ĽǶ�
ri0=Roc*sin(phi);%ri�ǵ�i����Ӧ�뾶
Sm=lambda/6;%ÿ��Բ����Ӧ��Բ�����Ƶ�ʹ̶�
ntheta=round(2*pi*ri0/Sm);  %��i��Բ����ɢ�Ļ��ȸ���
dtheta=Sm./ri0;  %ÿ��Բ����ɢ�Ļ���
% detheta1=dtheta-dtheta/2;
theta=dthetarepet(dtheta,ntheta);%���һ������
ri=repelem(ri0,ntheta);
Z0=Roc-Roc*cos(phi); %ÿ���뾶Ϊri��Բ����Ӧ��z������
Z=repelem(Z0,ntheta);

% ��y=0ʱxzƽ��
y=0;

%rayleigh���ּ���xzƽ�������  
tic
for ix=1:nx  %�۲������x���������
    for iz=1:nz  %�۲������z���������
        rn=sqrt(ri.^2+x(ix)^2+(Z-z(iz)).^2-2*ri.*x(ix).*cos(theta).*sign(x(ix)));%�۲�㵽��Դ�ľ���
        dS=dR*Sm; % ���Դ���=Բ���Ĳ������dR*��i��������ɢ�Ļ���
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%��������õ�ֵ�ۼ�
        p=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %������ز����õ���ѹp
        pr(ix,iz)=squeeze(p);%�õ�������϶�Ӧ��ֵ��ȡģ
    end
end
I_pr= acousticintensity(pr,medium.density,medium.soundspeed); %��ǿ����
toc

%����focus����rayleigh����
xdcr = get_spherical_shell(a,Roc); %�õ����滻����
delta = [dx dy dz];%����㲽��
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%�����Ļ���
ndiv = 80;%���ֵĵ���
dflag = 0;%���=1�������ʲô��ͬ��

tic
pre=rayleigh_cw(xdcr,ps,medium,ndiv,f0,dflag);%FOCUS�Դ���fnm����
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


