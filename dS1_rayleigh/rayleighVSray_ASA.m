%ʹ���Զ���rayleigh��������ά�ռ��е���ǿ�ֲ�,�Ա��Զ���rayleigh���ֺ��Զ���rayleigh����+ASA�ĸ���Ȥ�����  
clc;
clear all;
clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ�����->ˮ���ɸĳɶ�㣬�ú���set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
k=2*pi/lambda;%����

R = 5 *1* 2 * lambda;%ROC���ʰ뾶
a = 5 * lambda;%ע������a�ǿ׾���һ��
fnumber=R/(2*a);%����f-number=���ʰ뾶/�׾���2*a��
d = sqrt(R^2 - a^2);%���۽��㵽�׾����ĵľ���
%���������
xmin=-0.7*a;%�۲������ķ�Χ
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
toc 

I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
I_asa=rayleigh_ASA(R,a,u,f0);%�Զ���rayleigh+ASA����ǿ
error_I=abs(I_pr-I_asa)./I_asa;%������ά�ռ������ַ��������

%�ҵ�����Ȥ��������ǿ-6dB��Χ�����������Χ�ڵ�error;�ҵ����ֵ�����ض�Ӧ����ά���꣬�������Ӧ��xyƽ�棬�Ѹ���Ȥ���򻭳���
max_index=find_maxpoint(I_asa);%��������λ������
max_index2=find_maxpoint(I_pr);
x_index=max_index(1);%���ֵλ��x�±�
y_index=max_index(2);%���ֵλ��y�±�
z_index=max_index(3);%���ֵλ��z�±�

I_pr_max=max(I_pr(:));%�Զ���rayleigh 3D��ǿ���ֵ
I_asa_max=max(I_asa(:));

%ѡ����ά�ռ����Ȥ�������
[index]=find(I_pr>=0.25*I_pr_max);
error_zeros(index)=error_I(index);

%��������Ȥ�������
figure(1);
surf(y*1000,x*1000,error_zeros(:,:,z_index));
shading interp;
axis equal;
xlabel('y (mm) ');
ylabel('x (mm) ');
title('rayleigh VS rayleigh+ASA ���㴦xyƽ��������');

figure(2);
error_xz=squeeze(error_zeros(:,y_index,:));
surf(error_xz);
axis equal;
shading interp;
xlabel('z');
ylabel('x');
title('rayleigh VS rayleigh+ASA ���㴦xzƽ��������');

figure(3);
I_pr_xz=finddB(I_pr(:,y_index,:),nx,nz);
pcolor(I_pr_xz);
axis equal;
shading interp; 
colorbar;
xlabel('z');
ylabel('x');
title('rayleigh ���㴦xzƽ����ǿ�ֲ����');

figure(4);
I_asa_xz=finddB(I_asa(:,y_index,:),nx,nz);
pcolor(I_asa_xz);
axis equal
shading interp; 
colorbar;
xlabel('z');
ylabel('x');
title('rayleigh+ASA ���㴦xzƽ����ǿ�ֲ����');

figure(5);
plot(max(error_xz));
xlabel('z');
ylabel('error');
title('rayleigh VS rayleigh+ASA ����Ȥ��xz����������ı仯���');
% figure(7);
% plot(mean(error_xz));
% xlabel('z');
% ylabel('error');
% title('rayleigh VS rayleigh+ASA ����Ȥ��xy����ƽ�����ı仯���');

