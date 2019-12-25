%改进的rayleigh积分，计算在两层介质中的分布 水30mm+肌肉30mm-90mm
% clc;
clear all;
j=1i;
f0=1e6;
P=100;%定义频率和法向阵速
medium = set_layered_medium([0,30e-3,35e-3],[set_medium('water'),set_medium('fat'),set_medium('muscle')]);
lambda1 = medium(1).soundspeed/f0;%波长=c/f
lambda2 = medium(2).soundspeed/f0;%波长=c/f
lambda3 = medium(3).soundspeed/f0;%波长=c/f

dBperNeper = 20 * log10(exp(1));
attenuationNeperspermeter1=medium(1).attenuationdBcmMHz/dBperNeper*100*f0/1e6;
k1=2*pi/lambda1-j*attenuationNeperspermeter1;%波数
attenuationNeperspermeter2=medium(2).attenuationdBcmMHz/dBperNeper*100*f0/1e6;
k2=2*pi/lambda2-j*attenuationNeperspermeter2;%波数
attenuationNeperspermeter3=medium(3).attenuationdBcmMHz/dBperNeper*100*f0/1e6;
k3=2*pi/lambda3-j*attenuationNeperspermeter3;%波数

R = 75e-3;%ROC曲率半径
a = 30e-3;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,0,medium(1).density,medium(1).soundspeed);
%划分网格点
dx = lambda1/6; %网格点的步长
dy = lambda1/6;
dz = lambda1/6;

xmin=-1.5*a;%观察点坐标的范围
xmax=-xmin;
ymax=1.5*a;
ymin=-ymax;

zmin1=7.75e-3;
zmax1=30e-3;
zmin2=30e-3+dz;
zmax2=35e-3;
zmin3=35e-3+dz;
zmax3=75e-3;

x=xmin:dx:xmax;%网格点的分布
y=ymin:dy:ymax;
z2=zmin2:dz:zmax2;
z1=zmin1:dz:zmax1;
z3=zmin3:dz:zmax3;

nx=length(x);%网格点的点数
ny=length(y);
nz1=length(z1);
nz2=length(z2);
nz3=length(z3);

%换能器离散化为点声源  注意离散化之后点声源的大小
dr=lambda1/6;%球面分割成很多个圆环，对应的半径是r，dr是半径增加的步长
r_back=0:dr:a-dr;%点声源前一段弧长对应的r
r_after=dr:dr:a;%点声源后一段弧长对应的r
r=r_after-dr/2;%第i个环带对应的中心点的r
Sm=lambda1/6;%中间环带离散化后对应的dS的弧长等于lambda/6
median=length(r)/2+1;%取中间环带的索引
ntheta=round(2*pi*r(median)/Sm);%中间层一个环带的划分点数,取整
dtheta=2*pi./ntheta;%根据取整重新调整每个点声源的对应弧度
theta_after=dtheta:dtheta:2*pi;%每个环带离散成多个点对应的弧度数组
theta_back=0:dtheta:(2*pi-dtheta);%每个环带离散成多个点对应的弧度数组
theta=theta_after-dtheta/2;%第i个环带离散成点声源，dS中间的点对应的弧度数组
X=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
Y=sin(theta)'*r;%点源的y坐标
Z0=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z=repmat(Z0,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。

%rayleigh积分计算水中xz平面的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz1  %观察网格点z方向的坐标
        rn=sqrt((X-x(ix)).^2+(Y-0).^2+(Z-z1(iz)).^2);%观察点到点源的距离
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=dS.*exp(-1i.*k1.*rn)./rn;
        B=sum(sum(A));%对上述求得的值累加 
        pr1(ix,iz)=1i*medium(1).density*u*medium(1).soundspeed*k1/(2*pi)*B; %乘以相关参数得到声压p
    end
end

% interface
z_interface1=30e-3;
z_interface2=35e-3;
%rayleigh积分计算经过界面1（介质1和介质2之间）处的法向阵速
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iy=1:ny  %观察网格点z方向的坐标      
        rm1=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z_interface1).^2);%interface1到点源的距离
        theta_m1i=acos((z_interface1-Z)./rm1);
        theta_m2t=asin(medium(2).soundspeed/medium(1).soundspeed.*sin(theta_m1i));
        Tm1=2*medium(1).soundspeed*medium(1).density.*cos(theta_m2t)./(medium(2).soundspeed*medium(2).density.*cos(theta_m1i)+medium(1).soundspeed*medium(1).density.*cos(theta_m2t));
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS_m1=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=exp(-j.*k1.*rm1)./rm1.*(1-j./(k1.*rm1)).*abs(Tm1).*(cos(theta_m1i)).*dS_m1;
        B=sum(sum(A));%对上述求得的值累加 
        u_m2(ix,iy)=j*u*k1/(2*pi)*B; 
    end
end
toc
% x1=(xmin+dx/2):dx:(xmax-dx/2)
X1=repmat(x',1,ny);
Y1=repmat(y,nx,1);
Z1=repmat(z_interface1,nx,ny);
%计算介质2中的声压分布
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz2  %观察网格点z方向的坐标      
        rm2=sqrt((X1-x(ix)).^2+(Y1-0).^2+(Z1-z2(iz)).^2);%interface1到点源的距离
        dS_m2=dx*dy; % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        Am2=u_m2.*exp(-j.*k2.*rm2)./rm2*dS_m2;
        Bm2=sum(sum(Am2));%对上述求得的值累加 
        pr2(ix,iz)=j*k2*medium(2).soundspeed*medium(2).density/(2*pi)*Bm2; %乘以相关参数得到声压p
    end
end

%rayleigh积分计算经过界面2（介质2和介质3之间）处的法向阵速
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iy=1:ny  %观察网格点z方向的坐标      
        rm2=sqrt((X1-x(ix)).^2+(Y1-0).^2+(Z1-z2(iz)).^2);%interface1到点源的距离
        theta_m2i=theta_m2t;
        theta_m3t=asin(medium(3).soundspeed/medium(2).soundspeed.*sin(theta_m2i));
        Tm2=2*medium(2).soundspeed*medium(2).density.*cos(theta_m3t)./(medium(3).soundspeed*medium(3).density.*cos(theta_m2i)+medium(2).soundspeed*medium(2).density.*cos(theta_m3t));
        dS_m2=dx*dy;
        A=exp(-j.*k2.*rm2)./rm2.*(1-j./(k2.*rm2)).*abs(Tm2).*(cos(theta_m2i)).*dS_m2;
        B=sum(sum(A));%对上述求得的值累加 
        u_m3(ix,iy)=j*u*k2/(2*pi)*B; 
    end
end
toc

Z2=repmat(z_interface2,nx,ny);
%计算介质2中的声压分布
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz3  %观察网格点z方向的坐标      
        rm3=sqrt((X1-x(ix)).^2+(Y1-0).^2+(Z2-z3(iz)).^2);%interface1到点源的距离
        dS_m3=dx*dy; % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        Am3=u_m3.*exp(-j.*k3.*rm3)./rm3*dS_m3;
        Bm3=sum(sum(Am3));%对上述求得的值累加 
        pr3(ix,iz)=j*k3*medium(3).soundspeed*medium(3).density/(2*pi)*Bm3; %乘以相关参数得到声压p
    end
end
pr3_abs=abs(pr3);
pr_max=max(pr3_abs(:));

x_median=floor(length(x)/2)+1;
focus_index_1=find(pr2_abs==pr_max);
s=size(pr2_abs);
[x_index,z_index]=ind2sub(s,focus_index_1);%将最大值单下标转为三维多下标
focus_index=z_index;
focus_forward=R-z(focus_index);
%轴向-6dB
axial_dB_index=find(pr2_abs(x_median,:)>=0.5*pr_max);
axial_dB=z2(max(axial_dB_index))-z2(min(axial_dB_index));

radial_dB_index=find(pr2_abs(:,z_index)>=0.5*pr_max);
radial_dB=x(max(radial_dB_index))-x(min(radial_dB_index));

%两部分矩阵拼接
pr=[pr1 pr2 pr3];
pr_abs=abs(pr);
z=[z1 z2 z3];

figure(1);
surf(z*1000,x*1000,pr_abs);
shading interp;
xlabel('z (mm)');
ylabel('x (mm)');
axis equal;
shading flat;
title('Rayleigh ');


