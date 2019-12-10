%改进的rayleigh积分，计算在两层介质中的分布 水5mm+肌肉10mm-20mm
%ROC =15mm

clc;
clear all;
j=sqrt(-1);
f0=1e6;
P=100;%定义频率和法向阵速
% medium = set_medium('lossless');%定义介质（单层：水），可以改成多层set_layered_medium
medium = set_layered_medium([0,5e-3],[set_medium('water'),set_medium('muscle')]);
lambda1 = medium(1).soundspeed/f0;%波长=c/f
k1=2*pi/lambda1-j*medium(1).attenuationdBcmMHz;%波数
k2=2*pi/lambda1-j*medium(2).attenuationdBcmMHz;%波数

R = 15e-3;%ROC曲率半径
a = 7.5e-3;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,0,medium(1).density,medium(1).soundspeed);
%划分网格点
xmin=-a;%观察点坐标的范围
xmax=-xmin;
ymax=a;
ymin=-ymax;
zDiff=0.7*d;
zmin=5.1e-3;
zmax=R+zDiff;

dx = lambda1/6; %网格点的步长
dy = lambda1/6;
dz = lambda1/6;


x=xmin:dx:xmax;%网格点的分布
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);%网格点的点数
ny=length(y);
nz=length(z);


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

% interface
z_interface1=5e-3;

%rayleigh积分计算经过界面1（介质1和介质2之间）处的法向阵速
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iy=1:ny  %观察网格点z方向的坐标      
        rm1=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z_interface1).^2);%interface1到点源的距离
        theta_m1=acos((Z-z_interface1)./rm1);
        theta_m2=asin(medium(2).soundspeed/medium(1).soundspeed*sin(theta_m1));
        Tm1=2*medium(2).soundspeed*medium(2).density.*cos(theta_m1)./(medium(2).soundspeed*medium(2).density.*cos(theta_m1)+medium(1).soundspeed*medium(1).density.*cos(theta_m2));
        dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
        dS_m1=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        A=dS_m1.*exp(-j.*k1.*rm1)./rm1.*(1-j./(k1.*rm1)).*Tm1.*cos(theta_m1);
        B=sum(sum(A));%对上述求得的值累加 
        u_m2(ix,iy)=j*u*k1/(2*pi)*B; %乘以相关参数得到声压p
    end
end
toc

X1=repmat(x',1,ny);
Y1=repmat(y,nx,1);
Z1=repmat(z_interface1,nx,ny);
%计算介质2中的声压分布
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标      
        rm2=sqrt((X1-x(ix)).^2+(Y1-0).^2+(Z1-z(iz)).^2);%interface1到点源的距离
        dS_m2=dx*dy; % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        Am2=u_m2.*exp(-j.*k2.*rm2)./rm2*dS_m2;
        Bm2=sum(sum(Am2));%对上述求得的值累加 
        pr(ix,iz)=j*k1*medium(2).soundspeed*medium(2).density/(2*pi)*Bm2; %乘以相关参数得到声压p
    end
end
figure(1);
pr_nor=abs(pr)./abs(max(pr(:)));
surf(z*1000,x*1000,pr_nor);
shading interp;





