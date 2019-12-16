function [output] = layered_hole_rayleigh_2D(R,a,hole_a,P,medium,x,z)
%LAYERED_RAYLEIGH_2D 计算两层介质的声压分布（水+肌肉）
%   此处显示详细说明
j=sqrt(-1);
f0=1e6;
lambda1 = medium(1).soundspeed/f0;%波长=c/f
lambda2 = medium(2).soundspeed/f0;%波长=c/f
k1=2*pi/lambda1-j*medium(1).attenuationdBcmMHz;%波数
k2=2*pi/lambda2-j*medium(2).attenuationdBcmMHz;%波数

fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,hole_a,medium(1).density,medium(1).soundspeed);

dx1=lambda1/6;
dy1=lambda1/6;
x1=-a:dx1:a;
y1=x1;

nx1=length(x1);%网格点的点数
ny1=length(y1);
nx=length(x);
nz=length(z);

%换能器离散化为点声源  注意离散化之后点声源的大小
dr=lambda1/6;%球面分割成很多个圆环，对应的半径是r，dr是半径增加的步长
ndr=round((a-hole_a)/dr);
dr=(a-hole_a)/ndr;
r_back=(0+hole_a):dr:a-dr;%点声源前一段弧长对应的r
r_after=(dr+hole_a):dr:a;%点声源后一段弧长对应的r
r=r_after-dr/2;%第i个环带对应的中心点的r
Sm=lambda1/6;%中间环带离散化后对应的dS的弧长等于lambda/6
median=round(length(r)/2);%取中间环带的索引
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
for ix=1:nx1  %观察网格点x方向的坐标
    for iy=1:ny1  %观察网格点z方向的坐标      
        rm1=sqrt((X-x1(ix)).^2+(Y-y1(iy)).^2+(Z-z_interface1).^2);%interface1到点源的距离
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

X1=repmat(x1',1,ny1);
Y1=repmat(y1,nx1,1);
Z1=repmat(z_interface1,nx1,ny1);
%计算介质2中的声压分布
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标      
        rm2=sqrt((X1-x(ix)).^2+(Y1-0).^2+(Z1-z(iz)).^2);%interface1到点源的距离
        dS_m2=dx1*dy1; % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
        Am2=u_m2.*exp(-j.*k2.*rm2)./rm2*dS_m2;
        Bm2=sum(sum(Am2));%对上述求得的值累加 
        pr(ix,iz)=j*k2*medium(2).soundspeed*medium(2).density/(2*pi)*Bm2; %乘以相关参数得到声压p
    end
end

output = pr;

end

