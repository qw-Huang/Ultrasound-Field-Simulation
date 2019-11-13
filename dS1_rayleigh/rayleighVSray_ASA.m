%使用自定义rayleigh积分求三维空间中的声强分布,对比自定义rayleigh积分和自定义rayleigh积分+ASA的感兴趣体误差  
clc;
clear all;
clear all;
f0=1e6;u=1;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质：单层->水，可改成多层，用函数set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数

R = 5 *1* 2 * lambda;%ROC曲率半径
a = 5 * lambda;%注意这里a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
%划分网格点
xmin=-0.7*a;%观察点坐标的范围
xmax=-xmin;
ymax=xmax;
ymin=-ymax;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %网格点的步长
dy = lambda/6; 
dz = lambda/6;

x=xmin:dx:xmax;
y=ymin:dy:ymax;
z=zmin:dz:zmax;

nx=length(x);
ny=length(y);
nz=length(z);

error_zeros=zeros(nx,ny,nz);

%换能器离散化为点声源  注意离散化之后点声源的大小
dr=lambda/6;%球面分割成很多个圆环，对应的半径是r，dr是半径增加的步长
r_back=0:dr:a-dr;%点声源前一段弧长对应的r
r_after=dr:dr:a;%点声源后一段弧长对应的r
r=r_after-dr/2;%第i个环带对应的中心点的r
Sm=lambda/6;%最外层环带（2pi*a）离散化后对应的dS的弧长小于lambda/6
ntheta=round(2*pi*(a-dr/2)/Sm);%最外层一个环带的划分点数,取整
dtheta=2*pi./ntheta;%根据取整重新调整每个点声源的对应弧度
theta_after=dtheta:dtheta:2*pi;%每个环带离散成多个点对应的弧度数组
theta_back=0:dtheta:(2*pi-dtheta);%每个环带离散成多个点对应的弧度数组
theta=theta_after-dtheta/2;%第i个环带离散成点声源，dS中间的点对应的弧度数组
X=cos(theta)'*r;%点源的x坐标  矩阵形式 ntheta*nr
Y=sin(theta)'*r;%点源的y坐标
Z_ring=R-sqrt(R*R-r.*r);%点声源三维空间中的z坐标
Z=repmat(Z_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。

%自定义rayleigh积分计算三维空间的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标划分
    for iy=1:ny %观察网格点y方向的坐标划分
       for iz=1:nz  %观察网格点z方向的坐标划分
            rn=sqrt((X-x(ix)).^2+(Y-y(iy)).^2+(Z-z(iz)).^2);%观察点到点源的距离
            dS_ring=r.*dtheta.*R.*(asin(r_after./R)-asin(r_back./R));%每一环带离散的dS的大小，第i环离散的点声源面积dS=第i环的离散弧长*dr对应的短弧
            dS=repmat(dS_ring,ntheta,1); % repmat( A , m , n )：将向量／矩阵在垂直方向复制m次，在水平方向复制n次。
            A=dS.*exp(-1i.*k.*rn)./rn;
            B=sum(sum(A));%对上述求得的值累加
            pr(ix,iy,iz)=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p   
       end 
    end
end
toc 

I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
I_asa=rayleigh_ASA(R,a,u,f0);%自定义rayleigh+ASA求声强
error_I=abs(I_pr-I_asa)./I_asa;%计算三维空间中两种方法的误差

%找到感兴趣的区域声强-6dB范围，画出这个范围内的error;找到最大值，返回对应的三维坐标，在最大点对应的xy平面，把感兴趣区域画出来
max_index=find_maxpoint(I_asa);%返回最大点位置坐标
max_index2=find_maxpoint(I_pr);
x_index=max_index(1);%最大值位置x下标
y_index=max_index(2);%最大值位置y下标
z_index=max_index(3);%最大值位置z下标

I_pr_max=max(I_pr(:));%自定义rayleigh 3D声强最大值
I_asa_max=max(I_asa(:));

%选择三维空间感兴趣区域误差
[index]=find(I_pr>=0.25*I_pr_max);
error_zeros(index)=error_I(index);

%画出感兴趣区域误差
figure(1);
surf(y*1000,x*1000,error_zeros(:,:,z_index));
shading interp;
axis equal;
xlabel('y (mm) ');
ylabel('x (mm) ');
title('rayleigh VS rayleigh+ASA 焦点处xy平面误差情况');

figure(2);
error_xz=squeeze(error_zeros(:,y_index,:));
surf(error_xz);
axis equal;
shading interp;
xlabel('z');
ylabel('x');
title('rayleigh VS rayleigh+ASA 焦点处xz平面误差情况');

figure(3);
I_pr_xz=finddB(I_pr(:,y_index,:),nx,nz);
pcolor(I_pr_xz);
axis equal;
shading interp; 
colorbar;
xlabel('z');
ylabel('x');
title('rayleigh 焦点处xz平面声强分布情况');

figure(4);
I_asa_xz=finddB(I_asa(:,y_index,:),nx,nz);
pcolor(I_asa_xz);
axis equal
shading interp; 
colorbar;
xlabel('z');
ylabel('x');
title('rayleigh+ASA 焦点处xz平面声强分布情况');

figure(5);
plot(max(error_xz));
xlabel('z');
ylabel('error');
title('rayleigh VS rayleigh+ASA 感兴趣体xz截面最大误差的变化情况');
% figure(7);
% plot(mean(error_xz));
% xlabel('z');
% ylabel('error');
% title('rayleigh VS rayleigh+ASA 感兴趣体xy截面平均误差的变化情况');

