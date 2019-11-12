clc;
clear all;

f0=1e6;u=1;%定义频率和法向阵速
medium = set_medium('lossless');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数

R = 6 * 5 * 2 * lambda;%ROC曲率半径
a = 6 * 5 * lambda;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
phi0=asin(a/R);

%划分网格点
xmin=-a;%观察点坐标的范围
xmax=-xmin;
ymax=0;
ymin=0;
zDiff=0.7*d;
zmin=R-zDiff;
zmax=R+zDiff;

dx = lambda/6; %网格点的步长 
dz = lambda/6;

x=xmin:dx:xmax;
z=zmin:dz:zmax;

nx=length(x);
nz=length(z);
error_zeros_xz=zeros(nx,nz);
% 求y=0时xz平面
y=0;

%坐标转化
r=sqrt(x.^2+y.^2);
o_theta=acos(x./r);


%换能器离散化为点声源  注意离散化之后点声源的大小
dR=lambda/6;%每个环带的弧度相等
dphi=dR/R;%xz平面上对应的等角度
nphi=round(phi0/dphi); %环带的个数，
dR=R*phi0/nphi;%nphi取整之后重新调整弧度
dphi=dR/R;%调整之后dphi的值
phi_after=dphi:dphi:phi0;
phi_back=0:dphi:(phi0-dphi);
phi=phi_after-dphi/2;%ri对应的xz平面的角度
ri0=R*sin(phi);%ri是第i环对应半径
Sm=lambda/6;%每个圆环对应的圆弧间隔频率固定
ntheta=round(2*pi*ri0/Sm);  %第i个圆环离散的弧度个数
dtheta=2*pi./ntheta; %根据取整个数，重新调整每一环的dtheta
Sm0=dtheta.*ri0; %调整每一环离散的弧度
Sm=repelem(Sm0,ntheta);
theta=dthetarepet(dtheta,ntheta);%变成一个数组
ri=repelem(ri0,ntheta);
Z0=R-R*cos(phi); %每个半径为ri的圆环对应的z轴坐标
Z=repelem(Z0,ntheta);

%rayleigh积分计算xz平面的声场  
tic
for ix=1:nx  %观察网格点x方向的坐标
    for iz=1:nz  %观察网格点z方向的坐标
       if x(ix)==0
           rn=sqrt(ri.^2+r(ix).^2+(Z-z(iz)).^2-2*ri.*r(ix));%观察点到点源的距离
       else
           rn=sqrt(ri.^2+r(ix).^2+(Z-z(iz)).^2-2*ri.*r(ix).*cos(theta-o_theta(ix)));%观察点到点源的距离
       end
        dS=dR.*Sm; % 求点源面积=圆环的采样间隔dR*第i个环的离散的弧度
        A=dS.*exp(-1i.*k.*rn)./rn;
        B=sum(sum(A));%对上述求得的值累加
        p=1i*medium.density*u*medium.soundspeed*k/(2*pi)*B; %乘以相关参数得到声压p
        pr(ix,iz)=squeeze(p);%得到网格点上对应的值，取模
    end
end
toc

%比较dS的累加
S_discrete=sum(dS);%离散之后的全部面积累加
S_theory=2*pi*R*(R-d);%理论计算球面面积，即解析面积
dSi=Sm0.*dR.*ntheta;%求每一环面积的累加，也就是把每一列（70个数）累加
ri_back=R*sin(phi_back);%ri是第i环短弧端对应半径
ri_after=R*sin(phi_after);%ri是第i环场=长弧端对应半径
di_back=sqrt(R^2-ri_back.^2);%从0开始作为第i层环的d
di_after=sqrt(R^2-ri_after.^2);%从dr开始作为第i层环的d
dSi_theory=2*pi*R*(R-di_after)-2*pi*R*(R-di_back);%第i环的环带面积
error_dSi=abs(dSi-dSi_theory)./dSi;

%调用focus内置rayleigh函数rayleigh_cw()
omega = 2 * pi * f0;
phi0 = asin(a/R);
xdcr = get_spherical_shell(a,R); %得到球面换能器

delta = [dx 0 dz];%网格点步长
ps = set_coordinate_grid(delta, xmin, xmax, ymin, ymax, zmin, zmax);%网格点的划分
ndiv = 100;%积分的点数
dflag = 0;%如果=1，结果有什么不同？

tic
prs=rayleigh_cw(xdcr,ps,medium,ndiv,f0);%FOCUS自带的rayleigh计算
toc

I_prs=acousticintensity(prs,medium.density,medium.soundspeed); 
I_prs_nor=squeeze(I_prs./max(I_prs(:)));
I_pr=acousticintensity(pr,medium.density,medium.soundspeed); 
I_pr_nor=squeeze(I_pr./max(I_pr(:)));
error=abs(I_pr_nor-I_prs_nor)./I_pr_nor;


%画图
figure(1);
surf(z*1000,x*1000,I_pr_nor); %横坐标为z，纵坐标为x，颜色代表p1的大小  ，为什么矩阵是这样的对应关系？
shading interp
title('Rayleigh ');
axis equal;%定义坐标的比例相同
colorbar
xlabel('z ');
ylabel('x ');
zlabel('normalized pressure');

figure(2);
surf(z*1000, x*1000, I_prs_nor);%对prs归一化   为什么对应的坐标是z在前面？
shading interp;
colorbar;
axis equal;
title('rayleigh(FOCUS)');
xlabel('z ');
ylabel('x ');
zlabel('normalized pressure');

%选择感兴趣区域误差
[index]=find(I_prs_nor>=0.25);
error_zeros_xz(index)=error(index);

%对比两种方法画出来的结果
figure(3);
surf(z*1000,x*1000,error_zeros_xz); 
shading interp %
axis equal;
colorbar  %加颜色条
xlabel('z ');
ylabel('x');
title('error(R=120mm,a=60mm)');

figure(4);
% plot(z*1000,max(error_zeros_xz));
% datestr(now)
 histogram(dS,10);


figure(5);
plot(error_dSi);