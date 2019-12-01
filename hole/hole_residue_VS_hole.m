%使用dS1点声源划分方法，在2d平面xz面，加开圆孔的换能器声强计算，比较开孔部分、剩余部分和完整换能器的声压幅值、相位
clear all;
close all;
f0=1e6;%定义频率和声功率
P=100;
medium = set_medium('lossless');%定义介质（单层：水），可以改成多层set_layered_medium
lambda = medium.soundspeed/f0;%波长=c/f
k=2*pi/lambda;%波数
R = 1 * 5 * 2 * lambda;%ROC曲率半径
a = 1 * 5 * lambda;%注意这里的a是孔径的一半
fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离
u=normal_velocity(P,R,a,0,medium.density,medium.soundspeed);

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

x=xmin:dx:xmax;%网格点的分布
z=zmin:dz:zmax;

% 求y=0时xz平面
y=0;

%加开孔剩余部分
hole_a=15e-4;
pr_res=hole_rayleigh2(R,a,u,hole_a,x,z); %乘以相关参数得到声压p

%被挖掉的开孔
pr_hole=hole_rayleigh2(R,hole_a,u,0,x,z); %乘以相关参数得到声压p

%完整换能器
pr_all=hole_rayleigh2(R,a,u,0,x,z); %乘以相关参数得到声压p

% %选择感兴趣区域误差
% [index]=find(I_pr>=0.25*I_pr_max);
% contribution_res=I_pr_res./I_pr;
% contribution_hole=I_pr_hole./I_pr;
% res_zeros_xz(index)=contribution_res(index);
% hole_zeros_xz(index)=contribution_hole(index);

%比较声压的幅值和相位
pr_all_amplitude=abs(pr_all);
pr_all_phase=angle(pr_all);
pr_sum=pr_hole+pr_res;
pr_sum_amplitude=abs(pr_sum);
pr_sum_phase=angle(pr_sum);
error_phase=abs(pr_sum_phase-pr_all_phase)./pr_all_phase;
error_amplitude=abs(pr_sum_amplitude-pr_all_amplitude)./pr_all_amplitude;
I_all=abs(pr_all).^2./(2*medium.soundspeed*medium.density);
I_hole=abs(pr_hole).^2./(2*medium.soundspeed*medium.density);
I_res=abs(pr_res).^2./(2*medium.soundspeed*medium.density);
max_hole=find_maxpoint(I_hole);
max_res=find_maxpoint(I_res);
max_all=find_maxpoint(I_all);
% %画图
% % 对比剩余部分和去掉部分对声场分布的贡献
figure(1);
surf(z*1000, x*1000,pr_all_amplitude); 
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('pr_all_amplitude(R=15mm,a=7.5mm,20%a)');
figure(2);
surf(z*1000, x*1000,abs(pr_hole));
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('pr_hole_amplitude(R=15mm,a=7.5mm,20%a)');
figure(3);
surf(z*1000, x*1000,abs(pr_res));
axis equal;
shading interp;
colorbar;xlabel('z (mm) ');
ylabel('x (mm) ');
title('pr_res_amplitude(R=15mm,a=7.5mm,20%a)');
 figure(4);
surf(z*1000, x*1000,error_amplitude); 
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the error of amplitude(R=15mm,a=7.5mm,20%a)');
figure(5);
surf(z*1000, x*1000,error_phase);
axis equal;
shading interp;
colorbar;
xlabel('z (mm) ');
ylabel('x (mm) ');
title('the error of phase(R=15mm,a=7.5mm,20%a)');

% figure(6);
% plot(abs(pr_all(31,:)));
% hold on;
% plot(abs(pr_res(31,:)));
% hold on;
% plot(abs(pr_hole(31,:)));
% figure(7);
% plot(angle(pr_all(31,:)));
% hold on;
% plot(angle(pr_res(31,:)));
% hold on;
% plot(angle(pr_hole(31,:)));

% shading interp ;%去掉网格，平滑曲面
% axis equal;
% colorbar  %加颜色条
% xlabel('z (mm) ');
% ylabel('x (mm) ');
% title('(R=15mm,a=7.5mm,30%a)');
% save('R15_hole0.3a_holeVSresidue.mat','R','a','pr','pr_hole','pr_res','I_pr_hole','I_pr_res','I_pr','res_zeros_xz','hole_zeros_xz');
