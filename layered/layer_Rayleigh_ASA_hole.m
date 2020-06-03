%rayleigh积分+ASA  加开孔 pr_res=pr_all-pr_hole 
clc;
clear all;

% Set up the array
R=75e-3;
a=30e-3;
n=0:0.05:0.5;
% Use a layered medium
medium1 = set_medium('water');
medium2 = set_medium('muscle');
P=100;

% Center frequency and wavelength
f0 = 1e6;
lambda = medium1.soundspeed/f0;
d = sqrt(R^2 - a^2);%理论焦点到孔径中心的距离

% Set up the coordinate grid
xmin = -1.5*a;
xmax = -xmin;
ymin = -1.5*a;
ymax = -ymin;
zmin = 30e-3;
zmax = 90e-3;

dx = lambda/6;
dy = lambda/6;
dz = lambda/6;

x = xmin:dx:xmax;
y = ymin:dy:ymax;
z = zmin:dz:zmax;

% Determine where the source pressure will be calculated
z0=7.75e-3;
y_median = floor((ymax-ymin)/2/dy)+1;
x_median=floor(length(x)/2)+1;

zmin1=z0;
zmax1=30e-3;
z1=zmin1:dz:zmax1;

cg_3d1 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin1,zmax1);
cg_3d = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin,zmax);

for i=1:length(n)
    hole_a(i)=n(i)*a;
    u(i)=normal_velocity(P,R,a,hole_a(i),medium1.density,medium1.soundspeed);
    tic
    pr_all=rayleigh_2D_xy(R,a,f0,u(i),medium1,x,y,z0);
    pr_hole=rayleigh_2D_xy(R,hole_a(i),f0,u(i),medium1,x,y,z0);
    pr_res=pr_all-pr_hole;
    pr_res=single(pr_res);
    [p_asa1,p_interface]=layer_cw_angular_spectrum(pr_res,cg_3d1,medium1,f0,1024,'Pa');
    p_asa=cw_angular_spectrum(p_interface,cg_3d,medium2,f0,1024,'Pa');
   
    p_asa_abs=abs(p_asa);
    p_asa_max(i)=max(p_asa_abs(:));
    focus_index_1(i)=find(p_asa_abs==p_asa_max(i));
    s=size(p_asa_abs);
    [x_index,y_index,z_index]=ind2sub(s,focus_index_1(i));%将最大值单下标转为三维多下标
    focus_index(i)=z_index;
    focus_forward(i)=R-z(focus_index(i));
    %轴向-6dB
    axial_dB_index=find(p_asa_abs(x_median,y_median,:)>=0.5*p_asa_max(i));
    axial_dB(i)=z(max(axial_dB_index))-z(min(axial_dB_index));
    radial_dB_index=find(p_asa_abs(:,y_median,z_index)>=0.5*p_asa_max(i));
    radial_dB(i)=x(max(radial_dB_index))-x(min(radial_dB_index));

end
toc
figure(1);
scatter(n,focus_forward,25,'.','k');
hold on;
figure(2);
scatter(n,p_asa_max,25,'.','k');
hold on;
figure(3);
scatter(n,axial_dB,25,'.','k');
hold on;
figure(4);
scatter(n,radial_dB,25,'.','k');
hold on;

% %介质全是水
% for i=1:length(n)
%     hole_a(i)=n(i)*a;
%     u(i)=normal_velocity(P,R,a,hole_a(i),medium1.density,medium1.soundspeed);
%     tic
%     pr_all_w=rayleigh_2D_xy(R,a,f0,u(i),medium1,x,y,z0);
%     pr_hole_w=rayleigh_2D_xy(R,hole_a(i),f0,u(i),medium1,x,y,z0);
%     pr_res_w=pr_all_w-pr_hole_w;
%     [p_asa1,p_interface]=layer_cw_angular_spectrum(pr_res_w,cg_3d1,medium1,f0,1024,'Pa');
%     p_asa_w=cw_angular_spectrum(p_interface,cg_3d,medium1,f0,1024,'Pa');
%    
%     p_asa_abs_w=abs(p_asa_w);
%     p_asa_max_w(i)=max(p_asa_abs_w(:));
%     focus_index_1w(i)=find(p_asa_abs_w==p_asa_max_w(i));
%     s=size(p_asa_abs_w);
%     [x_index,y_index,z_index]=ind2sub(s,focus_index_1w(i));%将最大值单下标转为三维多下标
%     focus_index_w(i)=z_index;
%     focus_forward_w(i)=R-z(focus_index_w(i));
%     %轴向-6dB
%     axial_dB_index_w=find(p_asa_abs_w(x_median,y_median,:)>=0.5*p_asa_max_w(i));
%     axial_dB_w(i)=z(max(axial_dB_index_w))-z(min(axial_dB_index_w));
%     radial_dB_index_w=find(p_asa_abs_w(:,y_median,z_index)>=0.5*p_asa_max_w(i));
%     radial_dB_w(i)=x(max(radial_dB_index_w))-x(min(radial_dB_index_w));
% 
% end
% 
% figure(1);
% scatter(n,focus_forward_w,25,'.','b');
% xlabel('hole/a');
% ylabel('focus forward distance(m)');
% figure(2);
% scatter(n,p_res_max_w,25,'.','b');
% xlabel('hole/a');
% ylabel('pmax(Pa)');
% figure(3);
% scatter(n,axial_dBw,25,'.','b');
% xlabel('hole/a');
% ylabel('axial -6dB distance(m)');
% figure(4);
% scatter(n,radial_dBw,25,'.','b');
% xlabel('hole/a');
% ylabel('radial -6dB distance(m)');
