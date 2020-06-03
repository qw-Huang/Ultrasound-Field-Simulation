%使用内置FNM和ASA，计算多层组织 自变量fnumber
clc;
clear all;
close all;
%set up the tranducer

f0=1e6;
P=100;
R = 100e-3:5e-3:150e-3;
fnumber=0.8;

% Use a layered medium
medium1 = set_medium('water');
medium2 = set_medium('muscle');

% Center frequency and wavelength
lambda = medium1.soundspeed/f0;
z_interface=30e-3;
ndiv = 200;

tic
for i=1:length(R)
    % Create planar array of rectangular transducers
    a(i)=R(i)/(2*fnumber);
    u(i)=normal_velocity(P,R(i),a(i),0,medium1.density,medium1.soundspeed);
    d(i) = sqrt(R(i)^2 - a(i)^2);%理论焦点到孔径中心的距离
    z0=R(i)-d(i)+lambda;   % Determine where the source pressure will be calculated
    
    % Set up the coordinate grid
    xmin = -1.5*a(i);
    xmax = -xmin;
    ymin = -1.5*a(i);
    ymax = -ymin;
    zmin1=z0;
    zmax1= z_interface;
    zmin2 = z_interface;
    zmax2 = R(i)+20e-3;

    dx = lambda/6;
    dy = lambda/6;
    dz = lambda/6;

    x = xmin:dx:xmax;
    y = ymin:dy:ymax;
    z1=zmin1:dz:zmax1;
    z2 = zmin2:dz:zmax2;

    y_median = floor((ymax-ymin)/2/dy)+1;
    x_median = floor(length(x)/2)+1;

    cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
    cg_3d1 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin1,zmax1);
    cg_3d2 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin2,zmax2);
    
    xdcr_array = get_spherical_shell(a(i),R(i));% 球形换能器
    tic
    % Calculate the pressure
    p0 = fnm_call(xdcr_array,cg_p0,medium1,ndiv,f0,0);
    p0=p0.*u(i);%把法向阵速u=1转换为P=100
    [p_asa1,p_interface]=layer_cw_angular_spectrum(p0,cg_3d1,medium1,medium2,f0,1024,'Pa');
    p_asa2=cw_angular_spectrum(p_interface,cg_3d2,medium2,f0,1024,'Pa');
    toc
    p_asa_abs=abs(p_asa2);
    p_asa_max(i)=max(p_asa_abs(:));
    focus_index_max(i)=find(p_asa_abs==p_asa_max(i));
    s=size(p_asa_abs);
    [x_index,y_index,z_index]=ind2sub(s,focus_index_max(i));%将最大值单下标转为三维多下标
    focus_index_z(i)=z_index;
    focus_forward(i)=R(i)-z2(focus_index_z(i));
    %focal region -6dB
    axial_dB_index=find(p_asa_abs(x_median,y_median,:)>=0.5*p_asa_max(i));
    axial_dB(i)=z2(max(axial_dB_index))-z2(min(axial_dB_index));
    radial_dB_index=find(p_asa_abs(:,y_median,z_index)>=0.5*p_asa_max(i));
    radial_dB(i)=x(max(radial_dB_index))-x(min(radial_dB_index));
end
toc
figure(1);
plot(R,focus_forward,'k');
hold on;
figure(2);
plot(R,p_asa_max/10^6,'k');
hold on;
figure(3);
plot(R,axial_dB*1000,'k');
hold on;
figure(4);
plot(R,radial_dB*1000,'k');
hold on;
% 
% 
% % %介质全部是水

% tic
% tic
% for i=1:length(R)
%     z0=R(i)-d(i)+lambda;   % Determine where the source pressure will be calculated
%     
%     % Set up the coordinate grid
%     xmin = -1.5*a(i);
%     xmax = -xmin;
%     ymin = -1.5*a(i);
%     ymax = -ymin;
%     zmin1=z0;
%     zmax1= z_interface;
%     zmin2 = z_interface;
%     zmax2 = R(i)+20e-3;
% 
%     dx = lambda/6;
%     dy = lambda/6;
%     dz = lambda/6;
% 
%     x = xmin:dx:xmax;
%     y = ymin:dy:ymax;
%     z1 = zmin1:dz:zmax1;
%     z2 = zmin2:dz:zmax2;
% 
%     y_median = floor((ymax-ymin)/2/dy)+1;
%     x_median=floor(length(x)/2)+1;
% 
%     cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
%     cg_3d1 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin1,zmax1);
%     cg_3d2 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin2,zmax2);
%     
%     xdcr_array = get_spherical_shell(a(i),R(i));% 球形换能器
%     tic
%     % Calculate the pressure
%     p0 = fnm_call(xdcr_array,cg_p0,medium1,ndiv,f0,0);
%     p0=p0.*u(i);%把法向阵速u=1转换为P=100
%     [p_asa1,p_interface]=layer_cw_angular_spectrum(p0,cg_3d1,medium1,medium1,f0,1024,'Pa');
%     p_asa2_w=cw_angular_spectrum(p_interface,cg_3d2,medium1,f0,1024,'Pa');
%     toc
%     p_asa_abs_w=abs(p_asa2_w);
%     p_asa_max_w(i)=max(p_asa_abs_w(:));
%     focus_index_max_w(i)=find(p_asa_abs_w==p_asa_max_w(i));
%     s=size(p_asa_abs_w);
%     [x_index,y_index,z_index]=ind2sub(s,focus_index_max_w(i));%将最大值单下标转为三维多下标
%     focus_index_zw(i)=z_index;
%     focus_forward_w(i)=R(i)-z2(focus_index_zw(i));
%     %focal region -6dB
%     axial_dB_index=find(p_asa_abs_w(x_median,y_median,:)>=0.5*p_asa_max_w(i));
%     axial_dB_w(i)=z2(max(axial_dB_index))-z2(min(axial_dB_index));
%     radial_dB_index=find(p_asa_abs_w(:,y_median,z_index)>=0.5*p_asa_max_w(i));
%     radial_dB_w(i)=x(max(radial_dB_index))-x(min(radial_dB_index));
% end
% toc
% figure(1);
% plot(R,focus_forward_w*1000,'b');
% xlabel('ROC');
% ylabel('focus forward distance(mm)');
% figure(2);
% plot(R,p_asa_max_w/10^6,'b');
% xlabel('ROC');
% ylabel('pmax(MPa)');
% figure(3);
% plot(R,axial_dBw*1000,'b');
% xlabel('ROC');
% ylabel('axial -6dB distance(mm)');
% figure(4);
% plot(fnumber,radial_dBw*1000,'b');
% xlabel('ROC');
% ylabel('radial -6dB distance(mm)');
