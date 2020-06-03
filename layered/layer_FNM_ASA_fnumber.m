%使用内置FNM和ASA，计算多层组织 自变量fnumber
clc;
clear all;
close all;
%set up the tranducer

f0=1e6;
P= 100;
R = 140e-3;
fnumber=0.7:0.1:1.1;
depth=6e-2;
hole_a=2e-2;
% Use a layered medium
medium1 = set_medium('water');
medium2 = set_medium('water');

% Center frequency and wavelength
lambda = medium1.soundspeed/f0;
z_interface=R-depth;
ndiv = 200;

tic
for i=1:length(fnumber)
    % Create planar array of rectangular transducers
    a(i)=R/(2*fnumber(i));
    u(i)=normal_velocity(P,R,a(i),hole_a,medium1.density,medium1.soundspeed);
    d(i) = sqrt(R^2 - a(i)^2);%理论焦点到孔径中心的距离
    z0=R-d(i)+lambda;   % Determine where the source pressure will be calculated
    
    % Set up the coordinate grid
    xmin = -1.5*a(i);
    xmax = -xmin;
    ymin = -1.5*a(i);
    ymax = -ymin;
    zmin1=z0;
    zmax1= z_interface;
    zmin2 = z_interface;
    zmax2 = R+20e-3;

    dx = lambda/6;
    dy = lambda/6;
    dz = lambda/6;

    x = xmin:dx:xmax;
    y = ymin:dy:ymax;
    z1=zmin1:dz:zmax1;
    z2 = zmin2:dz:zmax2;

    y_median = floor((ymax-ymin)/2/dy)+1;
    x_median=floor(length(x)/2)+1;

    cg_p0 = set_coordinate_grid([dx dy 1], xmin,xmax,ymin,ymax,z0,z0);
    cg_3d1 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin1,zmax1);
    cg_3d2 = set_coordinate_grid([dx dy dz],xmin,xmax,ymin,ymax,zmin2,zmax2);
    
    xdcr_array = get_spherical_shell(a(i),R);% 球形换能器
    xdcr_array_hole = get_spherical_shell(hole_a,R);
    tic
    % Calculate the pressure
    p0_all = fnm_call(xdcr_array,cg_p0,medium1,ndiv,f0,0);
    p0_hole = fnm_call(xdcr_array_hole,cg_p0,medium1,ndiv,f0,0);
    p0_res=p0_all-p0_hole;
    p0=p0_res.*u(i);%把法向阵速u=1转换为P=100
    [p_asa1,p_interface]=layer_cw_angular_spectrum(p0,cg_3d1,medium1,medium2,f0,1024,'Pa');
    p_asa2=cw_angular_spectrum(p_interface,cg_3d2,medium2,f0,1024,'Pa');
    toc
    p_asa_abs=abs(p_asa2);
    p_asa_max(i)=max(p_asa_abs(:));
    focus_index_max(i)=find(p_asa_abs==p_asa_max(i));
    s=size(p_asa_abs);
    [x_index,y_index,z_index]=ind2sub(s,focus_index_max(i));%将最大值单下标转为三维多下标
    focus_index_z(i)=z_index;
    focus_forward(i)=R-z2(focus_index_z(i));
    %focal region -6dB
    axial_dB_index=find(p_asa_abs(x_median,y_median,:)>=0.5*p_asa_max(i));
    axial_dB(i)=z2(max(axial_dB_index))-z2(min(axial_dB_index));
    radial_dB_index=find(p_asa_abs(:,y_median,z_index)>=0.5*p_asa_max(i));
    radial_dB(i)=x(max(radial_dB_index))-x(min(radial_dB_index));
end
toc
