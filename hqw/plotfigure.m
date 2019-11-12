clc;
clear all;
load('Rayleigh_3DF.mat');
figure(1);
    axis equal;
    surf(abs(squeeze(I_prs(:,:,49)))/max(max(squeeze(I_prs(:,:,49)))));%对prs归一化   
     shading interp
    colorbar
    title('Rayleigh Sommerfeld Result');
    xlabel('y');
    ylabel('x ');
    zlabel('normalized pressure');
     figure(2);
    axis equal;
    surf(abs(squeeze(I_prs(x_index,:,:)))/abs(I_prs_max));%对prs归一化   
    shading interp
    colorbar
    title('Rayleigh Sommerfeld Result');
    xlabel('z');
    ylabel('y');
    zlabel('normalized pressure');
    
    figure(3);
    plot(squeeze(I_prs(x_index,y_index,:)/abs(I_prs_max)));
     figure(4);
    plot(squeeze(I_prs(:,y_index,z_index)/abs(I_prs_max)));
    figure(5);
    plot(squeeze(I_prs(x_index,:,z_index)/abs(I_prs_max)));
    
    clear all;
    load('Rayleigh_FNM_xy.mat');
    figure(6);
    surf(p2);
    shading interp