%改变f-number（固定a，改变R），看焦点声强和表面声强的关系
clc
clear all;
f0=1e6;%    fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
medium = set_medium('lossless');%定义介质：单层->水
lambda = medium.soundspeed/f0;%波长=c/f 1.5mm   lambda/6=0.25mm
fnumber=0.7:0.05:1.4;
P=10;
for i=1:length(fnumber)
    a = 5 *10* 1.5e-3;%注意这里a是孔径的一半
    R = 2*a*fnumber(i);%ROC曲率半径
    u=normal_velocity(P,R,a,medium.density,medium.soundspeed);
    I_surface(i)=surface_intensity(P,R,a,medium.density,medium.soundspeed);
    I=rayleigh_ASA(R,a,u,f0);
    I_focus(i)=max(I(:));
end

figure(2);
plot(fnumber,I_surface);
title('改变f-number,换能器表面声强变化情况');
xlabel('f-number');
ylabel('surface intensity');

figure(3);
plot(I_surface,I_focus);
title('改变f-number，焦点声强随表面声强的变化情况');
xlabel('surface intensity');
ylabel('focus intensity');


