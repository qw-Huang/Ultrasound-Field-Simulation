%改变f-number（固定a，改变R），看-6dB范围的绝对大小clc;clear all
clc
clear all;
f0=1e6;%    fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
medium = set_medium('lossless');%定义介质：单层->水
lambda = medium.soundspeed/f0;%波长=c/f 1.5mm   lambda/6=0.25mm
fnumber=0.7:0.05:1.4;
for i=1:length(fnumber)
    a = 5 *1* 1.5e-3;%注意这里a是孔径的一半
    R = 2*a*fnumber(i);%ROC曲率半径
    range=range_dB_ASA(R,a,f0);
    range_radius(i)=range(1);
    range_axial(i)=range(2);
    range_z(i)=range(3);
end
figure(1);
plot(fnumber,range_radius);
title('径向-6dB范围大小');
xlabel('f-number');
ylabel('radius');

figure(2);
plot(fnumber,range_axial);
title('轴向-6dB范围大小');
xlabel('f-number');
ylabel('axial');

figure(3);
plot(fnumber,range_z);
title('改变f-number，焦点前移距离变化情况');
xlabel('f-number');
ylabel('z');
