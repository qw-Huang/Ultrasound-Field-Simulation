%调用range_dB function 得到对应平面的轴向和径向-6dB的大小，本程序为，改变f-number，得到不同的-6dB
%range的值，与理论值作比较，计算误差结果
clc;clear all;
f0=1e6;%    fnumber=R/(2*a);%所以f-number=曲率半径/孔径（2*a）
medium = set_medium('lossless');%定义介质：单层->水
lambda = medium.soundspeed/f0;%波长=c/f 1.5mm   lambda/6=0.25mm
fnumber=0.7:0.05:1.35;
for i=1:length(fnumber)
    a = 5 * 1.5e-3;%注意这里a是孔径的一半
    R = 2*a*fnumber(i);%ROC曲率半径
    range=range_dB_ASA(R,a,f0);
    range_radius(i)=range(1);
    range_axial(i)=range(2);
end
radius_ref=lambda.*fnumber.^2;%径向-6dB参考值
axial_ref=7*lambda.*fnumber.^2;%轴向-6dB
error_radius=abs(range_radius-radius_ref)./radius_ref;
error_axial=abs(range_axial-axial_ref)./axial_ref;

figure(1);
plot(fnumber,error_radius);
title('改变f-number，径向-6dB范围大小');
xlabel('f-number');
ylabel('radius error');

figure(2);
plot(fnumber,error_axial);
title('改变f-number，轴向-6dB范围大小');
xlabel('f-number');
ylabel('axial error');

