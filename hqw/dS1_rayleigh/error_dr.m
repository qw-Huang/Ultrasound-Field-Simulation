clc;
clear all;
f0=1e6;u=1;%����Ƶ�ʺͷ�������
medium = set_medium('lossless');%������ʣ����㣺ˮ�������Ըĳɶ��set_layered_medium
lambda = medium.soundspeed/f0;%����=c/f
r=[lambda/60,lambda/66,lambda/72,lambda/78,lambda/84,lambda/90,lambda/96,lambda/102,lambda/108,lambda/114,lambda/120,lambda/126];
for i=2:12
    p=rayleigh_2D_xz(r(i));
    p0=rayleigh_2D_xz(r(i-1));
    error=abs((p-p0)./p0);
    error_max(i)=max(error(:));
end
plot(error_max);
