function [output] = finddB(I,nx,ny)
%FIND-6DB 找到二维平面内的-6dB的范围，将大于-6dB声强的值存放到输出二维矩阵中
%I是一个两维数组nx*ny，输入的矩阵I为计算出来的二维平面内声强分布的值，创建一个全零矩阵，感兴趣的点被赋值到全零矩阵中
%   返回选择好的二维数组
I_max=max(I(:));
output=zeros(nx,ny);
for ix=1:nx
    for iy=1:ny
        if I(ix,iy)>=0.25*I_max
           output(ix,iy)=I(ix,iy);
        else
           output(ix,iy)=0;
        end
    end
end
end

