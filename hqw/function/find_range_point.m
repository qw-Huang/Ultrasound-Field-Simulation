function [output] = find_range_point(matrix,nx,ny)
%FIND_RANGE_VALUE 输入是一个矩阵（这个矩阵是处理过的，只取声强大于-6dB的数值，存放到相同位置，小于则保持为0），找到每一行非零点对应的最小的列数和非零点对应的最大的列数,相减得到间隔的点数
%   此处显示详细说明
c1=0;c2=0;
for ix=1:nx
    for iy=1:ny
        if  matrix(ix,iy)>0&&matrix(ix,iy-1)==0 %如果出现零和非零的连续变化
            c1=iy;%非零点的列标被放到c1中
        else c1=c1;
        end
        if matrix(ix,iy)>0&& matrix(ix,iy+1)==0%如果在边界出现非零和零的连续变化，则认为是非零元素的另一个边界
            c2=iy;%非零点的列标被放到c1中，认为是列标的最大值
        else c2=c2;    
        end      
    end
    rangepoint(ix)=c2-c1;%列标之差，得到每一行非零数间隔的点数
end
output=max(rangepoint);
end

