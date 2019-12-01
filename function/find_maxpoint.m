function [max_index]=find_maxpoint(matrix_I)
%找到矩阵matrix_I中最大点的位置坐标，即对应的x,y,z三维下标，以数组的形式输出这三个下标的值
I_max=max(matrix_I(:));%三维空间中的声强最大值
s=size(matrix_I);%取数组的大小，三维的长度
max_index_s=find(matrix_I>=I_max);%计算最大值位置的单下标
[x_index,y_index,z_index]=ind2sub(s,max_index_s);%将最大值单下标转为三维多下标
max_index=[x_index,y_index,z_index];
end
