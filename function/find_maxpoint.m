function [max_index]=find_maxpoint(matrix_I)
%�ҵ�����matrix_I�������λ�����꣬����Ӧ��x,y,z��ά�±꣬���������ʽ����������±��ֵ
I_max=max(matrix_I(:));%��ά�ռ��е���ǿ���ֵ
s=size(matrix_I);%ȡ����Ĵ�С����ά�ĳ���
max_index_s=find(matrix_I>=I_max);%�������ֵλ�õĵ��±�
[x_index,y_index,z_index]=ind2sub(s,max_index_s);%�����ֵ���±�תΪ��ά���±�
max_index=[x_index,y_index,z_index];
end
