function [output] = find_range_point(matrix,nx,ny)
%FIND_RANGE_VALUE ������һ��������������Ǵ�����ģ�ֻȡ��ǿ����-6dB����ֵ����ŵ���ͬλ�ã�С���򱣳�Ϊ0�����ҵ�ÿһ�з�����Ӧ����С�������ͷ�����Ӧ����������,����õ�����ĵ���
%   �˴���ʾ��ϸ˵��
c1=0;c2=0;
for ix=1:nx
    for iy=1:ny
        if  matrix(ix,iy)>0&&matrix(ix,iy-1)==0 %���������ͷ���������仯
            c1=iy;%�������б걻�ŵ�c1��
        else c1=c1;
        end
        if matrix(ix,iy)>0&& matrix(ix,iy+1)==0%����ڱ߽���ַ������������仯������Ϊ�Ƿ���Ԫ�ص���һ���߽�
            c2=iy;%�������б걻�ŵ�c1�У���Ϊ���б�����ֵ
        else c2=c2;    
        end      
    end
    rangepoint(ix)=c2-c1;%�б�֮��õ�ÿһ�з���������ĵ���
end
output=max(rangepoint);
end

