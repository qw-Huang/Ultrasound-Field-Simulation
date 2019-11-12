function [output] = finddB(I,nx,ny)
%FIND-6DB �ҵ���άƽ���ڵ�-6dB�ķ�Χ��������-6dB��ǿ��ֵ��ŵ������ά������
%I��һ����ά����nx*ny������ľ���IΪ��������Ķ�άƽ������ǿ�ֲ���ֵ������һ��ȫ����󣬸���Ȥ�ĵ㱻��ֵ��ȫ�������
%   ����ѡ��õĶ�ά����
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

