function [output]=normal_velocity(P,R,a,hole_a,density,c)
%����Ϊ������P�������������R,a���ܶȣ����٣�����������Դ����ķ�������
H=R-sqrt(R^2 - a^2);
H_hole=R-sqrt(R^2 - hole_a^2);
S=2*pi*R*(H-H_hole);  %��������ı����
u=sqrt(P/(S*density*c));%���ݹ�ʽ���㷨������
%I_surface=P/S;
output=u;
end



