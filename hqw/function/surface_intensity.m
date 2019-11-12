function [output]=surface_intensity(P,R,a,hole_a,density,c)
%����Ϊ������P�������������R,a���ܶȣ����٣�����������Դ�������ǿ
H=R-sqrt(R^2 - a^2);
H_hole=R-sqrt(R^2 - hole_a^2);
S=2*pi*R*(H-H_hole);  %��������ı����
I_surface=P/S;
output=I_surface;
end

