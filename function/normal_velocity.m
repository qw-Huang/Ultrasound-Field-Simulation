function [output]=normal_velocity(P,R,a,hole_a,density,c)
%输入为声功率P，换能器表面的R,a，密度，声速，输出结果是声源表面的法向阵速
H=R-sqrt(R^2 - a^2);
H_hole=R-sqrt(R^2 - hole_a^2);
S=2*pi*R*(H-H_hole);  %计算球面的表面积
u=sqrt(P/(S*density*c));%根据公式计算法向阵速
%I_surface=P/S;
output=u;
end



