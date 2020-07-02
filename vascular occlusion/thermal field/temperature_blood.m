function [T2] = temperature_blood(m,n,q,s,T1,A,delta_r,delta_x,delta_t)
%TEMPERATURE_BLOOD �˴���ʾ�йش˺�����ժҪ��
%   ���Ѫ�����¶ȷֲ�

T2 =( m .* (T1 * A + q .* T1) + n .* (-Tp_2 + 16 * Tp_1 + 16 * Tm_1 - Tm_2)...
    -4*delta_t*V0*(1-(rb/r0)^2)*(1/(12*delta_r^2)*(-Tq_2 + 8 * Tq_1 - 8 * Tn_1 + Tn_2)+1/(12*delta_x^2)*(-Tp_2 + 8 * Tp_1 - 8 * Tm_1 + Tm_2))  ...
    +4*T1 -T0 + s .* abs(pressure_xr).^2)./3; %ȥ��Ѫ����ע������+o * Ta

end

