%% Authored by Nick LEE 2015 Dec
function [T,T_focus,D,Q_perfusion] = ThermalDose(A,B,D,T,Q,Q_perfusion,Qc,delta_t)
T = bicg(A,B*T+Q-Q_perfusion);
    T_focus = max(T);
    Q_perfusion = Qc * T;
    for i=1:length(Q)
        if T(i)<37
            T(i)=37;
        end
        if T(i)>43
            D(i)=D(i)+0.5^(43-T(i))*delta_t;
        else
            D(i)=D(i)+0.25^(43-T(i))*delta_t;
        end
        
    end
    