function Tol=IP_Stop(x,avt)
Tol=0;
global n Y_n_a Y_n_r Y_line_net reply voltage M M1 P_d_k ...
    P_k_Max P_k_Min Q_d_k Q_k_Max Q_k_Min V_k_M V_k_m lcount
for i=1:n
        Tol=Tol+(1/avt)*(log(P_k_Max(i)-P_d_k(i)-trace((Y_n_a{i})*x))+log(-P_k_Min(i)+P_d_k(i)+trace((Y_n_a{i})*x))+...
            log(Q_k_Max(i)-Q_d_k(i)-trace((Y_n_r{i})*x))+log(-Q_k_Min(i)+Q_d_k(i)+trace((Y_n_r{i})*x))+...
            log(V_k_M-trace((M{i})*x))+log(-V_k_m+trace((M{i})*x)));
end
for c=1:lcount
    Tol=Tol+(1/avt)*(log(reply-trace((Y_line_net{c})*x))+log((voltage^2)-trace((M1{c})*x)));
end
Tol=Tol+(1/avt)*(log(det(x)));