function Lag=Lagrangian(x,avt)
Lag=0;
j=1;
global n Y_n_a Y_n_r Y_line_net reply voltage M M1 P_d_k ...
    P_k_Max P_k_Min Q_d_k Q_k_Max Q_k_Min G_conn c2 c1 V_k_M V_k_m lcount
for i=1:n
    if G_conn(i)==1
        Lag=Lag+(c2(j))*((trace((Y_n_a{i})*x)+P_d_k(i))^2)+(c1(j))*(trace((Y_n_a{i})*x)+P_d_k(i))+...
            -(1/avt)*(log(P_k_Max(i)-P_d_k(i)-trace((Y_n_a{i})*x))+log(-P_k_Min(i)+P_d_k(i)+trace((Y_n_a{i})*x))+...
            log(Q_k_Max(i)-Q_d_k(i)-trace((Y_n_r{i})*x))+log(-Q_k_Min(i)+Q_d_k(i)+trace((Y_n_r{i})*x))+...
            log(V_k_M-trace((M{i})*x))+log(-V_k_m+trace((M{i})*x)));
        j=j+1;
    end
end
for c=1:lcount
    Lag=Lag-(1/avt)*(log(reply-trace((Y_line_net{c})*x))+log((voltage^2)-trace((M1{c})*x)));
end
Lag=Lag-(1/avt)*(log(det(x)));            
        