function Grad=Gradient(x,avt)
global n Y_n_a Y_n_r Y_line_net reply voltage M M1 P_d_k ...
    P_k_Max P_k_Min Q_d_k Q_k_Max Q_k_Min G_conn c2 c1 V_k_M V_k_m lcount
Grad=zeros(2*n,2*n);
j=1;
for e=1:n
    if G_conn(e)==1
        Grad=sparse(Grad)+2*c2(j)*sparse((trace((Y_n_a{e})*x)+P_d_k(e))*(Y_n_a{e}))+c1(j)*(Y_n_a{e});
        j=j+1;
    end
end
for i=1:n
    Grad=sparse(Grad)-(1/avt)*(sparse((Y_n_a{i})*(1/(trace((Y_n_a{i})*x)-P_k_Max(i)+P_d_k(i))))+...
        sparse((Y_n_a{i})*(1/(trace((Y_n_a{i})*x)-P_k_Min(i)+P_d_k(i))))+...
        sparse((Y_n_r{i})*(1/(trace((Y_n_r{i})*x)-Q_k_Max(i)+Q_d_k(i))))+...
        sparse((Y_n_r{i})*(1/(trace((Y_n_r{i})*x)-Q_k_Min(i)+Q_d_k(i))))+...
        sparse((M{i})*(1/(trace((M{i})*x)-V_k_M(i))))+...
        sparse((M{i})*(1/(trace((M{i})*x)-V_k_m(i)))));
end
for c=1:lcount
    Grad=sparse(Grad)-(1/avt)*(sparse((Y_line_net{c})*(1/(trace((Y_line_net{c})*x)-reply)))+...
        sparse((M1{c})*(1/(trace((M1{c})*x)-(voltage^2)))));
end
Grad=sparse(Grad)-(1/avt)*(sparse(inv(x)));