function Hess=Hessian(x,avt)
global n Y_n_a Y_n_r Y_line_net reply voltage M M1 P_d_k ...
    P_k_Max P_k_Min Q_d_k Q_k_Max Q_k_Min G_conn c2 V_k_M V_k_m lcount
Hess=zeros(((2*n)^2),((2*n)^2));
j=1;
for e=1:n
    if G_conn(e)==1
        Hess=sparse(Hess)+2*c2(j)*sparse(kron(Y_n_a{e},Y_n_a{e}));
        j=j+1;
    end
end
for i=1:n
    Hess=sparse(Hess)+(1/avt)*(sparse((kron(Y_n_a{i},Y_n_a{i}))*(1/((trace((Y_n_a{i})*x)-P_k_Max(i)+P_d_k(i))^2)))+...
        sparse((kron(Y_n_a{i},Y_n_a{i}))*(1/((trace((Y_n_a{i})*x)-P_k_Min(i)+P_d_k(i))^2)))+...
        sparse((kron(Y_n_r{i},Y_n_r{i}))*(1/((trace((Y_n_r{i})*x)-Q_k_Max(i)+Q_d_k(i))^2)))+...
        sparse((kron(Y_n_r{i},Y_n_r{i}))*(1/((trace((Y_n_r{i})*x)-Q_k_Min(i)+Q_d_k(i))^2)))+...
        sparse((kron(M{i},M{i}))*(1/((trace((M{i})*x)-V_k_M(i))^2)))+...
        sparse((kron(M{i},M{i}))*(1/((trace((M{i})*x)-V_k_m(i))^2))));
end
for c=1:lcount
    Hess=sparse(Hess)+(1/avt)*(sparse((kron(Y_line_net{c},Y_line_net{c}))*(1/((trace((Y_line_net{c})*x)-reply)^2)))+...
        sparse((kron(M1{c},M1{c}))*(1/((trace((M1{c})*x)-(voltage^2))^2))));
end
B=speye(2*n,2*n);
C=cell(2*n,2*n);
for i=1:2*n
    for j=1:2*n
        C{i,j}=sparse((x\(B(:,i)))*((B(j,:))/x));
    end
end
D=sparse(cell2mat(C));
Hess=sparse(Hess+D);