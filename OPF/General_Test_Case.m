clc
clear;
n=input('Please Enter the particular IEEE Test Case You wish to run this Simulation for\nby indicating the number of buses.\nWarning:Valid Test Cases consist of 30, 57, 118 & 300 buses.\n');
if n==30
    Y_bus_r=xlsread('30bus_Ybus.xlsx','Sheet2','A1:AD30');
    Y_bus_i=xlsread('30bus_Ybus.xlsx','Sheet3','A1:AD30');
    Yline_G=xlsread('30bus_G.xlsx','Sheet1','B2:AE31');
    Yline_B=xlsread('30bus_B.xlsx','Sheet1','B2:AE31');
    Yshunt_B=xlsread('30bus_Bshunt.xlsx','Sheet1','B2:AE31');
    Yconn=xlsread('30bus_Connectivity.xlsx','Sheet1','B2:AE31');
    P_k_Max=xlsread('Limits.xlsx','30_Bus','B2:B31');
    P_k_Min=xlsread('Limits.xlsx','30_Bus','C2:C31');
    Q_k_Max=xlsread('Limits.xlsx','30_Bus','D2:D31');
    Q_k_Min=xlsread('Limits.xlsx','30_Bus','E2:E31');
    P_d_k=xlsread('30bus_load.xlsx','Sheet1','D2:D31');
    Q_d_k=xlsread('30bus_load.xlsx','Sheet1','E2:E31');
    G_conn=xlsread('Limits.xlsx','30_Bus','F2:F31');
    c1=xlsread('Limits.xlsx','30','I2:I7');
    c2=xlsread('Limits.xlsx','30','H2:H7');
    c0=xlsread('Limits.xlsx','30','G2:G7');
    g=6;
elseif n==57
        Y_bus_r=xlsread('57bus_Ybus.xlsx','Sheet2','A1:BE57');
        Y_bus_i=xlsread('57bus_Ybus.xlsx','Sheet3','A1:BE57');
        Yline_G=xlsread('57bus_G.xlsx','Sheet1','B2:BF58');
        Yline_B=xlsread('57bus_B.xlsx','Sheet1','B2:BF58');
        Yshunt_B=xlsread('57bus_Bshunt.xlsx','Sheet1','A1:BE57');
        Yconn=xlsread('57bus_Connectivity.xlsx','Sheet1','B2:BF58');
        P_k_Max=xlsread('Limits.xlsx','57_Bus','B2:B58');
        P_k_Min=xlsread('Limits.xlsx','57_Bus','C2:C58');
        Q_k_Max=xlsread('Limits.xlsx','57_Bus','D2:D58');
        Q_k_Min=xlsread('Limits.xlsx','57_Bus','E2:E58');
        P_d_k=xlsread('57bus_load.xlsx','Load Records','L3:L59');
        Q_d_k=xlsread('57bus_load.xlsx','Load Records','M3:M59');
        G_conn=xlsread('Limits.xlsx','57_Bus','F2:F58');
        c1=xlsread('Limits.xlsx','57','I2:I8');
        c2=xlsread('Limits.xlsx','57','H2:H8');
        c0=xlsread('Limits.xlsx','57','G2:G8');
        g=7;
elseif n==118
            Y_bus_r=xlsread('118bus_Ybus.xlsx','Sheet2','A1:DN118');
            Y_bus_i=xlsread('118bus_Ybus.xlsx','Sheet3','A1:DN118');
            Yline_G=xlsread('118bus_G.xlsx','Sheet1','B2:DO119');
            Yline_B=xlsread('118bus_B.xlsx','Sheet1','B2:DO119');
            Yshunt_B=xlsread('118bus_Bshunt.xlsx','Sheet1','B2:DO119');
            Yconn=xlsread('118bus_Connectivity.xlsx','Sheet1','B2:DO119');
            P_k_Max=xlsread('Limits.xlsx','118_Bus','B2:B119');
            P_k_Min=xlsread('Limits.xlsx','118_Bus','C2:C119');
            Q_k_Max=xlsread('Limits.xlsx','118_Bus','D2:D119');
            Q_k_Min=xlsread('Limits.xlsx','118_Bus','E2:E119');
            P_d_k=xlsread('118bus_load.xlsx','Load Records','L3:L120');
            Q_d_k=xlsread('118bus_load.xlsx','Load Records','M3:M120');
            G_conn=xlsread('Limits.xlsx','118_Bus','F2:F119');
            c1=xlsread('Limits.xlsx','118','I2:I55');
            c2=xlsread('Limits.xlsx','118','H2:H55');
            c0=xlsread('Limits.xlsx','118','G2:G55');
            g=54;
elseif n==300
                Y_bus_r=xlsread('300bus_Ybus.xlsx','Sheet2','A1:KN300');
                Y_bus_i=xlsread('300bus_Ybus.xlsx','Sheet3','A1:KN300');
                Yline_G=xlsread('300bus_G.xlsx','Sheet1','B2:KO301');
                Yline_B=xlsread('300bus_B.xlsx','Sheet1','B2:KO301');
                Yshunt_B=xlsread('300bus_Bshunt.xlsx','Sheet1','B2:KO301');
                Yconn=xlsread('300bus_Connectivity.xlsx','Sheet1','B2:KO301');
                P_k_Max=xlsread('Limits.xlsx','300_Bus','B2:B301');
                P_k_Min=xlsread('Limits.xlsx','300_Bus','C2:C301');
                Q_k_Max=xlsread('Limits.xlsx','300_Bus','D2:D301');
                Q_k_Min=xlsread('Limits.xlsx','300_Bus','E2:E301');
                P_d_k=xlsread('300bus_load.xlsx','Load Records','L3:L302');
                Q_d_k=xlsread('300bus_load.xlsx','Load Records','M3:M302');
                G_conn=xlsread('Limits.xlsx','300_Bus','F2:F301');
                c1=xlsread('Limits.xlsx','300','I2:I70');
                c2=xlsread('Limits.xlsx','300','H2:H70');
                c0=xlsread('Limits.xlsx','300','G2:G70');
                g=69;
else
                fprintf('Test case with this particular number of buses is not supported by the Program');
end
V_k_Max=1.06*ones(n,1);
V_k_Min=0.94*ones(n,1);
V_k_M=V_k_Max.*V_k_Max;
V_k_m=V_k_Min.*V_k_Min;
B=eye(n);
for j=1:n
    Y_r(:,:,j)=B(:,j)*B(j,:)*Y_bus_r;
    Y_i(:,:,j)=B(:,j)*B(j,:)*Y_bus_i;
    Y_n_a(:,:,j)=0.5*[Y_r(:,:,j)+transpose(Y_r(:,:,j)) transpose(Y_i(:,:,j))-Y_i(:,:,j);
        Y_i(:,:,j)-transpose(Y_i(:,:,j)) Y_r(:,:,j)+transpose(Y_r(:,:,j))];
    Y_n_r(:,:,j)=-0.5*[Y_i(:,:,j)+transpose(Y_i(:,:,j)) Y_r(:,:,j)-transpose(Y_r(:,:,j));
        transpose(Y_r(:,:,j))-Y_r(:,:,j) Y_i(:,:,j)+transpose(Y_i(:,:,j))];
    M(:,:,j)=[B(:,j)*B(j,:) zeros(n,n);zeros(n,n) B(:,j)*B(j,:)];
end
c=0;
for l=1:n
    for m=1:n
        if Yconn(l,m)==1
            c=c+1;
            Y_line_r(:,:,c)=Yline_G(l,m)*B(:,l)*B(l,:)-Yline_G(l,m)*B(:,l)*B(m,:);
            Y_line_i(:,:,c)=(Yline_B(l,m)+Yshunt_B(l,m))*B(:,l)*B(l,:)-Yline_B(l,m)*B(:,l)*B(m,:);
            Y_line_net(:,:,c)=0.5*[Y_line_r(:,:,c)+transpose(Y_line_r(:,:,c)) transpose(Y_line_i(:,:,c))-Y_line_i(:,:,c);Y_line_i(:,:,c)-transpose(Y_line_i(:,:,c)) Y_line_r(:,:,c)+transpose(Y_line_r(:,:,c))];
        end
    end
end
reply=input('Enter the pu value of the real power line limit\n');
cvx_begin sdp
variables lambdak_m(n) lambdak_M(n)  lambda_k_m(n) lambda_k_M(n) mu_k_M(n) mu_k_m(n) lambda_lm(c);
variables r_1(g) r_2(g);
dual variable A;
expression S(2*n,2*n,n+1);
S(:,:,1)=zeros(2*n,2*n);
for count=1:n
    S(:,:,count+1)=S(:,:,count)+(lambdak_M(count)-lambdak_m(count))*Y_n_a(:,:,count)+(lambda_k_M(count)-lambda_k_m(count))*Y_n_r(:,:,count)+(mu_k_M(count)-mu_k_m(count))*M(:,:,count);
end
expression s(2*n,2*n,c+1);
s(:,:,1)=zeros(2*n,2*n);
for lcount=1:c
    s(:,:,lcount+1)=s(:,:,lcount)+(lambda_lm(lcount))*Y_line_net(:,:,lcount);
end
expression R(2*g,2*g,g+1);
expression T(2*g,2*g,g);
R(:,:,1)=zeros(2*g,2*g);
for i=1:g
    for j=1:2*g
        for k=1:2*g
            if (j==2*i-1)&&(k==2*i-1)
                T(j,k,i)=1;
            elseif (j==2*i-1)&&(k==2*i)
                T(j,k,i)=r_1(i);
            elseif (j==2*i)&&(k==2*i-1)
                T(j,k,i)=r_1(i);
            elseif (j==2*i)&&(k==2*i)
                T(j,k,i)=r_2(i);
            else
                T(j,k,i)=0;
            end
        end
    end
end          
for gcount=1:g
  R(:,:,gcount+1)=R(:,:,gcount)+T(:,:,gcount);
end
expression N(2*n,2*n,g+1);
N(:,:,1)=zeros(2*n,2*n);
j=1;
for i=1:n
    if G_conn(i)==1
        N(:,:,j+1)=N(:,:,j)+(c1(j)+2*sqrt(c2(j))*r_1(j))*Y_n_a(:,:,i);
        j=j+1;
    end
end
expression O(g+1);
O(1)=0;
j=1;
for i=1:n
    if G_conn(i)==1
        O(j+1)=O(j)+(c1(j)+2*sqrt(c2(j))*r_1(j))*P_d_k(i);
        j=j+1;
    end
end
maximize transpose(lambdak_m)*P_k_Min-transpose(lambdak_M)*P_k_Max+transpose(lambdak_M-lambdak_m)*P_d_k+transpose(lambda_k_m)*Q_k_Min-transpose(lambda_k_M)*Q_k_Max+transpose(lambda_k_M-lambda_k_m)*Q_d_k+transpose(mu_k_m)*(V_k_m)-transpose(mu_k_M)*(V_k_M)+(transpose(c0-r_2))*ones(g,1)+O(g+1)-reply*(transpose(lambda_lm)*ones(c,1));
subject to
 lambda_lm>=0;
 lambdak_m>=0;
 lambdak_M>=0;
 lambda_k_m>=0;
 lambda_k_M>=0;
 mu_k_M>=0;
 mu_k_m>=0;
 R(:,:,g+1)>=0;
S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1)>=0:A;
cvx_end
fprintf('\nLagrange Multiplier of line real power flow limit constraints\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambda_lm,'Sheet1','A2:A3000');
fprintf('\nLagrange Multiplier of lower real power limit constraints\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambdak_m,'Sheet2','A2:A3100');
fprintf('\nLagrange Multiplier of upper real power limit constraints\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambdak_M,'Sheet2','B2:B3100');
fprintf('\nLagrange Multiplier of lower reactive power limit constraints\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambda_k_m,'Sheet2','C2:C3100');
fprintf('\nLagrange Multiplier of upper reactive power limit constraints\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambda_k_M,'Sheet2','D2:D3100');
fprintf('\nLagrange Multiplier of upper voltage limit constraints\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',mu_k_M,'Sheet2','E2:E3100');
fprintf('\nLagrange Multiplier of lower voltage limit constraints\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',mu_k_m,'Sheet2','F2:F3100');
fprintf('\nNumber of Iterations\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',cvx_slvitr,'Sheet3','B2');
fprintf('\nSolution Tolerance\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',cvx_slvtol,'Sheet3','B3');
fprintf('\nOptimal Dual Value\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',cvx_optval,'Sheet3','B4');
fprintf('\nThe following are the eigenvalues of the network matrix\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',eig(S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1)),'Sheet4','A2:A6100');
fprintf('\nThe following is the minimum eigenvalue of the network matrix\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',min(eig(S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1))),'Sheet3','B5');
fprintf('\nThe following are the eigenvalues of the dual matrix\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',eig(A),'Sheet4','B2:B6100');
[V,D]=eig(S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1));
fprintf('\nThese are the Eigenvectors and Eigenvalues of the Network Matrix\n');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',[V,D],'Sheet5','A1:ZK3000');
for j=1:4
    fprintf('\nThese are the values corresponding to the %d th eigenvector\n',j);
for i=1:n
    if mu_k_M(i)>=0.001
        Si_2=(1.06*V(n+1,j))/sqrt((((V(i,j))^2)+((V(i+n,j))^2))*(((V(1,j))^2)+((V(n+1,j))^2)));
        Si_1=-((Si_2)*V(1,j))/V(n+1,j);
        for f=1:n
            V_comp(f)=(Si_1)*V(f,j)-(Si_2)*V(f+n,j);
            V_comp(f+n)=(Si_2)*V(f,j)+(Si_1)*V(f+n,j);
            V_mod(f)=sqrt(((V_comp(f))^2)+((V_comp(f+n))^2));
        end
        W=((V_comp)')*V_comp;
for s=1:n
P(s)=trace(Y_n_a(:,:,s)*W);
Q(s)=trace(Y_n_r(:,:,s)*W);
end
P_Sol=0;
d=1;
for e=1:n
    if G_conn(e)==1
        P_Sol=P_Sol+c2(d)*((P(e)+P_d_k(e))^2)+c1(d)*(P(e)+P_d_k(e));
        d=d+1;
    end
end
fprintf('\nValues corresponding to binding constraint at bus %d\n',i);
fprintf('\nThese are the real and imaginary components of bus voltages\n');
disp(V_comp');
fprintf('\nThese are the bus voltage magnitudes\n');
disp(V_mod');
fprintf('\nReal Power Injection\n');
disp(P);
fprintf('\nReactive Power Injection\n');
disp(Q);
fprintf('\nGeneration Cost\n');
disp(P_Sol);
for e=1:n
    if G_conn(e)==1
        fprintf('\nGeneration at Bus %d\n',e);
        disp(P(e)+P_d_k(e));
    end
end
d=1;
for t=1:n
    if G_conn(t)==0
        LMPP(t)=(lambdak_M(t)-lambdak_m(t))/100;
    elseif G_conn(t)==1
            LMPP(t)=(lambdak_M(t)-lambdak_m(t)+2*(P(t)+P_d_k(t))*c2(d)+c1(d))/100;
            d=d+1;
    end
    LMPQ(t)=(lambda_k_M(t)-lambda_k_m(t))/100;
end
fprintf('\nReal Power LMP\n');
disp(LMPP');
fprintf('\nReactive Power LMP\n');
disp(LMPQ');
lcount=0;
for l=1:n
    for m=1:n
        if Yconn(l,m)==1
           lcount=lcount+1;
           P_flow(lcount)=trace(Y_line_net(:,:,lcount)*W);
           fprintf('\nThe real power flow from bus %d to bus %d on line %d is the following\n',l,m,lcount);
           disp(P_flow(lcount));
        end
    end
end

    end
end
end