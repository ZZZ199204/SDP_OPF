clc
clear;
n=input('Please Enter the particular IEEE Test Case You wish to run this Simulation for\nby indicating the number of buses.\nWarning:Valid Test Cases consist of 14, 30, 57, 118 & 300 buses.\n');
%Read the Bus Admittance Matrix, Matrices of line conductance, susceptance
%& shunt charging susceptance, Adjacency Matrix, Real & Reactive Power
%upper & lower limits, Real & Reactive Load data, Generator Location data,
%Generator Cost Coefficients and Number of Generators.
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
elseif n==14
    Y_bus_r=xlsread('14bus_Ybus.xlsx','Real','A1:N14');
    Y_bus_i=xlsread('14bus_Ybus.xlsx','Imag','A1:N14');
    Yline_G=xlsread('14bus_G.xlsx','Original','B2:O15');
    Yline_B=xlsread('14bus_B.xlsx','Sheet1','B2:O15');
    Yshunt_B=xlsread('14bus_shuntB.xlsx','Sheet1','B2:O15');
    Yconn=xlsread('14bus_Connectivity.xlsx','Sheet1','B2:O15');
    P_k_Max=[3.324;1.4;0;0;0;0;0;0;0;0;0;0;0;0];
    P_k_Min=[0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    Q_k_Max=[.10;.5;.4;0;0;.24;0;.24;0;0;0;0;0;0];
    Q_k_Min=[-.2;-.4;0;0;0;-.06;0;-.06;0;0;0;0;0;0];
    P_d_k=[0;.217;.942;.478;.076;.112;0;0;.295;.09;.035;.061;.135;.149];
    Q_d_k=[0;.127;.190;-.039;.016;.075;0;0;-.024;.058;.018;.016;.058;.05];
    G_conn=[1;1;0;0;0;0;0;0;0;0;0;0;0;0];
    c1=[2000;2000];
    c2=[430.293;2500];
    c0=[0;0];
    g=2;
    
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
%Maximum and Minimum Bus Voltage Limits.
V_k_Max=1.06*ones(n,1);
V_k_Min=0.94*ones(n,1);
V_k_M=V_k_Max.*V_k_Max;
V_k_m=V_k_Min.*V_k_Min;
B=eye(n);
%Computation of Network Matrices.
B_Count=1:1:n;
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
            K_line(:,:,c)=[transpose(Y_line_r(:,:,c)) zeros(n,n);zeros(n,n) -transpose(Y_line_i(:,:,c))];
            K_line_p(:,:,c)=[transpose(Y_line_i(:,:,c)) zeros(n,n);zeros(n,n) transpose(Y_line_r(:,:,c))];
            Theta_line(:,:,c)=K_line(:,:,c)*[B(:,l); B(:,l)]*[B(l,:) B(l,:)]*transpose(K_line(:,:,c))+K_line_p(:,:,c)*[B(:,l); B(:,l)]*[B(l,:) B(l,:)]*transpose(K_line_p(:,:,c));
            M1(:,:,c)=[(B(:,l)-B(:,m))*(B(l,:)-B(m,:)) zeros(n,n);zeros(n,n) (B(:,l)-B(:,m))*(B(l,:)-B(m,:))];
        end
    end
end
reply=input('Enter the pu value of the real power line limit\n');
voltage=input('Enter the pu value of the magnitude of the voltage difference allowed.\n(The allowed range is 0pu to 3.2pu).\n');
current=input('Enter the pu value of the line current flow limit\n');
fprintf('Dont close the window or open the Excel Spreadsheets to view results\nTill the time you are prompted to do so.\nComputation is in Progress.\nPlease Wait.');
%Solution of the Rank Relaxed Dual OPF.
cvx_begin sdp
variables lambdak_m(n) lambdak_M(n)  lambda_k_m(n) lambda_k_M(n) mu_k_M(n) mu_k_m(n) lambda_lm(c) mu_lm(c) Theta_lm(c);
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
    s(:,:,lcount+1)=s(:,:,lcount)+(lambda_lm(lcount))*Y_line_net(:,:,lcount)+(mu_lm(lcount))*M1(:,:,lcount)+(Theta_lm(lcount))*Theta_line(lcount);
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
maximize transpose(lambdak_m)*P_k_Min-transpose(lambdak_M)*P_k_Max+transpose(lambdak_M-lambdak_m)*P_d_k+transpose(lambda_k_m)*Q_k_Min-transpose(lambda_k_M)*Q_k_Max+transpose(lambda_k_M-lambda_k_m)*Q_d_k+transpose(mu_k_m)*(V_k_m)-transpose(mu_k_M)*(V_k_M)+(transpose(c0-r_2))*ones(g,1)+O(g+1)-reply*(transpose(lambda_lm)*ones(c,1))-(voltage^2)*(transpose(mu_lm)*ones(c,1))-current*(transpose(Theta_lm)*ones(c,1));
subject to
 lambda_lm>=0;
 lambdak_m>=0;
 lambdak_M>=0;
 lambda_k_m>=0;
 lambda_k_M>=0;
 mu_k_M>=0;
 mu_k_m>=0;
 mu_lm>=0;
 Theta_lm>=0;
 R(:,:,g+1)>=0;
S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1)>=0:A;
cvx_end
%If Primal OPF is Feasible or the Dual OPF is not unbounded, generates
%Excel Spreadsheets containing the Lagrange Multipliers, Optimal Value,
%Iteration count, Tolerance, Eigenvalues and Eigenvectors of the Network
%Matrix and Dual Matrix at the Optimum.
if isinf(cvx_optval)==0
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambda_lm,'Sheet1','A2:A3000');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',mu_lm,'Sheet1','B2:B3000');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',Theta_lm,'Sheet1','C2:C3000');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambdak_m,'Sheet2','A2:A3100');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambdak_M,'Sheet2','B2:B3100');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambda_k_m,'Sheet2','C2:C3100');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambda_k_M,'Sheet2','D2:D3100');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',mu_k_M,'Sheet2','E2:E3100');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',mu_k_m,'Sheet2','F2:F3100');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',cvx_slvitr,'Sheet3','B2');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',cvx_slvtol,'Sheet3','B3');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',cvx_optval,'Sheet3','B4');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',eig(S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1)),'Sheet4','A2:A6100');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',min(eig(S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1))),'Sheet3','B5');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',eig(A),'Sheet4','B2:B6100');
[V,D]=eig(S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1));
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',[V,D],'Sheet5','A1:ZK3000');
%Recovers the Optimal values of the Primal OPF following the zero angle
%condition at the slack bus (bus#1) and the KKT Optimality Condition at the
%bus at which the upper voltage limit constraint is tight based on the four
%smallest Eigenvalues of the Network Matrix.
z=0;
for j=1:4
for i=1:n
    if mu_k_M(i)>=0.001
        z=z+1;
        Si_2=(1.06*V(n+1,j))/sqrt((((V(i,j))^2)+((V(i+n,j))^2))*(((V(1,j))^2)+((V(n+1,j))^2)));
        Si_1=-((Si_2)*V(1,j))/V(n+1,j);
        for f=1:n
            %real and imaginary components of bus voltages
            V_comp(f,z,j)=(Si_1)*V(f,j)-(Si_2)*V(f+n,j);
            V_comp(f+n,z,j)=(Si_2)*V(f,j)+(Si_1)*V(f+n,j);
            %Absolute Values of Bus Voltages
            V_mod(f,z,j)=sqrt(((V_comp(f,z,j))^2)+((V_comp(f+n,z,j))^2));
        end
        W=(V_comp(:,z,j))*(V_comp(:,z,j))';
for s=1:n
    %Real and Reactive Bus Power Injections
P(s,z,j)=trace(Y_n_a(:,:,s)*W);
Q(s,z,j)=trace(Y_n_r(:,:,s)*W);
end
%Total Generation Cost at Optimum (Primal Optimal Objective Value)
P_Sol(z,j)=0;
d=1;
for e=1:n
    if G_conn(e)==1
        P_Sol(z,j)=P_Sol(z,j)+c2(d)*((P(e,z,j)+P_d_k(e))^2)+c1(d)*(P(e,z,j)+P_d_k(e));
        d=d+1;
    end
end
for e=1:n
    if G_conn(e)==1
        Gen(e,z,j)=(P(e,z,j)+P_d_k(e));
    elseif G_conn(e)==0
        Gen(e,z,j)=0;
    end
end
%Computation of Real & Reactive Power LMPs.
d=1;
for t=1:n
    if G_conn(t)==0
        LMPP(t,z,j)=(lambdak_M(t)-lambdak_m(t))/100;
    elseif G_conn(t)==1
            LMPP(t,z,j)=(lambdak_M(t)-lambdak_m(t)+2*(P(t,z,j)+P_d_k(t))*c2(d)+c1(d))/100;
            d=d+1;
    end
    LMPQ(t,z,j)=(lambda_k_M(t)-lambda_k_m(t))/100;
end
%Computation of Line Real Power Flows.
lcount=0;
for l=1:n
    for m=1:n
        if Yconn(l,m)==1
           lcount=lcount+1;
           P_flow(lcount,z,j)=trace(Y_line_net(:,:,lcount)*W);
           L(lcount,:)=[l m lcount];
   
        end
    end
end

    end
end
end

%Generation of Excel Spreadsheets for the computed primal optimal values.
v_comp=reshape(V_comp,2*n,j*z);
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',v_comp,'Sheet1','A1:ZK3000');
v_mod=reshape(V_mod,n,j*z);
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',v_mod,'Sheet2','A1:ZK3000');
p=reshape(P,n,j*z);
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',p,'Sheet3','A1:ZK3000');
q=reshape(Q,n,j*z);
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',q,'Sheet4','A1:ZK3000');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',P_Sol,'Sheet5','A1:ZK3000');
lmpp=reshape(LMPP,n,j*z);
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',lmpp,'Sheet6','A1:ZK3000');
lmpq=reshape(LMPQ,n,j*z);
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',lmpq,'Sheet7','A1:ZK3000');
gen=reshape(Gen,n,j*z);
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',gen,'Sheet8','A1:ZK3000');
p_flow=reshape(P_flow,lcount,j*z);
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',p_flow,'Sheet9','A1:ZK3000');
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',L,'Sheet10','A1:ZK3000');
%Graphing LMPs, Bus Voltage Magnitudes, Real and Reactive Demands and Bus
%Injections.
subplot(2,2,1);
plot(B_Count,lmpp(:,1),B_Count,lmpq(:,1));
title('Real & Reactive Power LMPs($/MWh)');
xlabel('Bus Number');
ylabel('LMP($/MWh)');
subplot(2,2,2);
plot(B_Count,v_mod(:,1),B_Count,V_k_Max,B_Count,V_k_Min);
title('Bus Voltage Limits & Magnitudes(pu)');
xlabel('Bus Number');
ylabel('pu Voltages');
subplot(2,2,3);
plot(B_Count,q(:,1),B_Count,Q_d_k);
title('Reactive Power Demands & Injections(pu)');
xlabel('Bus Number');
ylabel('Reactive Powers(pu)');
subplot(2,2,4);
plot(B_Count,p(:,1),B_Count,P_d_k);
title('Real Power Demands & Injections(pu)');
xlabel('Bus Number');
ylabel('Real Powers(pu)');
fprintf('All Computations done!!! Now you can open the Excel Spreadsheets to view the results.');
elseif (isinf(cvx_optval)==1)
    fprintf('Primal OPF is Infeasible for this problem instance');
end