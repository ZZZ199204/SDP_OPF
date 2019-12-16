clc
clear all;
global n Y_n_a Y_n_r Y_line_net reply voltage M M1 P_d_k ...
    P_k_Max P_k_Min Q_d_k Q_k_Max Q_k_Min G_conn c2 c1 V_k_M V_k_m lcount
t0=tic;
n=input('Please Enter the particular IEEE Test Case You wish to run this Simulation for\nby indicating the number of buses.\nWarning:Valid Test Cases consist of 3, 14, 30, 57, 118 & 300 buses.\n');
%% Read the Bus Admittance Matrix, Matrices of line conductance, susceptance
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
 elseif n==3
    Y_bus_r=xlsread('3bus_Ybus.xlsx','Sheet2','A1:C3');
    Y_bus_i=xlsread('3bus_Ybus.xlsx','Sheet3','A1:C3');
    Yline_G=xlsread('3bus_G.xlsx','Sheet1','A1:C3');
    Yline_B=xlsread('3bus_B.xlsx','Sheet1','A1:C3');
    Yshunt_B=xlsread('3bus_Bshunt.xlsx','Sheet1','A1:C3');
    Yconn=xlsread('3bus_Connectivity.xlsx','Sheet1','A1:C3');
    P_k_Max=xlsread('Limits.xlsx','3_Bus','B2');
    P_k_Min=xlsread('Limits.xlsx','3_Bus','C2');
    Q_k_Max=xlsread('Limits.xlsx','3_Bus','D2');
    Q_k_Min=xlsread('Limits.xlsx','3_Bus','E2');
    P_d_k=xlsread('3bus_load.xlsx','Sheet1','D2:D4');
    Q_d_k=xlsread('3bus_load.xlsx','Sheet1','E2:E4');
    G_conn=xlsread('Limits.xlsx','3_Bus','F2:F4');
    c1=xlsread('Limits.xlsx','3','I2');
    c2=xlsread('Limits.xlsx','3','H2');
    c0=xlsread('Limits.xlsx','3','G2');
    g=1;
else
    fprintf('Test case with this particular number of buses is not supported by the Program\n');
end
fprintf('Time taken to read %3.5f',toc(t0));
%Maximum and Minimum Bus Voltage Limits.
V_k_Max=1.06*ones(n,1);
V_k_Min=0.94*ones(n,1);
V_k_M=V_k_Max.*V_k_Max;
V_k_m=V_k_Min.*V_k_Min;
B=speye(n);
y_bus_r=sparse(Y_bus_r);
y_bus_i=sparse(Y_bus_i);
yline_G=sparse(Yline_G);
yline_B=sparse(Yline_B);
yshunt_B=sparse(Yshunt_B);
%% Computation of Network Matrices.
B_Count=1:1:n;
Y_n_a = cell(n,1);
Y_n_r = cell(n,1);
M=cell(n,1);
for j=1:n
    Y_r=B(:,j)*B(j,:)*y_bus_r;
    Y_i=B(:,j)*B(j,:)*y_bus_i;
    Y_n_a{j}= 0.5*[Y_r+transpose(Y_r) transpose(Y_i)-Y_i;
        Y_i-transpose(Y_i) Y_r+transpose(Y_r)];
    Y_n_r{j}=-0.5*[Y_i+transpose(Y_i) Y_r-transpose(Y_r);
        transpose(Y_r)-Y_r Y_i+transpose(Y_i)];
    M{j}=[B(:,j)*B(j,:) zeros(n,n);zeros(n,n) B(:,j)*B(j,:)];
end
fprintf('Node Matrix Computations Completed\n');
Y_line_net = cell(n*n, 1);
M1=cell(n*n,1);
c=0;
for l=1:n
    for m=1:n
        if Yconn(l,m)==1
            c=c+1;
            Y_line_r=yline_G(l,m)*B(:,l)*B(l,:)-yline_G(l,m)*B(:,l)*B(m,:);
            Y_line_i=(yline_B(l,m)+yshunt_B(l,m))*B(:,l)*B(l,:)-yline_B(l,m)*B(:,l)*B(m,:);
            Y_line_net{c}=0.5*[Y_line_r+transpose(Y_line_r) transpose(Y_line_i)-Y_line_i;Y_line_i-transpose(Y_line_i) Y_line_r+transpose(Y_line_r)];
            M1{c}=[(B(:,l)-B(:,m))*(B(l,:)-B(m,:)) zeros(n,n);zeros(n,n) (B(:,l)-B(:,m))*(B(l,:)-B(m,:))];
        end
    end
end
lcount=c;
fprintf('Line Matrix Computations Completed\n');
reply=input('Enter the pu value of the real power line limit\n');
voltage=input('Enter the pu value of the magnitude of the voltage difference allowed.\n(The allowed range is 0pu to 2.12pu).\n');
fprintf('Dont close the window or open the Excel Spreadsheets to view results\nTill the time you are prompted to do so.\nComputation is in Progress.\nPlease Wait.\n');
%% Solution of the Rank Relaxed Dual OPF.
cvx_begin sdp
variables lambdak_m(n) lambdak_M(n)  lambda_k_m(n) lambda_k_M(n) mu_k_M(n) mu_k_m(n) ...
    lambda_lm(lcount) mu_lm(lcount);
variables r_1(g) r_2(g);
dual variable A;
S=cell(n+1,1);
S{1}=zeros(2*n,2*n);
for count=1:n
    S{count+1}=sparse(S{count})+...
        sparse((lambdak_M(count)-lambdak_m(count))*Y_n_a{count})+...
        sparse((lambda_k_M(count)-lambda_k_m(count))*Y_n_r{count})+...
        sparse((mu_k_M(count)-mu_k_m(count))*M{count});
end
fprintf('End of 1st cvx loop\n');
s=cell(lcount+1,1);
s{1}=zeros(2*n,2*n);
for c=1:lcount
    s{c+1}=sparse(s{c})+...
        sparse((lambda_lm(c))*Y_line_net{c})+...
        sparse((mu_lm(c))*M1{c});
end
fprintf('End of 2nd cvx loop\n');
R=cell(g+1,1);
expression T(2*g,2*g,g);
R{1}=zeros(2*g,2*g);
for i=1:g
    T(:,:,i)=zeros(2*g,2*g);
    T(2*i-1,2*i-1,i)=T(2*i-1,2*i-1,i)+1;
    T(2*i-1,2*i,i)=T(2*i-1,2*i,i)+r_1(i);
    T(2*i,2*i-1,i)=T(2*i,2*i-1,i)+r_1(i);
    T(2*i,2*i,i)=T(2*i,2*i,i)+r_2(i);
end
fprintf('End of 3rd cvx loop\n');
for gcount=1:g
  R{gcount+1}=sparse(R{gcount})+sparse(T(:,:,gcount));
end
fprintf('End of 4th cvx loop\n');
N=cell(g+1,1);
N{1}=zeros(2*n,2*n);
j=1;
for i=1:n
    fprintf('%3d\n',i);
    if G_conn(i)==1
        N{j+1}=sparse(N{j})+sparse((c1(j)+2*sqrt(c2(j))*r_1(j))*Y_n_a{i});
        j=j+1;
        fprintf('%3d\t%3d\n',j,i);
    end
end
fprintf('End of 5th cvx loop\n');
expression O(g+1);
O(1)=0;
j=1;
for i=1:n
    if G_conn(i)==1
        O(j+1)=O(j)+(c1(j)+2*sqrt(c2(j))*r_1(j))*P_d_k(i);
        j=j+1;
    end
end
fprintf('End of 6th cvx loop\n');
maximize transpose(lambdak_m)*P_k_Min-transpose(lambdak_M)*P_k_Max+transpose(lambdak_M-lambdak_m)*P_d_k+...
    transpose(lambda_k_m)*Q_k_Min-transpose(lambda_k_M)*Q_k_Max+transpose(lambda_k_M-lambda_k_m)*Q_d_k+...
    transpose(mu_k_m)*(V_k_m)-transpose(mu_k_M)*(V_k_M)+(transpose(c0-r_2))*ones(g,1)+O(g+1)-...
    reply*(transpose(lambda_lm)*ones(lcount,1))-(voltage^2)*(transpose(mu_lm)*ones(lcount,1));
subject to
 lambda_lm>=0;
 lambdak_m>=0;
 lambdak_M>=0;
 lambda_k_m>=0;
 lambda_k_M>=0;
 mu_k_M>=0;
 mu_k_m>=0;
 mu_lm>=0;
 R{g+1}>=0;
 S{n+1}+s{lcount+1}+N{g+1}>=0:A;
cvx_end
%% If Primal OPF is Feasible or the Dual OPF is not unbounded, generates
%Excel Spreadsheets containing the Lagrange Multipliers, Optimal Value,
%Iteration count, Tolerance, Eigenvalues and Eigenvectors of the Network
%Matrix and Dual Matrix at the Optimum.
if isinf(cvx_optval)==0
    Mult=[lambdak_m';lambdak_M';lambda_k_m';lambda_k_M';mu_k_M';mu_k_m'];
    Eig=[(eig(full(S{n+1}+s{lcount+1}+N{g+1})))';(eig(full(A)))'];
    fileID=fopen('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30_volt.txt','w');
    fprintf(fileID,'%60s\t','Number of Bus,MW Limit,pu Voltage Diff. Magnitude Limit:');
    fprintf(fileID,'%3d\t%4.2f\t%2.7f\r\n',n,reply*100,voltage);
    fprintf(fileID,'%11s\t%11s\t%11s\t%11s\t%11s\t%11s\r\n','lambdak_m','lambdak_M','lambda_k_m','lambda_k_M','mu_k_M','mu_k_m');
    fprintf(fileID,'%5.5f\t%5.5f\t%5.5f\t%5.5f\t%5.5f\t%5.5f\r\n',Mult);
    fprintf(fileID,'%11s\t%12s\r\n','N/W Mat Eig','Dual Mat Eig');
    fprintf(fileID,'%6.5f\t%6.5f\r\n',Eig);
    fprintf(fileID,'%12s','Iterations:');
    fprintf(fileID,'%3d\r\n',cvx_slvitr);
    fprintf(fileID,'%12s','Tolerance:');
    fprintf(fileID,'%6.7f\r\n',cvx_slvtol);
    fprintf(fileID,'%12s','Dual Optval:');
    fprintf(fileID,'%7.8f\r\n',cvx_optval);
    fprintf(fileID,'%12s','Min Eigen:');
    fprintf(fileID,'%7.8f\r\n',min(eig(full(S{n+1}+s{lcount+1}+N{g+1}))));
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambda_lm,'Sheet1','A2:A300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',mu_lm,'Sheet1','B2:B300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambdak_m,'Sheet2','A2:A310');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambdak_M,'Sheet2','B2:B310');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambda_k_m,'Sheet2','C2:C310');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',lambda_k_M,'Sheet2','D2:D310');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',mu_k_M,'Sheet2','E2:E310');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',mu_k_m,'Sheet2','F2:F310');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',cvx_slvitr,'Sheet3','B2');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',cvx_slvtol,'Sheet3','B3');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',cvx_optval,'Sheet3','B4');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',eig(S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1)),'Sheet4','A2:A610');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',min(eig(S(:,:,n+1)+s(:,:,c+1)+N(:,:,g+1))),'Sheet3','B5');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',eig(A),'Sheet4','B2:B610');
    [V,D]=eig(full(S{n+1}+s{lcount+1}+N{g+1}));
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\out30.xlsx',[V,D],'Sheet5','A1:AZ300');
    %% Recovers the Optimal values of the Primal OPF following the zero angle
    %condition at the slack bus (bus#1) and the KKT Optimality Condition at the
    %bus at which the upper voltage limit constraint is tight based on the four
    %smallest Eigenvalues of the Network Matrix.
    for j=1:4
        [C,I]=max(mu_k_M);
        Si_2=(1.06*V(n+1,j))/sqrt((((V(I,j))^2)+((V(I+n,j))^2))*(((V(1,j))^2)+((V(n+1,j))^2)));
        Si_1=-((Si_2)*V(1,j))/V(n+1,j);
        for f=1:n
            %real and imaginary components of bus voltages
            V_comp(f,j)=(Si_1)*V(f,j)-(Si_2)*V(f+n,j);
            V_comp(f+n,j)=(Si_2)*V(f,j)+(Si_1)*V(f+n,j);
            %Absolute Values of Bus Voltages
            V_mod(f,j)=sqrt(((V_comp(f,j))^2)+((V_comp(f+n,j))^2));
        end
        W=(V_comp(:,j))*(V_comp(:,j))';
        for s_count=1:n
            %Real and Reactive Bus Power Injections
            P(s_count,j)=trace(Y_n_a{s_count}*W);
            Q(s_count,j)=trace(Y_n_r{s_count}*W);
        end
        %Total Generation Cost at Optimum (Primal Optimal Objective Value)
        P_Sol(j)=0;
        d=1;
        for e=1:n
            if G_conn(e)==1
                P_Sol(j)=P_Sol(j)+c2(d)*((P(e,j)+P_d_k(e))^2)+c1(d)*(P(e,j)+P_d_k(e));
                d=d+1;
            end
        end
        for e=1:n
            if G_conn(e)==1
                Gen(e,j)=(P(e,j)+P_d_k(e));
            elseif G_conn(e)==0
                Gen(e,j)=0;
            end
        end
        %Computation of Real & Reactive Power LMPs.
        d=1;
        for t=1:n
            if G_conn(t)==0
                LMPP(t,j)=(lambdak_M(t)-lambdak_m(t))/100;
            elseif G_conn(t)==1
                LMPP(t,j)=(lambdak_M(t)-lambdak_m(t)+2*(P(t,j)+P_d_k(t))*c2(d)+c1(d))/100;
                d=d+1;
            end
            LMPQ(t,j)=(lambda_k_M(t)-lambda_k_m(t))/100;
        end
        %Computation of Line Real Power Flows.
        c=0;
        for l=1:n
            for m=1:n
                if Yconn(l,m)==1
                    c=c+1;
                    P_flow(c,j)=trace(Y_line_net{c}*W);
                    Ang(c,j)=(180/pi)*(abs(atan(V_comp(l+n,j)/V_comp(l,j))-atan(V_comp(m+n,j)/V_comp(m,j))));
                    L(c,:)=[l m c];  
                end
            end
        end
    end
    %% Generation of Excel Spreadsheets for the computed primal optimal values.
    fprintf(fileID,'%50s\t','Binding Const. Upper Voltage Limit at Bus:');
    fprintf(fileID,'%3d\r\n',I);
    fprintf(fileID,'%40s\r\n','Bus Voltage Components');
    fprintf(fileID,'%7s\t%7s\t%7s\t%7s\r\n','1st Eig','2nd Eig','3rd Eig','4th Eig');
    fprintf(fileID,'%2.3f\t%2.3f\t%2.3f\t%2.3f\r\n',V_comp');
    fprintf(fileID,'%40s\r\n','Bus Voltage Magnitudes');
    fprintf(fileID,'%7s\t%7s\t%7s\t%7s\r\n','1st Eig','2nd Eig','3rd Eig','4th Eig');
    fprintf(fileID,'%2.3f\t%2.3f\t%2.3f\t%2.3f\r\n',V_mod');
    fprintf(fileID,'%40s\r\n','Real Power Injections');
    fprintf(fileID,'%7s\t%7s\t%7s\t%7s\r\n','1st Eig','2nd Eig','3rd Eig','4th Eig');
    fprintf(fileID,'%3.3f\t%3.3f\t%3.3f\t%3.3f\r\n',P');
    fprintf(fileID,'%40s\r\n','Reactive Power Injections');
    fprintf(fileID,'%7s\t%7s\t%7s\t%7s\r\n','1st Eig','2nd Eig','3rd Eig','4th Eig');
    fprintf(fileID,'%3.3f\t%3.3f\t%3.3f\t%3.3f\r\n',Q');
    fprintf(fileID,'%40s\r\n','Real Power LMP');
    fprintf(fileID,'%7s\t%7s\t%7s\t%7s\r\n','1st Eig','2nd Eig','3rd Eig','4th Eig');
    fprintf(fileID,'%4.3f\t%4.3f\t%4.3f\t%4.3f\r\n',LMPP');
    fprintf(fileID,'%40s\r\n','Reactive Power LMP');
    fprintf(fileID,'%7s\t%7s\t%7s\t%7s\r\n','1st Eig','2nd Eig','3rd Eig','4th Eig');
    fprintf(fileID,'%4.3f\t%4.3f\t%4.3f\t%4.3f\r\n',LMPQ');
    fprintf(fileID,'%40s\r\n','Generation');
    fprintf(fileID,'%7s\t%7s\t%7s\t%7s\r\n','1st Eig','2nd Eig','3rd Eig','4th Eig');
    fprintf(fileID,'%4.3f\t%4.3f\t%4.3f\t%4.3f\r\n',Gen');
    Flow=[L';lambda_lm';mu_lm';P_flow';Ang'];
    fprintf(fileID,'%10s\r\n','Flows');
    fprintf(fileID,'%4s\t%3s\t%4s\t%10s\t%10s\t%17s\t%17s\t%17s\t%17s\t%20s\t%20s\t%20s\t%20s\r\n','From','To',...
        'Line','Lambda_lm','mu_lm','Flow(1st Eig)','Flow(2nd Eig)','Flow(3rd Eig)','Flow(4th Eig)','Volt_Ang(1st Eig)',...
        'Volt_Ang(2nd Eig)','Volt_Ang(3rd Eig)','Volt_Ang(4th Eig)');
    fprintf(fileID,'%3d\t%3d\t%3d\t%5.3f\t%5.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\t%4.3f\r\n',Flow);
    fprintf(fileID,'%40s\r\n','Optimal Primal');
    fprintf(fileID,'%9s\t%9s\t%9s\t%9s\r\n','1st Eig','2nd Eig','3rd Eig','4th Eig');
    fprintf(fileID,'%6.3f\t%6.3f\t%6.3f\t%6.3f\r\n',P_Sol');
    fclose(fileID);
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',V_comp,'Sheet1','A1:Z300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',V_mod,'Sheet2','A1:Z300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',P,'Sheet3','A1:Z300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',Q,'Sheet4','A1:Z300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',P_Sol,'Sheet5','A1:Z300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',LMPP,'Sheet6','A1:Z300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',LMPQ,'Sheet7','A1:Z300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',Gen,'Sheet8','A1:Z300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',P_flow,'Sheet9','A1:Z300');
    %xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\deux.xlsx',L,'Sheet10','A1:Z300');
    %% Graphing LMPs, Bus Voltage Magnitudes, Real and Reactive Demands and Bus
    %Injections.
    L_Count=1:1:lcount;
    ang_const=90*ones(lcount,1);
    Reply=reply*ones(lcount,1);
    neg_reply=(-reply)*ones(lcount,1);
    subplot(4,2,5);
    plot(B_Count,LMPP(:,1),B_Count,LMPQ(:,1));
    title('Real & Reactive Power LMPs($/MWh)');
    xlabel('Bus Number');
    ylabel('LMP($/MWh)');
    subplot(4,2,6);
    plot(B_Count,V_mod(:,1),B_Count,V_k_Max,B_Count,V_k_Min);
    title('Bus Voltage Limits & Magnitudes(pu)');
    xlabel('Bus Number');
    ylabel('pu Voltages');
    subplot(4,2,7);
    plot(B_Count,Q(:,1),B_Count,Q_d_k);
    title('Reactive Power Demands & Injections(pu)');
    xlabel('Bus Number');
    ylabel('Reactive Powers(pu)');
    subplot(4,2,8);
    plot(B_Count,P(:,1),B_Count,P_d_k);
    title('Real Power Demands & Injections(pu)');
    xlabel('Bus Number');
    ylabel('Real Powers(pu)');
    subplot(4,2,3:4);
    plot(L_Count,Ang(:,1),L_Count,ang_const);
    title('Voltage Angle Differences in Degrees');
    xlabel('Line Number');
    ylabel('Voltage Angle Difference(Degrees)');
    fprintf('All Computations done!!! Now you can open the Excel Spreadsheets to view the results.');
    subplot(4,2,1);
    plot(L_Count,P_flow(:,1),L_Count,Reply,L_Count,neg_reply);
    title('Line Flows in pu');
    xlabel('Line Number');
    ylabel('Line Flows(pu)');
    fprintf('All Computations done!!! Now you can open the Excel Spreadsheets to view the results.');
    subplot(4,2,2);
    plot(B_Count,Gen(:,1),B_Count,P_k_Max,B_Count,P_k_Min);
    title('Real Power Generation (pu)');
    xlabel('Bus Number');
    ylabel('Real Power Generation (pu)');
elseif (isinf(cvx_optval)==1)
    fprintf('Primal OPF is Infeasible for this problem instance');
end
%%Primal Interior Point Algorithm with Closest Low Rank Approximation 
%Calculation for Tuning Parameter
%A=(V_comp(:,1))*(V_comp(:,1))';
for i=1:n
    t1(i)=-1/((lambdak_M(i))*(trace((Y_n_a{i})*A)-P_k_Max(i)+P_d_k(i)));
    t2(i)=1/((lambdak_m(i))*(trace((Y_n_a{i})*A)-P_k_Min(i)+P_d_k(i)));
    t3(i)=-1/((lambda_k_M(i))*(trace((Y_n_r{i})*A)-Q_k_Max(i)+Q_d_k(i)));
    t4(i)=1/((lambda_k_m(i))*(trace((Y_n_r{i})*A)-Q_k_Min(i)+Q_d_k(i)));
    t5(i)=-1/((mu_k_M(i))*(trace((M{i})*A)-V_k_M(i)));
    t6(i)=1/((mu_k_m(i))*(trace((M{i})*A)-V_k_m(i)));
end
t_nw=[t1' t2' t3' t4' t5' t6'];
for c=1:lcount
    t7(c)=-1/((lambda_lm(c))*(trace(Y_line_net{c}*A)-reply));
    t8(c)=-1/((lambda_lm(c))*(trace(M1{c}*A)-(voltage^2)));
end
t_lin=[t7' t8'];
t9=A*(S{n+1}+s{lcount+1}+N{g+1});
disp(t9);
t10=1/(trace(A*(S{n+1}+s{lcount+1}+N{g+1}))/(2*n));
t_avg=((t1+t2+t3+t4+t5+t6)*ones(n,1)+(t7+t8)*ones(lcount,1)+t10)/(6*n+2*lcount+1);
fprintf(fileID,'%25s\r\n','Tuning Parameters_nw');
fprintf(fileID,'%6s\t%6s\t%6s\t%6s\t%6s\t%6s\r\n','PMax','PMin','QMax','QMin','VMax','VMin');
fprintf(fileID,'%7.8f\t%7.8f\t%7.8f\t%7.8f\t%7.8f\t%7.8f\r\n',t_nw);
fprintf(fileID,'%35s\r\n','Tuning Parameters_line');
fprintf(fileID,'%6s\t%6s\r\n','Plm','Vlm');
fprintf(fileID,'%7.8f\t%7.8f\r\n',t_lin);
fprintf(fileID,'%25s\r\n','Tuning Parameters_mat');
fprintf(fileID,'%7.8f\t%7.8f\r\n',t10,t_avg);
%Calculation of Gradient of the Lagrangian
delL=zeros(2*n,2*n);
j=1;
for e=1:n
    if G_conn(e)==1
        delL=delL+2*c2(j)*(trace((Y_n_a(:,:,e))*A)+P_d_k(e))*(Y_n_a(:,:,e))+c1(j)*(Y_n_a(:,:,e));
        j=j+1;
    end
end
for i=1:n
    delL=delL-(1/t_avg)*((Y_n_a(:,:,i))*(1/(trace((Y_n_a(:,:,i))*A)-P_k_Max(i)+P_d_k(i)))+...
        (Y_n_a(:,:,i))*(1/(trace((Y_n_a(:,:,i))*A)-P_k_Min(i)+P_d_k(i)))+...
        (Y_n_r(:,:,i))*(1/(trace((Y_n_r(:,:,i))*A)-Q_k_Max(i)+Q_d_k(i)))+...
        (Y_n_r(:,:,i))*(1/(trace((Y_n_r(:,:,i))*A)-Q_k_Min(i)+Q_d_k(i)))+...
        (M(:,:,i))*(1/(trace((M(:,:,i))*A)-V_k_M(i)))+...
        (M(:,:,i))*(1/(trace((M(:,:,i))*A)-V_k_m(i))));
end
for c=1:lcount
    delL=delL-(1/t_avg)*((Y_line_net(:,:,c))*(1/(trace((Y_line_net(:,:,c))*A)-reply))+...
        (M1(:,:,c))*(1/(trace((M1(:,:,c))*A)-(voltage^2))));
end
delL=delL-(1/t_avg)*(inv(A));
%Calculation of Hessian of the Lagrangian
delsqL=zeros(((2*n)^2),((2*n)^2));
j=1;
for e=1:n
    if G_conn(e)==1
        delsqL=delsqL+2*c2(j)*kron(Y_n_a(:,:,e),Y_n_a(:,:,e));
        j=j+1;
    end
end
for i=1:n
    delsqL=delsqL+(1/t_avg)*((kron(Y_n_a(:,:,i),Y_n_a(:,:,i)))*(1/((trace((Y_n_a(:,:,i))*A)-P_k_Max(i)+P_d_k(i))^2))+...
        (kron(Y_n_a(:,:,i),Y_n_a(:,:,i)))*(1/((trace((Y_n_a(:,:,i))*A)-P_k_Min(i)+P_d_k(i))^2))+...
        (kron(Y_n_r(:,:,i),Y_n_r(:,:,i)))*(1/((trace((Y_n_r(:,:,i))*A)-Q_k_Max(i)+Q_d_k(i))^2))+...
        (kron(Y_n_r(:,:,i),Y_n_r(:,:,i)))*(1/((trace((Y_n_r(:,:,i))*A)-Q_k_Min(i)+Q_d_k(i))^2))+...
        (kron(M(:,:,i),M(:,:,i)))*(1/((trace((M(:,:,i))*A)-V_k_M(i))^2))+...
        (kron(M(:,:,i),M(:,:,i)))*(1/((trace((M(:,:,i))*A)-V_k_m(i))^2)));
end
for c=1:lcount
    delsqL=delsqL+(1/t_avg)*((kron(Y_line_net(:,:,c),Y_line_net(:,:,c)))*(1/((trace((Y_line_net(:,:,c))*A)-reply)^2))+...
        (kron(M1(:,:,c),M1(:,:,c)))*(1/((trace((M1(:,:,c))*A)-(voltage^2))^2)));
end
delsqL=delsqL+DiffInv(A);
%Calculation of Improving Step Direction
dsqL=transpose(im2col(delsqL,[(2*n) (2*n)],'distinct'));
dL=reshape(delL,((2*n)^2),1);
deltaA=-(inv(dsqL))*dL;
%Newton Decrement Calculation
Newt_decr_sq=dL'*(inv(dsqL))*dL;
dA=reshape(deltaA,(2*n),(2*n));
%while (Newt_decr_sq/2)>10^(-13)
%Calculation of Step Length: Backtracking Line Search
stlength=1;
while Lagrangian((A+(stlength*dA)),t_avg)>Lagrangian(A,t_avg)+0.25*stlength*trace(delL*dA)
    stlength=stlength*0.5;
end
A=A+stlength*dA;
%delL=Gradient(A,t_avg);
%delsqL=Hessian(A,t_avg);
%Calculation of Improving Step Direction
%dsqL=transpose(im2col(delsqL,[(2*n) (2*n)],'distinct'));
%dL=reshape(delL,((2*n)^2),1);
%deltaA=-(inv(dsqL))*dL;
%Newton Decrement Calculation
%Newt_decr_sq=dL'*(inv(dsqL))*dL;
%dA=reshape(deltaA,(2*n),(2*n));
%end
%Calculation of the closest rank-1 Approximate
[V,D]=eig(A);
D1=zeros((2*n),(2*n));
D1((2*n),(2*n))=D((2*n),(2*n));
A=V*D1*V';
t_avg=10*t_avg;
disp(Newt_decr_sq);
disp(IP_Stop(A,t_avg));
fclose(fileID);