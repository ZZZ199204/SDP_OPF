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
%%Primal Interior Point Algorithm with Closest Low Rank Approximation
it_ind=0;
t_avg=10000;
a_in=rand((2*n),(2*n));
A=a_in'*a_in;
outer_it_count=0;
while abs(IP_Stop(A,t_avg))>it_ind*(10^(-5))
    stop=IP_Stop(A,t_avg);
    t_avg=(it_ind*10+1)*t_avg;
    outer_it_count=outer_it_count+1;
    fprintf('Tolerance for outer loop is %2.8f\n',abs(stop));
    disp(t_avg);
    %Calculation of Gradient of the Lagrangian
    delL=zeros(2*n,2*n);
    j=1;
    for e=1:n
        if G_conn(e)==1
            delL=sparse(delL)+sparse(2*c2(j)*(trace((Y_n_a{e})*A)+P_d_k(e))*(Y_n_a{e})+c1(j)*(Y_n_a{e}));
            j=j+1;
        end
    end
    for i=1:n
        delL=sparse(delL)-(1/t_avg)*(sparse((Y_n_a{i})*(1/(trace((Y_n_a{i})*A)-P_k_Max(i)+P_d_k(i))))+...
            sparse((Y_n_a{i})*(1/(trace((Y_n_a{i})*A)-P_k_Min(i)+P_d_k(i))))+...
            sparse((Y_n_r{i})*(1/(trace((Y_n_r{i})*A)-Q_k_Max(i)+Q_d_k(i))))+...
            sparse((Y_n_r{i})*(1/(trace((Y_n_r{i})*A)-Q_k_Min(i)+Q_d_k(i))))+...
            sparse((M{i})*(1/(trace((M{i})*A)-V_k_M(i))))+...
            sparse((M{i})*(1/(trace((M{i})*A)-V_k_m(i)))));
    end
    for c=1:lcount
        delL=sparse(delL)-(1/t_avg)*(sparse((Y_line_net{c})*(1/(trace((Y_line_net{c})*A)-reply)))+...
            sparse((M1{c})*(1/(trace((M1{c})*A)-(voltage^2)))));
    end
    delL=sparse(delL)-sparse((1/t_avg)*(inv(A)));
    %Calculation of Hessian of the Lagrangian
    delsqL=zeros(((2*n)^2),((2*n)^2));
    j=1;
    for e=1:n
        if G_conn(e)==1
            delsqL=sparse(delsqL)+sparse(2*c2(j)*kron(Y_n_a{e},Y_n_a{e}));
            j=j+1;
        end
    end
    for i=1:n
        delsqL=sparse(delsqL)+(1/t_avg)*(sparse((kron(Y_n_a{i},Y_n_a{i}))*(1/((trace((Y_n_a{i})*A)-P_k_Max(i)+P_d_k(i))^2)))+...
            sparse((kron(Y_n_a{i},Y_n_a{i}))*(1/((trace((Y_n_a{i})*A)-P_k_Min(i)+P_d_k(i))^2)))+...
            sparse((kron(Y_n_r{i},Y_n_r{i}))*(1/((trace((Y_n_r{i})*A)-Q_k_Max(i)+Q_d_k(i))^2)))+...
            sparse((kron(Y_n_r{i},Y_n_r{i}))*(1/((trace((Y_n_r{i})*A)-Q_k_Min(i)+Q_d_k(i))^2)))+...
            sparse((kron(M{i},M{i}))*(1/((trace((M{i})*A)-V_k_M(i))^2)))+...
            sparse((kron(M{i},M{i}))*(1/((trace((M{i})*A)-V_k_m(i))^2))));
    end
    for c=1:lcount
        delsqL=sparse(delsqL)+(1/t_avg)*(sparse((kron(Y_line_net{c},Y_line_net{c}))*(1/((trace((Y_line_net{c})*A)-reply)^2)))+...
            sparse((kron(M1{c},M1{c}))*(1/((trace((M1{c})*A)-(voltage^2))^2))));
    end
    delsqL=sparse(delsqL)+sparse(DiffInv(A));
    %Calculation of Improving Step Direction
    dsqL=sparse(transpose(im2col(delsqL,[(2*n) (2*n)],'distinct')));
    dL=sparse(reshape(delL,((2*n)^2),1));
    deltaA=sparse(-(inv(dsqL))*dL);
    %Newton Decrement Calculation
    Newt_decr_sq=dL'*(inv(dsqL))*dL;
    %fprintf('Tolerance for inner loop is %2.8f\n',(Newt_decr_sq/2));
    dA=sparse(reshape(deltaA,(2*n),(2*n)));
    inner_it_count=0;
    while (Newt_decr_sq/2)>10^(-3)
        %Calculation of Step Length: Backtracking Line Search
        stlength=1;
        while Lagrangian((A+(stlength*dA)),t_avg)>Lagrangian(A,t_avg)+0.25*stlength*trace(delL*dA)
                stlength=stlength*0.5;
        end
        A=sparse(A+stlength*dA);
        inner_it_count=inner_it_count+1;
        delL=sparse(Gradient(A,t_avg));
        delsqL=sparse(Hessian(A,t_avg));
        %Calculation of Improving Step Direction
        dsqL=sparse(transpose(im2col(delsqL,[(2*n) (2*n)],'distinct')));
        dL=sparse(reshape(delL,((2*n)^2),1));
        deltaA=sparse(-(inv(dsqL))*dL);
        %Newton Decrement Calculation
        Newt_decr_sq=dL'*(inv(dsqL))*dL;
        dA=sparse(reshape(deltaA,(2*n),(2*n)));
        disp(inner_it_count);
        disp((Newt_decr_sq/2));
    end
    it_ind=1;
    disp(outer_it_count);
    disp(abs(stop));
end
%Calculation of the closest rank-1 Approximate
[V,D]=eig(full(A));
D1=zeros((2*n),(2*n));
D1((2*n),(2*n))=D((2*n),(2*n));
A=V*D1*V';
%disp(A);
for s_count=1:n
    %Real and Reactive Bus Power Injections
    P(s_count)=trace(Y_n_a{s_count}*A);
    Q(s_count)=trace(Y_n_r{s_count}*A);
end
disp(P);
disp(Q);
        %Total Generation Cost at Optimum (Primal Optimal Objective Value)
        P_Sol=0;
        d=1;
        for e=1:n
            if G_conn(e)==1
                P_Sol=P_Sol+c2(d)*((P(e)+P_d_k(e))^2)+c1(d)*(P(e)+P_d_k(e));
                d=d+1;
            end
        end
        disp(P_Sol);
        for e=1:n
            if G_conn(e)==1
                Gen(e)=(P(e)+P_d_k(e));
            elseif G_conn(e)==0
                Gen(e)=0;
            end
        end
        disp(Gen);
        