clc
clear;
Y_Bus=[6.03-19.45i -5+15.26i 0 0 -1.03+4.23i 0 0 0 0 0 0 0 0 0;
    -5+15.26i 9.52-30.27i -1.14+4.78i -1.69+5.12i -1.7+5.19i 0 0 0 0 0 0 0 0 0;
    0 -1.14+4.78i 3.12-9.82i -1.99+5.07i 0 0 0 0 0 0 0 0 0 0;
    0 -1.69+5.12i -1.99+5.07i 10.51-38.65i -6.84+21.58i 0 0+4.89i 0 0+1.86i 0 0 0 0 0;
    -1.03+4.23i -1.7+5.19i 0 -6.84+21.58i 9.57-35.53i 0+4.26i 0 0 0 0 0 0 0 0;
    0 0 0 0 0+4.26i 6.58-17.34i 0 0 0 0 -1.96+4.09i -1.53+3.18i -3.1+6.10i 0;
    0 0 0 0+4.89i 0 0 0-19.55i 0+5.68i 0+9.09i 0 0 0 0 0;
    0 0 0 0 0 0 0+5.68i 0-5.68i 0 0 0 0 0 0;
    0 0 0 0+1.86i 0 0 0+9.09i 0 5.33-24.09i -3.9+10.37i 0 0 0 -1.42+3.03i;
    0 0 0 0 0 0 0 0 -3.9+10.37i 5.78-14.77i -1.88+4.4i 0 0 0;
    0 0 0 0 0 -1.96+4.09i 0 0 0 -1.88+4.4i 3.84-8.5i 0 0 0;
    0 0 0 0 0 -1.53+3.18i 0 0 0 0 0 4.01-5.43i -2.49+2.25i 0;
    0 0 0 0 0 -3.1+6.10i 0 0 0 0 0 -2.49+2.25i 6.72-10.67i -1.14+2.31i;
    0 0 0 0 0 0 0 0 -1.42+3.03i 0 0 0 -1.14+2.31i 2.56-5.34i];
Yline_G=xlsread('14bus_G.xlsx','Original', 'B2:O15');
Yline_B=xlsread('14bus_B.xlsx', 'B2:O15');
Yshunt_B=xlsread('14bus_shuntB.xlsx', 'B2:O15');
Yconn=xlsread('14bus_Connectivity.xlsx', 'B2:O15');
%disp(Yline_G);
%disp(Yline_B);
%disp(Yshunt_B);
%disp(Yconn);
B=eye(14);
for j=1:14
    Y(:,:,j)=B(:,j)*B(j,:)*Y_Bus;
    %disp(Y(:,:,j));
    Y_n_a(:,:,j)=0.5*[real(Y(:,:,j)+transpose(Y(:,:,j))) imag(transpose(Y(:,:,j))-Y(:,:,j));
        imag(Y(:,:,j)-transpose(Y(:,:,j))) real(Y(:,:,j)+transpose(Y(:,:,j)))];
    Y_n_r(:,:,j)=-0.5*[imag(Y(:,:,j)+transpose(Y(:,:,j))) real(Y(:,:,j)-transpose(Y(:,:,j)));
        real(transpose(Y(:,:,j))-Y(:,:,j)) imag(Y(:,:,j)+transpose(Y(:,:,j)))];
    %disp(Y_net_act(:,:,j));
    %disp(Y_net_react(:,:,j));
    M(:,:,j)=[B(:,j)*B(j,:) zeros(14,14);zeros(14,14) B(:,j)*B(j,:)];
end
c=0;
for l=1:14
    for m=1:14
        if Yconn(l,m)==1
            c=c+1;
            Y_line_r(:,:,c)=Yline_G(l,m)*B(:,l)*B(l,:)-Yline_G(l,m)*B(:,l)*B(m,:);
            Y_line_i(:,:,c)=(Yline_B(l,m)+Yshunt_B(l,m))*B(:,l)*B(l,:)-Yline_B(l,m)*B(:,l)*B(m,:);
            Y_line_net(:,:,c)=0.5*[Y_line_r(:,:,c)+transpose(Y_line_r(:,:,c)) transpose(Y_line_i(:,:,c))-Y_line_i(:,:,c);Y_line_i(:,:,c)-transpose(Y_line_i(:,:,c)) Y_line_r(:,:,c)+transpose(Y_line_r(:,:,c))];
           fprintf('From Bus To Bus Line Index\n');
            disp([l m c]);
            %disp(m);
            %disp(c);
        end
    end
end  
%disp(c);
P_k_Max=[3.324;1.4;0;0;0;0;0;0;0;0;0;0;0;0];
P_k_Min=[0;0;0;0;0;0;0;0;0;0;0;0;0;0];
Q_k_Max=[.10;.5;.4;0;0;.24;0;.24;0;0;0;0;0;0];
Q_k_Min=[-.2;-.4;0;0;0;-.06;0;-.06;0;0;0;0;0;0];
V_k_Max=[1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06;1.06];
V_k_Min=[0.94;0.94;0.94;0.94;0.94;0.94;0.94;0.94;0.94;0.94;0.94;0.94;0.94;0.94];
V_k_M=V_k_Max.*V_k_Max;
V_k_m=V_k_Min.*V_k_Min;
P_d_k=[0;.217;.942;.478;.076;.112;0;0;.295;.09;.035;.061;.135;.149];
Q_d_k=[0;.127;.190;-.039;.016;.075;0;0;-.024;.058;.018;.016;.058;.05];
c12=0.0430293;c11=20;c10=0;c22=.25;c21=20;c20=0;
reply=input('Enter the pu value of the real power line limit\n');
cvx_begin sdp
variables lambdak_m(14) lambdak_M(14)  lambda_k_m(14) lambda_k_M(14) mu_k_M(14) mu_k_m(14) lambda_lm(c);
variables r_1_1 r_1_2 r_2_1 r_2_2;
dual variable A;
expression S(28,28,15);
S(:,:,1)=zeros(28,28);
for count=1:14
    S(:,:,count+1)=S(:,:,count)+(lambdak_M(count)-lambdak_m(count))*Y_n_a(:,:,count)+(lambda_k_M(count)-lambda_k_m(count))*Y_n_r(:,:,count)+(mu_k_M(count)-mu_k_m(count))*M(:,:,count);
end
expression s(28,28,c+1);
s(:,:,1)=zeros(28,28);
for lcount=1:c
    s(:,:,lcount+1)=s(:,:,lcount)+(lambda_lm(lcount))*Y_line_net(:,:,lcount);
end
maximize transpose(lambdak_m)*P_k_Min-transpose(lambdak_M)*P_k_Max+transpose(lambdak_M-lambdak_m)*P_d_k+transpose(lambda_k_m)*Q_k_Min-transpose(lambda_k_M)*Q_k_Max+transpose(lambda_k_M-lambda_k_m)*Q_d_k+transpose(mu_k_m)*(V_k_m)-transpose(mu_k_M)*(V_k_M)+c10+c20-r_1_2-r_2_2+(c11+2*sqrt(c12)*r_1_1)*P_d_k(1)+(c21+2*sqrt(c22)*r_2_1)*P_d_k(2)-reply*(transpose(lambda_lm)*ones(c,1));
subject to
 lambda_lm>=0;
 lambdak_m>=0;
 lambdak_M>=0;
 lambda_k_m>=0;
 lambda_k_M>=0;
 mu_k_M>=0;
 mu_k_m>=0;
 [1 r_1_1;r_1_1 r_1_2]>=0;
 [1 r_2_1;r_2_1 r_2_2]>=0;
S(:,:,15)+s(:,:,c+1)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2)>=0:A;
cvx_end
fprintf('\nLagrange Multiplier of line real power flow limit constraints\n');
disp(lambda_lm);
fprintf('\nLagrange Multiplier of lower real power limit constraints\n');
disp(lambdak_m);
fprintf('\nLagrange Multiplier of upper real power limit constraints\n');
disp(lambdak_M);
fprintf('\nLagrange Multiplier of lower reactive power limit constraints\n');
disp(lambda_k_m);
fprintf('\nLagrange Multiplier of upper reactive power limit constraints\n');
disp(lambda_k_M);
fprintf('\nLagrange Multiplier of upper voltage limit constraints\n');
disp(mu_k_M);
fprintf('\nLagrange Multiplier of lower voltage limit constraints\n');
disp(mu_k_m);
fprintf('\nNumber of Iterations\n');
disp(cvx_slvitr);
fprintf('\nSolution Tolerance\n');
disp(cvx_slvtol);
fprintf('\nOptimal Dual Value\n');
disp(cvx_optval);
%disp([1 r_1_1;r_1_1 r_1_2]);
%disp([1 r_2_1;r_2_1 r_2_2]);
%disp(S(:,:,15)+s(:,:,c+1)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2));
[V,D]=eig(S(:,:,15)+s(:,:,c+1)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2));
fprintf('\nThese are the Eigenvectors and Eigenvalues of the Network Matrix\n');
disp([V,D]);
fprintf('\nThe following are the eigenvalues of the network matrix\n');
disp(eig(S(:,:,15)+s(:,:,c+1)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2)));
fprintf('\nThe following is the minimum eigenvalue of the network matrix\n');
disp(min(eig(S(:,:,15)+s(:,:,c+1)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2))));
fprintf('\nThe following are the eigenvalues of the dual matrix\n');
disp(eig(A));
for j=1:4
    fprintf('\nThese are the values corresponding to the %d th eigenvector\n',j);
for i=1:14
    if mu_k_M(i)>=0.0001
        Si_2=(1.06*V(15,j))/sqrt((((V(i,j))^2)+((V(i+14,j))^2))*(((V(1,j))^2)+((V(15,j))^2)));
        Si_1=-((Si_2)*V(1,j))/V(15,j);
        for f=1:14
            V_comp(f)=(Si_1)*V(f,j)-(Si_2)*V(f+14,j);
            V_comp(f+14)=(Si_2)*V(f,j)+(Si_1)*V(f+14,j);
            V_mod(f)=sqrt(((V_comp(f))^2)+((V_comp(f+14))^2));
        end
        W=((V_comp)')*V_comp;
for s=1:14
P(s)=trace(Y_n_a(:,:,s)*W);
Q(s)=trace(Y_n_r(:,:,s)*W);
end
P_Sol=c12*((P(1)+P_d_k(1))^2)+c11*(P(1)+P_d_k(1))+c22*((P(2)+P_d_k(2))^2)+c21*(P(2)+P_d_k(2));
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
fprintf('\nGeneration at Bus 1\n');
disp(P(1)+P_d_k(1));
fprintf('\nGeneration at Bus 2\n');
disp(P(2)+P_d_k(2));
for t=1:14
    if ((t~=1)&&(t~=2))
        LMPP(t)=lambdak_M(t)-lambdak_m(t);
    else if t==1
            LMPP(t)=lambdak_M(t)-lambdak_m(t)+2*(P(1)+P_d_k(1))*c12+c11;
        end
        if t==2
            LMPP(t)=lambdak_M(t)-lambdak_m(t)+2*(P(2)+P_d_k(2))*c22+c21;
        end
    end
    LMPQ(t)=lambda_k_M(t)-lambda_k_m(t);
end
fprintf('\nReal Power LMP\n');
disp(LMPP');
fprintf('\nReactive Power LMP\n');
disp(LMPQ');
lcount=0;
for l=1:14
    for m=1:14
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






        
    
