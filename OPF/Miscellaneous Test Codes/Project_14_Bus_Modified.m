clc
clear;
Y_Bus=[6.03-19.45i -5+15.26i 0 0 -1.03+4.23i 0 0 0 0 0 0 0 0 0;
    -5+15.26i 9.52-30.27i -1.14+4.78i -1.69+5.12i -1.7+5.19i 0 0 0 0 0 0 0 0 0;
    0 -1.14+4.78i 3.12-9.82i -1.99+5.07i 0 0 0 0 0 0 0 0 0 0;
    0 -1.69+5.12i -1.99+5.07i 10.5157-38.6542i -6.8410+21.5786i 0 -0.0023+4.8895i 0 -0.0003+1.8555i 0 0 0 0 0;
    -1.03+4.23i -1.7+5.19i 0 -6.84+21.58i 9.5698-35.5336i -0.0017+4.2574i 0 0 0 0 0 0 0 0;
    0 0 0 0 -0.0017+4.2574i 6.5815-17.3407i 0 0 0 0 -1.96+4.09i -1.53+3.18i -3.1+6.10i 0;
    0 0 0 -0.0023+4.8895i 0 0 0.0138-19.5490i -0.0032+5.6770i -0.0083+9.0901i 0 0 0 0 0;
    0 0 0 0 0 0 -0.0032+5.6770i 0.0032-5.6770i 0 0 0 0 0 0;
    0 0 0 -0.0003+1.8555i 0 0 -0.0083+9.0901i 0 5.3346-24.0925i -3.9+10.37i 0 0 0 -1.42+3.03i;
    0 0 0 0 0 0 0 0 -3.9+10.37i 5.78-14.77i -1.88+4.4i 0 0 0;
    0 0 0 0 0 -1.96+4.09i 0 0 0 -1.88+4.4i 3.84-8.5i 0 0 0;
    0 0 0 0 0 -1.53+3.18i 0 0 0 0 0 4.01-5.43i -2.49+2.25i 0;
    0 0 0 0 0 -3.1+6.10i 0 0 0 0 0 -2.49+2.25i 6.72-10.67i -1.14+2.31i;
    0 0 0 0 0 0 0 0 -1.42+3.03i 0 0 0 -1.14+2.31i 2.56-5.34i];
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
Y_act=[Y_n_a(:,:,1) Y_n_a(:,:,2) Y_n_a(:,:,3) Y_n_a(:,:,4) Y_n_a(:,:,5) Y_n_a(:,:,6) Y_n_a(:,:,7) Y_n_a(:,:,8) Y_n_a(:,:,9) Y_n_a(:,:,10) Y_n_a(:,:,11) Y_n_a(:,:,12) Y_n_a(:,:,13) Y_n_a(:,:,14)];
Y_react=[Y_n_r(:,:,1) Y_n_r(:,:,2) Y_n_r(:,:,3) Y_n_r(:,:,4) Y_n_r(:,:,5) Y_n_r(:,:,6) Y_n_r(:,:,7) Y_n_r(:,:,8) Y_n_r(:,:,9) Y_n_r(:,:,10) Y_n_r(:,:,11) Y_n_r(:,:,12) Y_n_r(:,:,13) Y_n_r(:,:,14)];
M_net=[M(:,:,1) M(:,:,2) M(:,:,3) M(:,:,4) M(:,:,5) M(:,:,6) M(:,:,7) M(:,:,8) M(:,:,9) M(:,:,10) M(:,:,11) M(:,:,12) M(:,:,13) M(:,:,14)];
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
cvx_begin
variables lambdak_m(14) lambdak_M(14)  lambda_k_m(14) lambda_k_M(14) mu_k_M(14) mu_k_m(14);
variables r_1_1 r_1_2 r_2_1 r_2_2;
dual variable A;
maximize transpose(lambdak_m)*P_k_Min-transpose(lambdak_M)*P_k_Max+transpose(lambdak_M-lambdak_m)*P_d_k+transpose(lambda_k_m)*Q_k_Min-transpose(lambda_k_M)*Q_k_Max+transpose(lambda_k_M-lambda_k_m)*Q_d_k+transpose(mu_k_m)*V_k_m-transpose(mu_k_M)*V_k_M+c10+c20-r_1_2-r_2_2+(c11+2*sqrt(c12)*r_1_1)*P_d_k(1)+(c21+2*sqrt(c22)*r_2_1)*P_d_k(2);
%minimize -(transpose(lambda_k_M)*Q_k_Max);
subject to
 lambdak_m>=0;
 lambdak_M>=0;
 lambda_k_m>=0;
 lambda_k_M>=0;
 mu_k_M>=0;
 mu_k_m>=0;
 [1 r_1_1;r_1_1 r_1_2]==semidefinite(2);
 [1 r_2_1;r_2_1 r_2_2]==semidefinite(2);
 M_net*kron((mu_k_M-mu_k_m),eye(28))+Y_react*kron((lambda_k_M-lambda_k_m),eye(28))+Y_act*kron((lambdak_M-lambdak_m),eye(28))+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2)==semidefinite(28):A;
cvx_end
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
%disp(S(:,:,15)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2));
[V,D]=eig(S(:,:,15)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2));
disp([V,D]);
for j=1:4
    fprintf('\nThese are the values corresponding to the %d th eigenvector\n',j);
for i=1:14
    if mu_k_M(i)>=0.001
        Si_2=(1.06*V(15,j))/sqrt((((V(i,j))^2)+((V(i+14,j))^2))*(((V(1,j))^2)+((V(15,j))^2)));
        Si_1=-((Si_2)*V(1,j))/V(15,j);
        for c=1:14
            V_comp(c)=(Si_1)*V(c,j)-(Si_2)*V(c+14,j);
            V_comp(c+14)=(Si_2)*V(c,j)+(Si_1)*V(c+14,j);
            V_mod(c)=sqrt(((V_comp(c))^2)+((V_comp(c+14))^2));
        end
        W=((V_comp)')*V_comp;
for c=1:14
P(c)=trace(Y_n_a(:,:,c)*W);
Q(c)=trace(Y_n_r(:,:,c)*W);
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
for c=1:14
    if ((c~=1)&&(c~=2))
        LMPP(c)=(lambdak_M(c)-lambdak_m(c))/100;
    else if c==1
            LMPP(c)=(lambdak_M(c)-lambdak_m(c)+2*(P(1)+P_d_k(1))*c12+c11)/100;
        end
        if c==2
            LMPP(c)=(lambdak_M(c)-lambdak_m(c)+2*(P(2)+P_d_k(2))*c22+c21)/100;
        end
    end
    LMPQ(c)=(lambda_k_M(c)-lambda_k_m(c))/100;
end
fprintf('\nReal Power LMP\n');
disp(LMPP');
fprintf('\nReactive Power LMP\n');
disp(LMPQ');
    end
end
end
fprintf('\nThe following are the eigenvalues of the network matrix\n');
disp(eig(S(:,:,15)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2)));
fprintf('\nThe following is the minimum eigenvalue of the network matrix\n');
disp(min(eig(S(:,:,15)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2))));
fprintf('\nThe following are the eigenvalues of the dual matrix\n');
disp(eig(A));




        
    
