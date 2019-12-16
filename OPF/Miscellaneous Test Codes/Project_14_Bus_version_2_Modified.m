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
cvx_begin sdp
variables lambdak_m(14) lambdak_M(14)  lambda_k_m(14) lambda_k_M(14) mu_k_M(14) mu_k_m(14);
variables r_1_1 r_1_2 r_2_1 r_2_2;
dual variable A;
expression S(28,28,15);
S(:,:,1)=zeros(28,28);
for count=1:14
    S(:,:,count+1)=S(:,:,count)+(lambdak_M(count)-lambdak_m(count))*Y_n_a(:,:,count)+(lambda_k_M(count)-lambda_k_m(count))*Y_n_r(:,:,count)+(mu_k_M(count)-mu_k_m(count))*M(:,:,count);
end
maximize transpose(lambdak_m)*P_k_Min-transpose(lambdak_M)*P_k_Max+transpose(lambdak_M-lambdak_m)*P_d_k+transpose(lambda_k_m)*Q_k_Min-transpose(lambda_k_M)*Q_k_Max+transpose(lambda_k_M-lambda_k_m)*Q_d_k+transpose(mu_k_m)*(V_k_m.*V_k_m)-transpose(mu_k_M)*(V_k_M.*V_k_M)+c10+c20-r_1_2-r_2_2+(c11+2*sqrt(c12)*r_1_1)*P_d_k(1)+(c21+2*sqrt(c22)*r_2_1)*P_d_k(2);
subject to
 lambdak_m>=0;
 lambdak_M>=0;
 lambda_k_m>=0;
 lambda_k_M>=0;
 mu_k_M>=0;
 mu_k_m>=0;
 [1 r_1_1;r_1_1 r_1_2]>=0;
 [1 r_2_1;r_2_1 r_2_2]>=0;
S(:,:,15)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2)>=0:A;
cvx_end
disp(lambdak_m);
disp(lambdak_M);
disp(lambda_k_m);
disp(lambda_k_M);
disp(mu_k_M);
disp(mu_k_m);
disp(cvx_slvitr);
disp(cvx_slvtol);
%disp([1 r_1_1;r_1_1 r_1_2]);
%disp([1 r_2_1;r_2_1 r_2_2]);
%disp(S(:,:,15)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2));
disp(eig(S(:,:,15)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2)));
disp(min(eig(S(:,:,15)+(c11+2*sqrt(c12)*r_1_1)*Y_n_a(:,:,1)+(c21+2*sqrt(c22)*r_2_1)*Y_n_a(:,:,2))));
disp(rank(A));




        
    
