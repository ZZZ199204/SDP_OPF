clc
clear;
A=[12 0 0 0 56;35 67 0 0 43;0 0 0 56 0;0 0 0 67 0;98 0 89 0 0];
B=[2;9;3;5;11];
A1=sparse(A);
B1=sparse(B);
C=eye(5);
C1=speye(5);
M1=cell(5,1);
J=cell(5,1);
for i=1:5
K(:,:,i)=C(:,i)*C(i,:)*A;
K1=C1(:,i)*C1(i,:)*A1;
%K2=C(:,3)*C(3,:)*A1;
M(:,:,i)=[K(:,:,i)+transpose(K(:,:,i)) K(:,:,i)-transpose(K(:,:,i));-K(:,:,i)+transpose(K(:,:,i)) K(:,:,i)+transpose(K(:,:,i))];
M1{i}=[K1+transpose(K1) K1-transpose(K1);-K1+transpose(K1) K1+transpose(K1)];
J{i}=[K1+transpose(K1) zeros(5,5);zeros(5,5) K1+transpose(K1)];
%disp(J{i});
%K2(:,:,i)=full(K1);
M2(:,:,i)=full(M1{i});
%disp(K(:,:,i));
%disp(K2(:,:,i));
%disp(K2);
%disp(M1{i});
%disp(M(:,:,i));
%disp(M2(:,:,i));
end
F(:,:,1)=zeros(10,10);
for i=1:5
    F(:,:,i+1)=F(:,:,i)+M1{i};
end 
%disp(F(:,:,6));
F1=cell(6,1);
F1{1}=sparse(zeros(10,10));
for i=1:5
    F1{i+1}=F1{i}+M1{i};
end
F2=full(F1{6});
%disp(F2);
%disp(F1{6});
F=sparse(zeros(10,10))+M1{4};
%disp(F);
cvx_begin sdp
variables x(5);
minimize transpose(B)*x;
N1=cell(6,1);
N1{1}=zeros(10,10);
for i=1:5
    N1{i+1}=sparse(N1{i})+sparse(M1{i}*x(i));
end
subject to
N1{6}>=0;
x>=0;
cvx_end
disp(x);
[V,D]=eig(full(N1{6}));
disp([V,D]);
cvx_begin sdp
variables z(5);
minimize transpose(B)*z;
subject to
full(M1{1})*z(1)+full(M1{2})*z(2)+full(M1{3})*z(3)+full(M1{4})*z(4)+full(M1{5})*z(5)>=0;
z>=0;
cvx_end
disp(z);
[V1,D1]=eig(full(M1{1})*z(1)+full(M1{2})*z(2)+full(M1{3})*z(3)+full(M1{4})*z(4)+full(M1{5})*z(5));
disp([V1,D1]);
cvx_begin sdp
variables y(5);
minimize transpose(B)*y;
subject to
M(:,:,1)*y(1)+M(:,:,2)*y(2)+M(:,:,3)*y(3)+M(:,:,4)*y(4)+M(:,:,5)*y(5)>=0;
y>=0;
cvx_end
disp(y);
[V2,D2]=eig(M(:,:,1)*y(1)+M(:,:,2)*y(2)+M(:,:,3)*y(3)+M(:,:,4)*y(4)+M(:,:,5)*y(5));
disp([V2,D2]);