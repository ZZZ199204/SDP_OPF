function D=DiffInv(A)
global n
B=eye(2*n,2*n);
C=cell(2*n,2*n);
for i=1:2*n
    for j=1:2*n
        C{i,j}=(A\(B(:,i)))*((B(j,:))/A);
    end
end
D=cell2mat(C);