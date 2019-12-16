clc
clear;
A=[0;13;2;3;4;5;6;7;8;9];
A1=[1;2;3;4;5;67;6;8;34;29];
A2=[A';A1'];
fileID=fopen('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\test.txt','w');
fprintf(fileID,'%11s\t%11s\r\n','Test Values','T_Values');
fprintf(fileID,'%11d\t%11d\r\n',A2);
fclose(fileID);
%xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\test.xlsm',A1,'Sheet1','B1:B10');
%for i=1:10
    %for j=1:10
        %R(:,j,i)=randn(5,1);
    %end
%end
%xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\test.xlsm',R,'Sheet2','A1:Z100');
%A2=randn(10,2,8);
%A21=reshape(A2,10,2*8);
%B=1:1:10;
%A3=randn(10,3,9);
%A31=reshape(A3,10,3*9);
%plot(B,A21(:,1),'o',B,A31(:,1),'-');
%A=[1;2;3;4;5;6;2;8;12;54;9;5;17];
%[C,i]=max(A);
%disp(C);
%disp(i);
