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
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\14bus_Ybus.xlsx',Y_Bus,'Sheet1','A1:N14');
A=imag(Y_Bus);
xlswrite('\\austin.utexas.edu\disk\engrstu\ece\sc32589\My Documents\MATLAB\New folder\OPF\14bus_Ybus.xlsx',A,'Sheet2','A1:N14');