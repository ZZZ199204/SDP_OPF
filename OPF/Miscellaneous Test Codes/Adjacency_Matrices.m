clc
clear;
A14=xlsread('14bus_Connectivity.xlsx', 'B2:O15');
A30=xlsread('30bus_Connectivity.xlsx', 'B2:AE31');
A57=xlsread('57bus_Connectivity.xlsx', 'B2:BF58');
A118=xlsread('118bus_Connectivity.xlsx', 'B2:DO119');
A300=xlsread('300bus_Connectivity.xlsx', 'B2:KO301');
Bus_14=A14+eye(14,14);
Bus_30=A30+eye(30,30);
Bus_57=A57+eye(57,57);
Bus_118=A118+eye(118,118);
Bus_300=A300+eye(300,300);
disp(Bus_14);
disp(Bus_30);
disp(Bus_57);
disp(Bus_118);
disp(Bus_300);
