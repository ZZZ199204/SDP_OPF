%%fix the tap issue?????????????????????
%%%fix the capacitance at 9 (19)
%%%% mpc.gen(3) not important...

clear
clc


%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc_bus = [
	1	3	0	0	0	0	1	1.06	0	0	1	1.06	0.94;
	2	2	21.7	12.7	0	0	1	1.045	-4.98	0	1	1.06	0.94;
	3	2	94.2	19	0	0	1	1.01	-12.72	0	1	1.06	0.94;
	4	1	47.8	-3.9	0	0	1	1.019	-10.33	0	1	1.06	0.94;
	5	1	7.6	1.6	0	0	1	1.02	-8.78	0	1	1.06	0.94;
	6	2	11.2	7.5	0	0	1	1.07	-14.22	0	1	1.06	0.94;
	7	1	0	0	0	0	1	1.062	-13.37	0	1	1.06	0.94;
	8	2	0	0	0	0	1	1.09	-13.36	0	1	1.06	0.94;
	9	1	29.5	16.6	0	19	1	1.056	-14.94	0	1	1.06	0.94;
	10	1	9	5.8	0	0	1	1.051	-15.1	0	1	1.06	0.94;
	11	1	3.5	1.8	0	0	1	1.057	-14.79	0	1	1.06	0.94;
	12	1	6.1	1.6	0	0	1	1.055	-15.07	0	1	1.06	0.94;
	13	1	13.5	5.8	0	0	1	1.05	-15.16	0	1	1.06	0.94;
	14	1	14.9	5	0	0	1	1.036	-16.04	0	1	1.06	0.94;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc_gen = [
	1	232.4	-16.9	10	0	1.06	100	1	332.4	0	0	0	0	0	0	0	0	0	0	0	0;
	2	40	42.4	50	-40	1.045	100	1	140	0	0	0	0	0	0	0	0	0	0	0	0;
	3	0	23.4	40	0	1.01	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	6	0	12.2	24	-6	1.07	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
	8	0	17.4	24	-6	1.09	100	1	100	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc_branch = [
	1	2	0.01938	0.05917	0.0528	9900	0	0	0	0	1	-360	360;
	1	5	0.05403	0.22304	0.0492	9900	0	0	0	0	1	-360	360;
	2	3	0.04699	0.19797	0.0438	9900	0	0	0	0	1	-360	360;
	2	4	0.05811	0.17632	0.034	9900	0	0	0	0	1	-360	360;
	2	5	0.05695	0.17388	0.0346	9900	0	0	0	0	1	-360	360;
	3	4	0.06701	0.17103	0.0128	9900	0	0	0	0	1	-360	360;
	4	5	0.01335	0.04211	0	9900	0	0	0	0	1	-360	360;
	4	7	0	0.20912	0	9900	0	0	0.978	0	1	-360	360;
	4	9	0	0.55618	0	9900	0	0	0.969	0	1	-360	360;
	5	6	0	0.25202	0	9900	0	0	0.932	0	1	-360	360;
	6	11	0.09498	0.1989	0	9900	0	0	0	0	1	-360	360;
	6	12	0.12291	0.25581	0	9900	0	0	0	0	1	-360	360;
	6	13	0.06615	0.13027	0	9900	0	0	0	0	1	-360	360;
	7	8	0	0.17615	0	9900	0	0	0	0	1	-360	360;
	7	9	0	0.11001	0	9900	0	0	0	0	1	-360	360;
	9	10	0.03181	0.0845	0	9900	0	0	0	0	1	-360	360;
	9	14	0.12711	0.27038	0	9900	0	0	0	0	1	-360	360;
	10	11	0.08205	0.19207	0	9900	0	0	0	0	1	-360	360;
	12	13	0.22092	0.19988	0	9900	0	0	0	0	1	-360	360;
	13	14	0.17093	0.34802	0	9900	0	0	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc_gencost = [
	2	0	0	3	0.0430293	20	0;
	2	0	0	3	0.25	20	0;
	2	0	0	3	0.01	40	0;
	2	0	0	3	0.01	40	0;
	2	0	0	3	0.01	40	0;
];


% mpc_gencost = [
% 	2	0	0	3	0	1	0;
% 	2	0	0	3	0	1	0;
% 	2	0	0	3	0	1	0;
% 	2	0	0	3	0	1	0;
% 	2	0	0	3	0	1	0;
% ];


n=14;

Y=zeros(n,n);
for k=1:size(mpc_branch,1),
    Y(mpc_branch(k,1),mpc_branch(k,2))=-1/(mpc_branch(k,3)+mpc_branch(k,4)*i);
    Y(mpc_branch(k,2),mpc_branch(k,1))=-1/(mpc_branch(k,3)+mpc_branch(k,4)*i);
    C(mpc_branch(k,1),mpc_branch(k,2))=mpc_branch(k,5)*i/2;
    C(mpc_branch(k,2),mpc_branch(k,1))=mpc_branch(k,5)*i/2;
end

for k=1:n,
    Y(k,k)=-sum(Y(k,:))+sum(C(k,:));
end

YY=cell(n,1);
Y_aux=cell(n,1);
Y_real=cell(n,1);
Y_imag=cell(n,1);
M=cell(n,1);

for k=1:size(mpc_branch,1),
    if mpc_branch(k,9)~=0,
        x2=mpc_branch(k,1);
        x1=mpc_branch(k,2);
        a_aux=mpc_branch(k,9);
        Y(x1,x2)=Y(x1,x2)+1/(mpc_branch(k,4)*i)-1/(mpc_branch(k,4)*i)/a_aux;
        Y(x2,x2)=Y(x2,x2)-1/(mpc_branch(k,4)*i)+1/(mpc_branch(k,4)*i)/a_aux^2;
        Y(x2,x1)=Y(x2,x1)+1/(mpc_branch(k,4)*i)-1/(mpc_branch(k,4)*i)/a_aux;
    end
end


for k=1:n;
    T=eye(n);
    Y_aux{k}=T(:,k)*T(k,:)*Y;
    Y_real{k}=[real(Y_aux{k}) -imag(Y_aux{k});imag(Y_aux{k}) real(Y_aux{k})];
    Y_imag{k}=[imag(Y_aux{k}) real(Y_aux{k}); -real(Y_aux{k}) imag(Y_aux{k})];
    Y_real{k}=(Y_real{k}+Y_real{k}')/2;
    Y_imag{k}=(Y_imag{k}+Y_imag{k}')/2;
    
end


for k=1:n;
    T=eye(n);
    M{k}=diag([T(k,:) T(k,:)]);    
end


mpc_bus(:,4)=mpc_bus(:,4)-mpc_bus(:,6);

 v_min=sdpvar(n,1);
 v_max=sdpvar(n,1);
 r_min=sdpvar(n,1);
 r_max=sdpvar(n,1);
 a_min=sdpvar(n,1);
 a_max=sdpvar(n,1);
 for k = 1:size(mpc_gen,1),
     H{k} = sdpvar(2,2);
 end
 
 R_obj=0;
 A_obj=0;
 V_obj=0;
 L_obj=0;
 
 mpc_bus_m=mpc_bus;
 mpc_bus_m=[mpc_bus_m(:,1:3) mpc_bus_m(:,3:4) mpc_bus_m(:,4:size(mpc_bus_m,2))]; 
 
 for k=1:size(mpc_gen,1),
     t=mpc_gen(k,1);
      mpc_bus_m(t,4)=mpc_bus_m(t,4)-mpc_gen(k,9);
      mpc_bus_m(t,5)=mpc_bus_m(t,5)-mpc_gen(k,4);
      mpc_bus_m(t,6)=mpc_bus_m(t,6)-mpc_gen(k,5);

 end
 
 mpc_bus_m(:,15)=10*mpc_bus_m(:,15);
 mpc_bus_m(:,14)=10*mpc_bus_m(:,14);

 for k=1:n,
     V_obj=V_obj+v_min(k)*mpc_bus_m(k,15)^2-v_max(k)*mpc_bus_m(k,14)^2;
     A_obj=A_obj+a_min(k)*mpc_bus_m(k,4)-a_max(k)*mpc_bus_m(k,3);
     R_obj=R_obj+r_min(k)*mpc_bus_m(k,5)-r_max(k)*mpc_bus_m(k,6);
     
 end
 

V_cons=zeros(2*n,2*n);
R_cons=zeros(2*n,2*n);
A_cons=zeros(2*n,2*n);
L_cons=zeros(2*n,2*n);

F = set([]) 

 for k=1:n,
     V_cons=V_cons+(-v_min(k)+v_max(k))*M{k};
     R_cons=R_cons-(r_min(k)-r_max(k))*Y_imag{k};
     A_cons=A_cons+(a_min(k)-a_max(k))*Y_real{k};
 end
 
for k=1:size(mpc_gen,1);
    T_aux=H{k};
    L_cons=L_cons+(mpc_gencost(k,6)+2*sqrt(mpc_gencost(k,5))*T_aux(1,2))*Y_real{mpc_gen(k,1)};
    L_obj=L_obj-T_aux(2,2)+(mpc_gencost(k,6)+2*sqrt(mpc_gencost(k,5))*T_aux(1,2))*mpc_bus(mpc_gen(k,1),3);
 end
    
 F = F + set(L_cons+ A_cons+R_cons+V_cons>=0); 
 F = F + set(v_min>=0); 
 F = F + set(v_max>=0); 
 F = F + set(r_min>=0); 
 F = F + set(r_max>=0); 
 F = F + set(a_min>=0); 
 F = F + set(a_max>=0); 
  
  for k = 1:size(mpc_gen,1),
     F = F + set(H{k}>=0); 
     T_aux=H{k};
     F = F + set(T_aux(1,1)==1);
 end

 
g=solvesdp(F,-(A_obj+R_obj+V_obj+L_obj),sdpsettings('solver','sedumi',...,
    'sedumi.stepdif',1,'sedumi.alg',1,'sedumi.numtol',1e-10,'sedumi.eps',1e-10));


cons=full( double(L_cons)+ double(A_cons)+double(R_cons)+double(V_cons));

 eig(cons)

[d,dd]=eig(cons);
d=d(:,1);
v=( d(1:n)+d((n+1):2*n)*i)/(d(8)+d(n+8)*i)*mpc_bus_m(8,14);
I=Y*v;
Power=v.*transpose(I');
Power=[Power Power];
Power(:,1)=real(Power(:,1))+mpc_bus(:,3);
Power(:,2)=imag(Power(:,2))+mpc_bus(:,4);
[[1:n]' Power]

%[double(r_min-r_max) double(a_min-a_max) double(v_min-v_max)]

%double((A_obj+R_obj+V_obj+L_obj))
% 
% for k=1:size(mpc_gen,1);
%     T_aux=H{k};
%     double(mpc_gencost(k,6)+2*sqrt(mpc_gencost(k,5))*T_aux(1,2))
%     
%  end
% 