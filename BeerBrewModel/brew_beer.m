clear all; close all; clc; format long;

G = [10 20 30 40 50 60 70 80]; % glucose concentration (mol/m^3)
M = [60 80 100 120 140 160 180 200]; % maltose concentration (mol/m^3)
N = [5 10 15 20 25 30 35 40]; % maltotriose concentration (mol/m^3)

% Specific rates of sugar uptake
%miu1 as a function of G hr^-1
miu1 = [0.02989 0.03005 0.03010 0.03012 0.03014 0.03014 0.03016 0.03016];
%miu2 as a function of M and G hr^-1
miu2 = [0.00502 0.00620 0.00711 0.00799 0.00860 0.00903 0.00949 0.01001];
%miu3 as a function of N,M and G hr^-1
miu3 = [0.00025 0.00053 0.00084 0.00104 0.00132 0.00146 0.00171 0.00188];

% PART A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LINEAR FITS

% miu 1
n_miu1 = length(miu1);
St_miu1 = var(miu1)*(n_miu1-1);

Z_mat_miu1 = [ones(n_miu1,1) 1./G'];
miu1_vec = 1./miu1';
a_1_miu1 = (Z_mat_miu1'*Z_mat_miu1)\[Z_mat_miu1'*miu1_vec];
% a(1) = 1/Vg
% a(2) = Kg/Vg
inv_miu1 = Z_mat_miu1*a_1_miu1;
miu1fit_1 = 1./inv_miu1';
Sr_miu1_lin = sum((miu1-miu1fit_1).^2);
r2_1_miu1 = (St_miu1-Sr_miu1_lin)/St_miu1 % coefficient of determination

% output guess parameters for nonlinear fit
Vg_lin = 1/a_1_miu1(1);
Kg_lin = a_1_miu1(2)*Vg_lin;



% miu 2
n_miu2 = length(miu2);
St_miu2 = var(miu2)*(n_miu2-1);

Z_mat_miu2 = [ones(n_miu2,1) G' (1./M)' (G./M)'];
miu2_vec = 1./miu2';
a_1_miu2 = (Z_mat_miu2'*Z_mat_miu2)\[Z_mat_miu2'*miu2_vec];
% a(1) = 1/Vm
% a(2) = 1/(Kgp*Vm)
% a(3) = Km/Vm
% a(4) = Km/(Kg'*Vm)  (overspecifies system and <0; --> disregard)
inv_miu2 = Z_mat_miu2*a_1_miu2;
miu2fit_1 = 1./inv_miu2';
Sr_miu2_lin = sum((miu2-miu2fit_1).^2);
r2_1_miu2 = (St_miu2-Sr_miu2_lin)/St_miu2 % coefficient of determination

% output guess parameters for nonlinear fit
Vm_lin = 1/a_1_miu2(1);
Kgp_lin = 1/(a_1_miu2(2)*Vm_lin);
Km_lin = Vm_lin*a_1_miu2(3);

% miu 3
n_miu3 = length(miu3);
St_miu3 = var(miu3)*(n_miu3-1);

Z_mat_miu3 = [(1 + G./Kgp_lin)' (1./N + G./(N.*Kgp_lin))' (M + G.*M./Kgp_lin)'...
                (G.*M./(N.*Kgp_lin) + M./N)' ];
miu3_vec = 1./miu3';
a_1_miu3 = (Z_mat_miu3'*Z_mat_miu3)\[Z_mat_miu3'*miu3_vec];
% a(1) = 1/Vn
% a(3) = 1/(Kmp*Vn)
% a(2) = Kn/Vn
% a(4) = Kn/(Km'*Vn) (overspecifies system and <0; --> disregard)
inv_miu3 = Z_mat_miu3*a_1_miu3;
miu3fit_1 = 1./inv_miu3';
Sr_miu3_lin = sum((miu3-miu3fit_1).^2);
r2_1_miu3 = (St_miu3-Sr_miu3_lin)/St_miu3 % coefficient of determination

% output guess parameters for nonlinear fit
Vn_lin = 1/a_1_miu3(1)
Kmp_lin = 1/(Vn_lin*a_1_miu3(3))
Kn_lin = Vn_lin*a_1_miu3(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NON-LINEAR FITS

options = optimset('MaxFunEvals', 2000, 'TolFun', 10^-7,'TolX',10^-5, 'MaxIter', 2000);

% miu 1
a0_miu1 = [Vg_lin Kg_lin];
[a_3_miu1 Sr_miu1_nonlin] = fminsearch(@muone, a0_miu1',...
    options, G,miu1);

Vg = a_3_miu1(1)
Kg = a_3_miu1(2)


% miu 2 
a0_miu2 = [Vm_lin Kgp_lin Km_lin]
[a_3_miu2 Sr_miu2_nonlin] = fminsearch(@mutwo, a0_miu2',...
    options, G,M,miu2);

Vm = a_3_miu2(1)
Kgp = a_3_miu2(2)
Km = a_3_miu2(3)


% miu 3
a0_miu3 = [Vn_lin Kmp_lin Kn_lin]
[a_3_miu3 Sr_miu3_nonlin] = fminsearch(@muthree, a0_miu3',...
    options, G,M,N,miu3,Kgp);

Vn = a_3_miu3(1)
Kmp = a_3_miu3(2)
Kn = a_3_miu3(3)


% Plotting
G_plot = 10:80;
% M_plot = 10:80;
% N_plot = 10:80;

figure; p = plot(G,miu1,'o', G_plot, Vg_lin.*G_plot./(Kg_lin+G_plot),'--d',...
    G_plot, Vg.*G_plot./(Kg+G_plot), '--+');
set(p(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel('G (mol/m^3)')
ylabel('miu1 (hr^{-1})')
title('Specific rate of glucose uptake vs. glucose concentration')
legend('data', 'linear fit', 'nonlinear fit', 'location', 'best')


figure; pp = plot(G,miu2,'o',...
    G,Vm_lin.*M.*Kgp_lin./(Km_lin + M)./(Kgp_lin + G),'--d',...
    G,Vm.*M.*Kgp./(Km+M)./(Kgp + G), '--+');
set(pp(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel('G (mol/m^3)')
ylabel('miu2 (hr^{-1})')
title('Specific rate of sugar uptake vs. glucose concentration')
legend('data', 'linear fit', 'nonlinear fit', 'location', 'best')


figure; ppp = plot(G,miu3,'o',...
    G,Vn_lin.*N.*Kgp_lin.*Kmp_lin./(Kn_lin+N)./(Kgp_lin+G)./(Kmp_lin+M),'--d',...
    G,Vn.*N.*Kgp.*Kmp./(Kn+N)./(Kgp+G)./(Kmp+M), '--+');
set(ppp(1),'MarkerEdgeColor','k','MarkerFaceColor','k')
xlabel('G (mol/m^3)')
ylabel('miu3 (hr^{-1})')
title('Specific rate of sugar uptake vs. glucose concentration')
legend('data', 'linear fit', 'nonlinear fit', 'location', 'best')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART B


%Fructose #carbons = 6, called C(1) = 6, Re(1), Rx(1),delHf(1)
%Stachyose #carbons = 24, called C(2) = 24, Re(2), Rx(2), delHf(2)
%Fructosylnystose #carbons = 28, called C(3) = 28, Re(3), Rx(3),delHf(3)

%number of carbons
C = [6 24 28];

%stoichiometric yield of ethanol per mole of sugar reacted
Re = [1.92 7.68 9.60];

%stoichiometric yield of biomass per mole of sugar reacted
Rx = [0.13 0.54 0.67];

%heats of reaction
delHf = [-91.20 -496.36 -631.41];

%plot
figure;
%number of carbons vs. stoichiometric yield of ethanol
subplot(3,1,1);
plot(C,Re,'-o');
title('Effect of Number of Carbons on Properties of Sugars')
grid;
xlabel('Number of carbons');
ylabel('stoichiometric yield of ethanol');
%number of carbons vs. stoichiometric yield of biomass
subplot(3,1,2);
plot(C,Rx,'-o');
grid;
xlabel('Number of carbons');
ylabel('stoichiometric yield of biomass');
%number of carbons vs. heats of reaction
subplot(3,1,3);
plot(C,delHf,'-o');
grid;
xlabel('Number of carbons');
ylabel('heats of reaction');

%based on the plots, between 6 and 24 carbons, there is a direct positive
%linear correlation between number of carbons and stoichiometric yield of 
%biomass and ethanol; and there is a inverse negative linear correlation 
%between number of carbons and heats of reaction. However, the fit deviates 
%from linearity for number of carbons greater than 24. 
%We chose the best fit to be linear because the data has a
%strong linear correlation. 

%linear model, output(1) = slope, ouput(2)= y-intercept 
%number of carbons
CPrime = [6 24 28];

%stoichiometric yield of ethanol per mole of sugar reacted
RePrime = [1.92 7.68 9.60];

%stoichiometric yield of biomass per mole of sugar reacted
RxPrime = [0.13 0.54 0.67];

%heats of reaction
delHfPrime = [-91.20 -496.36 -631.41];

%Perform linear regression on ethanol
ethanol = polyfit(CPrime,RePrime,1)
ethanol_fit = polyval(ethanol,CPrime);

%Perform linear regression on biomass
biomass = polyfit(CPrime,RxPrime,1)
biomass_fit = polyval(biomass,CPrime);

%Perform linear regression on heats of reaction
heat = polyfit(CPrime,delHfPrime,1)
heat_fit = polyval(heat,CPrime);

%plot
figure;
%number of carbons vs. stoichiometric yield of ethanol, data and fit
subplot(3,1,1);
plot(CPrime,RePrime,'o',CPrime,ethanol_fit);
title('Effect of Number of Carbons on Properties of Sugars')
grid;
xlabel('Number of carbons');
ylabel('stoichiometric yield of ethanol');
legend ('data','fit')
%number of carbons vs. stoichiometric yield of biomass, data and fit
subplot(3,1,2);
plot(CPrime,RxPrime,'o',CPrime,biomass_fit);
grid;
xlabel('Number of carbons');
ylabel('stoichiometric yield of biomass');
legend ('data','fit')
%number of carbons vs. heats of reaction, data and fit
subplot(3,1,3);
plot(CPrime,delHfPrime,'o',CPrime,heat_fit);
grid;
xlabel('Number of carbons');
ylabel('heats of reaction');
legend ('data','fit')

%Test our fits
% Calculate standard error of the estimate and compare with
% standard deviation of Ethanol
SyEthanol = std(RePrime);
nEthanol = length(RePrime);
SrEthanol = 0.0;
for i=1:nEthanol
SrEthanol = SrEthanol+(RePrime(i)-ethanol_fit(i))^2;
end
SyxEthanol = sqrt(SrEthanol/(nEthanol-2));

%Calculate coefficient of determination
StEthanol = var(RePrime)*(nEthanol-1);
r2Ethanol = (StEthanol-SrEthanol)/StEthanol

% Calculate standard error of the estimate and compare with
%   standard deviation of Biomass
SyBiomass = std(RxPrime);
nBiomass = length(RxPrime);
SrBiomass = 0.0;
for i=1:nBiomass
SrBiomass = SrBiomass+(RxPrime(i)-biomass_fit(i))^2;
end
SyxBiomass = sqrt(SrBiomass/(nBiomass-2));

%Calculate coefficient of determination
StBiomass = var(RxPrime)*(nBiomass-1);
r2Biomass = (StBiomass-SrBiomass)/StBiomass

% Calculate standard error of the estimate and compare with
%   standard deviation of Heat of reaction
SyHeat = std(delHfPrime);
nHeat = length(delHfPrime);
SrHeat = 0.0;
for i=1:nEthanol
SrHeat = SrHeat+(delHfPrime(i)-heat_fit(i))^2;
end
SyxHeat = sqrt(SrHeat/(nHeat-2));
%Calculate coefficient of determination
StHeat = var(delHfPrime)*(nHeat-1);
r2Heat = (StHeat-SrHeat)/StHeat

%Glucose; number of carbons = 6
%Calculate stoichiometric yield of ethanol per mole of Glucose,ReGlucose
ReGlucose = ethanol(1)*6+ethanol(2)
%Calculate stoichiometric yield of biomass per mole of Glucose,RxGlucose
RxGlucose = biomass(1)*6+biomass(2)
%Calculate heat of reaction of Glucose,delHfGlucose
delHfGlucose = heat(1)*6+heat(2)

%Maltose; number of carbons = 12
%Calculate stoichiometric yield of ethanol per mole of Maltose,ReMaltose
ReMaltose = ethanol(1)*12+ethanol(2)
%Calculate stoichiometric yield of biomass per mole of Maltose,RxMaltose
RxMaltose = biomass(1)*12+biomass(2)
%Calculate heat of reaction of Maltose,delHfMaltose
delHfMaltose = heat(1)*12+heat(2)

%Maltotriose; number of carbons = 18
%Calculate stoichiometric yield of ethanol per mole of Maltotriose,Maltotriose
ReMaltotriose = ethanol(1)*18+ethanol(2)
%Calculate stoichiometric yield of biomass per mole of Maltotriose,RxMaltotriose
RxMaltotriose = biomass(1)*18+biomass(2)
%Calculate heat of reaction of Maltose,delHfMaltotriose
delHfMaltotriose = heat(1)*18+heat(2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PART C

% creating constant arrays

c = [Vg,Kg,Vm,Km,Kgp,Vn,Kn,Kmp];
delHf_c = [delHfGlucose delHfMaltose delHfMaltotriose];
Rx_c = [RxGlucose RxMaltose RxMaltotriose];

c = [c delHf_c Rx_c];

c2 = [0.030200040063171 0.102980133548209 0.031895563023332 ...
      3.067526184942582e+02 2.998597665119866e+02 ...
      5.041431090255855e+11 6.913941597894742e+16 ...
      8.444321344291284e+15 delHf_c Rx_c];

G0 = 105; %(mol/m^3)
M0 = 150; %(mol/m^3)
N0 = 80; %(mol/m^3)
T0 = 7; %(deg C)

tspan = [0 600];
x0 = [G0 M0 N0 T0];


A = 0.1; %m^2

thisthing = 1;

while thisthing ~=0
    A = A+0.01;
    [t_out,x_out] = ode45(@dTdt_version3,tspan,x0,odeset('RelTol',1E-7),c2,A);%,delHf_c,Rx_c,Vg,Kg,Vm,Km,Kgp,Vn,Kn,Kmp);
    T = x_out(:,4);
    t_index = find(t_out>50,1);
    T = T(t_index:end);
    badvalues = find((T < 15) | (T > 17));
    thisthing2 = sum(badvalues);
    if thisthing2 == 0
        thisthing = 0;
    end
    disp(A)
end

figure('name','Temperature Profile');
plot(t_out,x_out(:,4));
axis([t_out(1) t_out(end) 15 17]);
title('Reactor Temperature Profile');
ylabel('Temperature (°C)'), xlabel('Time (hr)');

figure('name','Glucose');
plot(t_out,x_out(:,1));
title('Glucose Concentration Profile');
ylabel('Concentration (mol/m^3)'), xlabel('Time (hr)');

figure('name','Maltose');
plot(t_out,x_out(:,2));
title('Maltose Concentration Profile');
ylabel('Concentration (mol/m^3)'), xlabel('Time (hr)');

figure('name','Maltotriose');
plot(t_out,x_out(:,3));
title('Maltotriose Concentration Profile');
ylabel('Concentration (mol/m^3)'), xlabel('Time (hr)');

figure('name','Biomass');
biomass = 75 + RxGlucose.*(G0-x_out(:,1)) + RxMaltose.*(M0-x_out(:,2)) + ...
               RxMaltotriose.*(N0-x_out(:,3));
plot(t_out,biomass);
title('Biomass Concentration Profile');
ylabel('Concentration (mol/m^3)'), xlabel('Time (hr)');

figure('name','Ethanol');
ethanol = 0 + ReGlucose.*(G0-x_out(:,1)) + ReMaltose.*(M0-x_out(:,2)) + ...
               ReMaltotriose.*(N0-x_out(:,3));
plot(t_out, 100.*(ethanol.*46.06844.*0.001.*(1/789).*0.1)./0.1);
title('Ethanol Percentage');
ylabel('Percent (%)'), xlabel('Time (hr)');

figure('name','All Concentrations');
plot(t_out,x_out(:,1:3),t_out,biomass);
title('Concentrations of Components');
ylabel('Concentration (mol/m^3)'), xlabel('Time (hr)');
legend('Glucose','Maltose','Maltotriose','Biomass');
