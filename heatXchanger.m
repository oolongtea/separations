% heat exchanger calculation
%
% Th1 -> COOLER -> Th2
% Tc2 <- DESAL. <- Tc1
%
% Q inout in MW

% returns cost (dollars) of floating head heat exchanger

function Cp = heatXchanger(Q_MW, Th2, Th1, Tc2, Tc1, U, a, b)

% MW to Btu/hr
Q = Q_MW * 3.6e9 * 0.00094781712; % Btu/hr

% pressure out of turbine
P = 14.6; %psig

% deltaTh - cooler input/output
deltaT1 = Th1 - Tc2;
% deltaTc - cooler input/output
deltaT2 = Th2 - Tc1;

logMeanTemp = (deltaT1 - deltaT2)/log(deltaT1/deltaT2);

% overall heat transfer coefficient estimate
% between Air/N2 and water/brine
%U = 30; % Btu / (degF-sqft-hr)

A = Q/U/logMeanTemp;

% floating head
Cb_float = exp(11.667 - 0.8709*log(A) + 0.09005*(log(A))^2);
% fixed head
Cb_fixed = exp(11.0545 - 0.9228*log(A) + 0.09861*(log(A))^2);
% U-tube
Cb_Utube = exp(11.147 - 0.9186*log(A) + 0.09790*(log(A))^2);
% kettle vaporizer
Cb_kettle = exp(11.967 -     0.8709*log(A) + 0.09005*(log(A))^2);
Cb = [Cb_float Cb_fixed Cb_Utube Cb_kettle];

% stainless steel
% a = 2.70;
% b = 0.07;
Fm = a + (A./100).^b;

% for tube length of 20 ft
Fl = 1.0;

%
Fp = 0.9803 + 0.018*(P/100) + 0.0017*(P/100)^2;

Cp_all = Fp.*Fm.*Fl.*Cb;
Cp = Cp_all(2);

end
