% pump calculation
% D distillate mass flow - kg/s
% Q flow rate
% Pout & Pin - torr

% returns cost (dollars) of floating head heat exchanger

function Cp  = pump(D, Pin_torr)

Pout_torr = 760; %torr
Pout = Pout_torr .* 133.322368; % Pa
Pin = Pin_torr .* 133.322368; % Pa 

densH2O = 1; %kg/m^3
Q = D ./ densH2O; %m^3/s

% W = Q * (Pout - Pin); %

densH2O = 1000; %kg/m^3
gravity = 9.81; %m/s^2
H = (gravity.*densH2O).^(-1) .* (Pout - Pin); % meters

S = Q.*sqrt(H);

Cb = exp(9.7171 - 0.6017.*log(S) + 0.0519.*(log(S)).^2);

Ft = 8.90; % for 2+ stages, range 100-1500 gpm

Fm = 2.00; %for stainless steel

Cp_all = Cb .* Ft .* Fm;

Cp = sum(Cp_all);

%volFlow = Q .* 15852; %gpm


end