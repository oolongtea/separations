function [x y V L Q] = rachford_rice(F, z, K, hV, hL, hF)

% rachford_rice: Rachford-Rice procedure for Isothermal-Flash 
%                   when K-values are independent of composition
%       [x y V L Q] = rachford_rice(F, z, K, hV, hL, hF)
% input:
%  F = total feed flow
%  z = array of feed composition
%  K = array of K-values
% (hV = vapor enthalpy)  [optional]
% (hL = liquid enthalpy) [optional]
% (hF = feed enthalpy)   [optional]
%
% output:
%  x = array of bottom liquid composition
%  y = array of top vapor composition
%  V = vapor flow out of flash drum (top)
%  L = liquid flow out of flash drum (bottom)
%  Q = heat needed for nonadiabatic flash

if nargin<3,error('at least 3 input arguments required (F, z, K)'),end
if nargin>6,error('at most 8 input arguments allowed'),end
if nargin<4, hV = 0; hL = 0; hF = 0; end
if nargin<5, hL = 0; hF = 0; end
if nargin<6, hF = 0; end

comp = length(z); % number of components

% step 1 | thermal equilibrium %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T_liq = T_vap; % included just for your information
% step 2 | mechanical equilibrium%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P_liq = P_vap; % included just for your information

% step 3 | solve with Newton's method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxiter = 50;             % max iterations
tolerance = 0.001;        % relative error tolerance
psi = zeros(1, maxiter);  % psi = V/F, array preallocation
psi(1) = 0.5;             % initial guess for psi
relError = tolerance + 1; % value needed to start the while loop
k = 1;                    % iteration counter
f = zeros(1, maxiter);    % R-R equation, array preallocation
fp = zeros(1, maxiter);   % derivative of R-R equation, array preallocation

while ((relError > tolerance) || (k > 50))
  
  % calculate f{psi}
  for i = 1:comp 
    f(k) = f(k) + z(i)*(1 - K(i)) / (1 + psi(k)*(K(i) - 1)); 
  end
  
  % calculate f'{psi}
  for i = 1:comp 
    fp(k) = fp(k) + z(i)*(1 - K(i))^2 / (1 + psi(k)*(K(i) - 1))^2; 
  end
    
  % recursive solution for psi
  psi(k+1) = psi(k) - f(k)/fp(k);
  
  % calculate relative error of psi's
  relError = abs(psi(k+1) - psi(k))/psi(k);
  
  k = k + 1;
end

% step 4.1 | determine vapor flow %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V = F * psi(k-1); 
% step 4.2 | determine liquid flow %
L = F - V; 

% step 5 | find equilibrium compositions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psi = psi(k-1);
x = zeros(1, comp); y = zeros(1, comp);
for i = 1:comp
   x(i) = z(i) / (1 + psi*(K(i) -1));
   y(i) = x(i) * K(i);
end

% step 6 | find heat needed to flash %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q = hV*V + hL*L - hF*F;