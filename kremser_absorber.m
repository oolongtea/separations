% Multi-Component Vapor-Liquid Casacdes %
close all; clear all; clc; format shortG;

% absorber %

%%%  component name       i %%%
s = ['C1 ';             % 1
     'C2 ';             % 2
     'C3 ';             % 3
     'nC4';             % 4
     'nC5';             % 5
     'Oil';             % 6
     'Tot';             % total
    ];

% average entering temperature
T_in = (90 + 105)/2; % degF

% number of stages
N = 6;

% K-values at 97.5 degF & 400 psia
K = [6.65 1.64 0.584 0.195 0.0713 0.0001];

% absorbent liquid (oil in this example)
% L_in = [total 1 2 3 4 5 6] % input absorber flow
L_in = [0 0 0 0.05 0.78 164.17]; % lbmol/h
L_in = cat(2, sum(L_in), L_in);
% feed gas
% V_in = [1 2 3 4 5 6 total] % feed gas
V_in = [160 370 240 25 5 0]; % lbmol/h
V_in(N+1) = sum(V_in);

% absorption factors
A = L_in(1) ./ K ./ V_in(N+1); 
Ae = A; % effective absorption factor for Kremser eqn
% stripping factors
S = 1 ./ A;
Se = S; % effective stripping factor for Kremser eqn

% Kremser equations
phiA = (Ae - 1) ./ (Ae.^(N+1) -1); phiA(end) = 0;
phiS = (Se - 1) ./ (Se.^(N+1) -1); phiS(1:3) = 0;

% absorber equation (5-54) page 186 to find v
V = zeros(1, N);
for i = 1:N
  V(i) = V_in(i).*phiA(i) + L_in(i+1).*(1 - phiS(i));
end
  
% absorber overall mass balance to find l
L = zeros(1, N);
for i = 1:N
  L(i) = L_in(i+1) + V_in(i) - V(i);  
end

% print results
fprintf('\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n', ... 
s(1,:), s(2,:), s(3,:), s(4,:), s(5,:), s(6,:), s(7,:));
fprintf('K\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n', ... 
K(1), K(2), K(3), K(4), K(5), K(6));
fprintf('A\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n', ... 
A(1), A(2), A(3), A(4), A(5), A(6));
fprintf('S\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n', ... 
S(1), S(2), S(3), S(4), S(5), S(6));
fprintf('pA\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n', ... 
phiA(1), phiA(2), phiA(3), phiA(4), phiA(5), phiA(6));
fprintf('pS\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n', ... 
phiS(1), phiS(2), phiS(3), phiS(4), phiS(5), phiS(6));
fprintf('S\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n', ... 
S(1), S(2), S(3), S(4), S(5), S(6));
fprintf('v\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n', ... 
V(1), V(2), V(3), V(4), V(5), V(6), sum(V));
fprintf('l\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\t%1.2f\n', ... 
L(1), L(2), L(3), L(4), L(5), L(6), sum(L));
