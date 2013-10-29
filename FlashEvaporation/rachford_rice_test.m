% Test script for rachford_rice.m
close all; clear all; clc;

% feed conditions
F = 100; %kmol/h
% compositions (add z(i) for more components)
z = [0.0079 0.1321 0.0849 0.2690 0.0589 0.1321 0.3151];
% K-values
K = [16.2   5.2    2.6    1.98   0.91   0.72   0.28];

[x y V L Q] = rachford_rice(F, z, K);

% print results ----------------------------------------------------------%
fprintf('Component\t\t  x\t\t  y\n');
for i = 1:length(x)
   fprintf('%d\t\t\t\t%1.3f\t%1.3f\n', i, x(i), y(i));
end
fprintf('Total\t\t\t%1.2f\t%1.2f\n', sum(x), sum(y));

fprintf('\nV = %1.3f\tL = %1.3f\n', V, L);
fprintf('Q = %1.3f\n', Q);