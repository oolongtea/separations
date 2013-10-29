% Optimizing Distillation Column - Heptane/Pentane
%
% Prints conditions and compositions of entire column
%
% Obed Espina

clc; clear all; close all; format long;

% D = distillate (top) & B = bottom (of distillation column)
% p/P/5 = pentane (C5H12) & h/H/7 = heptane (C7H16)

P = 1.01325; %bar, atmospheric pressure

Cp7 = 256.1; %J / K mol, heat capacity of heptane @ boiling T.
Hvap7 = 31.77 * 1000; %J / mol, heat of vaporization for heptane
Cp5 = 125; %J / K mol, heat capacity of pentane @ boiling T.
Hvap5 = 25.79 * 1000; %J / mol, heat of vaporization for heptane

R = 8.314; %J / mol K
bpHep = 370.38076; %K
bpPen = 309.284033; %K
Q = @(flow1, flow2, Tout) ...
     (flow1*(bpHep-Tout)*Cp7 + flow2*Hvap7)/1e6/60/60; %MW
Q2 = @(flow1, flow2, Tout) ...
      (flow1*(bpPen-Tout)*Cp5 - flow2*Hvap5)/1e6/60/60; %MW

T_range = linspace(bpPen,bpHep,100); %K

%vapor pressures
Pvp5 = @(T) exp(10.422 - 26799./(R.*T)); % vapor pressure of pentane
Pvp7 = @(T) exp(11.431 - 35200./(R.*T)); % vapor pressure of heptane
%pentane
x5 = @(T) (P - Pvp7(T)) ./ (Pvp5(T) - Pvp7(T));
y5 = @(T) Pvp5(T).*x5(T)./P;
%heptane
x7 = @(T) P.*(1 - y5(T))./Pvp7(T);
y7 = @(T) 1 - y5(T);
plot(x5(T_range), y5(T_range), x5(T_range), x5(T_range));
grid on; title('Vapor-Liquid Equilibrium of Heptane & Pentane'); 
xlabel('x5'), ylabel('y5');
axis([0 1 0 1]);

% Feed
F = 5000; %mol/hr, feed stream
xH0 = 0.5; % 50% liquid H in feed stream
xP0 = 0.5; % 50% liquid P in feed stream
feedH = F*xH0; %mol, moles of H in feed
feedP = F*xP0; %mol, moles of P in feed

% Distillate
yHD = 0.98; % 95% H gas 
yPD = 1 - yHD; % 05% P gas 
distH = feedH * yHD; %mol, moles of H 
distP = feedP * yPD; %mol, moles of P 

xHD = x7(fzero(@(T) (y7(T) - yHD), 360));
xPD = x5(fzero(@(T) (y5(T) - yPD), 360));

% Bottom
xHB = 0.05; % 05% H gas 
xPB = 0.95; % 95% P gas 
botmH = feedH * xHB; %mol, moles of H 
botmP = feedP * xPB; %mol, moles of P 



D = 2500;
B = 2500;
% 
recycleRatio = [1 5 10 100]; %recycle ratio
f = [2 3 4 5]; % feed entry at stage 2, 3, 4, or 5

fprintf('\n');
    fprintf('NOTE: All Temperatures in Kelvin\n');
    fprintf('& Heat Duty in MegaWatts\n');
for q = recycleRatio
    fprintf('\n');
    fprintf('------------------------------------------------\n');
    fprintf('---------------recycle ratio: %d --------------\n', q);
    L = q*D; 
    V = L + D;
    
    for feedEntry = f

if (feedEntry == 2) 
%% Entry at Stage 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x(1) = 0.95; %x1 = x(1)
y(1) = 0.95; y1 = y(1);
t(1) = fzero(@(T) (y5(T) - 0.95), 360);
fprintf('T1 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(1), bpPen-t(1), Q2(V, L, t(1)));
% Stage 1 %
t(2) = fzero(@(T) (y5(T) - y(1)), 360);
x(2) = x5(t(2)); 
y(2) = q*x(2)/(q+1) + y(1)/(q+1); 
if(x(2) < 0.05), 
    fprintf('T2 = %f\t\tdeltaT = %f\n', t(2), bpHep-t(2));
end
% Stage 2: Feed Entry %
t(3) = fzero(@(T) (y5(T) - y(2)), 360);
x(3) = x5(t(3)); 
y(3) = (V*y(2) + (L+F)*x(3) - F*xP0 - L*x(2))/V; 
if(x(3) < 0.05), 
    fprintf('T3 = %f\t\tdeltaT = %f\n', t(3), bpHep-t(3)); 
end
% Stage 3 %
t(4) = fzero(@(T) (y5(T) - y(3)), 360);
x(4) = x5(t(4)); 
y(4) = (V*y(3) + (L+F)*x(4) - (L+F)*x(3))/V;
if(x(4) < 0.05), 
    fprintf('T4 = %f\t\tdeltaT = %f\n', t(4), bpHep-t(4)); 
end
% Stage 4 %
t(5) = fzero(@(T) (y5(T) - y(4)), 360);
x(5) = x5(t(5)); 
y(5) = (V*y(4) + (L+F)*x(5) - (L+F)*x(4))/V;
if(x(5) < 0.05), 
    fprintf('T5 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(5), bpHep-t(5), Q((L+F), V, t(5))); 
end
% Stage 5 %
t(6) = fzero(@(T) (y5(T) - y(5)), 360);
x(6) = x5(t(6)); 
y(6) = (V*y(5) + (L+F)*x(6) - (L+F)*x(5))/V;
if(x(6) < 0.05), 
    fprintf('T6 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(6), bpHep-t(6), Q((L+F), V, t(6))); 
end

stgNum = 6;
fprintf('Pentane Balances:\n\n');
for k = 1:stgNum
    fprintf('\t\tx: %f\ty: %f\n', x(k), y(k));
    if k < stgNum, fprintf('---Stage %d----------------------\n', k); end
end
fprintf('\n\twith feed of %d mol/hr', F);
fprintf('\n\tand feed entrance at Stage 2\n');

end

if (feedEntry == 3) 
%% Entry at Stage 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n');
x(1) = 0.95; %x1 = x(1)
y(1) = 0.95; y1 = y(1);
t(1) = fzero(@(T) (y5(T) - 0.95), 360);
fprintf('T1 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(1), bpPen-t(1), Q2(V, L, t(1)));
% Stage 1 %
t(2) = fzero(@(T) (y5(T) - y(1)), 360);
x(2) = x5(t(2)); 
y(2) = q*x(2)/(q+1) + y(1)/(q+1);
if(x(2) < 0.05), 
    fprintf('T2 = %f\t\tdeltaT = %f\n', t(2), bpHep-t(2));
end
% Stage 2 %
t(3) = fzero(@(T) (y5(T) - y(2)), 360);
x(3) = x5(t(3)); 
y(3) = (V*y(2) + (L)*x(3) - (L)*x(2))/V;
if(x(3) < 0.05), 
    fprintf('T3 = %f\t\tdeltaT = %f\n', t(2), bpHep-t(3));
end
% Stage 3: Feed Entry %
t(4) = fzero(@(T) (y5(T) - y(3)), 360);
x(4) = x5(t(4)); 
y(4) = (V*y(3) + (L+F)*x(4) - F*xP0 - L*x(3))/V; 
if(x(4) < 0.05), 
    fprintf('T4 = %f\t\tdeltaT = %f\n', t(4), bpHep-t(4));
end
% Stage 4 % 
t(5) = fzero(@(T) (y5(T) - y(4)), 360);
x(5) = x5(t(5)); 
y(5) = (V*y(4) + (L+F)*x(5) - (L+F)*x(4))/V; 
if(x(5) < 0.05), 
    fprintf('T5 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(5), bpHep-t(5), Q((L+F), V, t(5))); 
end
% Stage 5 %
t(6) = fzero(@(T) (y5(T) - y(5)), 360);
x(6) = x5(t(6)); 
y(6) = (V*y(5) + (L+F)*x(6) - (L+F)*x(5))/V; 
if(x(6) < 0.05), 
    fprintf('T6 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(6), bpHep-t(6), Q((L+F), V, t(6))); 
end
stgNum = 6;
fprintf('Pentane Balances:\n\n');
for k = 1:stgNum
    fprintf('\t\tx: %f\ty: %f\n', x(k), y(k));
    if k < stgNum, fprintf('---Stage %d----------------------\n', k); end
end
fprintf('\n\twith feed of %d mol/hr', F);
fprintf('\n\tand feed entrance at Stage 3\n');
end

if (feedEntry == 4) 
%% Entry at Stage 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n');
x(1) = 0.95; %x1 = x(1)
y(1) = 0.95; y1 = y(1);
t(1) = fzero(@(T) (y5(T) - 0.95), 360);
fprintf('T1 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(1), bpPen-t(1), Q2(V, L, t(1)));
% Stage 1 %
t(2) = fzero(@(T) (y5(T) - y(1)), 360);
x(2) = x5(t(2)); 
y(2) = q*x(2)/(q+1) + y(1)/(q+1); 
if(x(2) < 0.05), 
    fprintf('T2 = %f\t\tdeltaT = %f\n', t(2), bpHep-t(2));
end
% Stage 2 %
t(3) = fzero(@(T) (y5(T) - y(2)), 360);
x(3) = x5(t(3)); 
y(3) = (V*y(2) + (L)*x(3) - (L)*x(2))/V;
if(x(3) < 0.05), 
    fprintf('T3 = %f\t\tdeltaT = %f\n', t(3), bpHep-t(3));
end
% Stage 3 %
t(4) = fzero(@(T) (y5(T) - y(3)), 360);
x(4) = x5(t(4)); 
y(4) = (V*y(3) + (L)*x(4) - (L)*x(3))/V; 
if(x(4) < 0.05), 
    fprintf('T4 = %f\t\tdeltaT = %f\n', t(4), bpHep-t(4));
end
% Stage 4: Feed Entry %
t(5) = fzero(@(T) (y5(T) - y(4)), 360);
x(5) = x5(t(5)); 
y(5) = (V*y(4) + (L+F)*x(5) - F*xP0 - L*x(4))/V; 
if(x(5) < 0.05), 
    fprintf('T5 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(5), bpHep-t(5), Q((L+F), V, t(5))); 
end
% Stage 5 %
t(6) = fzero(@(T) (y5(T) - y(5)), 360);
x(6) = x5(t(6)); 
y(6) = (V*y(5) + (L+F)*x(6) - (L+F)*x(5))/V;
if(x(6) < 0.05), 
    fprintf('T6 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(6), bpHep-t(6), Q((L+F), V, t(6))); 
end
% Stage 6 %
t(7) = fzero(@(T) (y5(T) - y(6)), 360);
x(7) = x5(t(7)); 
y(7) = (V*y(6) + (L+F)*x(7) - (L+F)*x(6))/V; 
if(x(7) < 0.05), 
    fprintf('T7 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(7), bpHep-t(7), Q((L+F), V, t(7))); 
end

stgNum = 7;
fprintf('Pentane Balances:\n\n');
for k = 1:stgNum
    fprintf('\t\tx: %f\ty: %f\n', x(k), y(k));
    if k < stgNum, fprintf('---Stage %d----------------------\n', k); end
end
fprintf('\n\twith feed of %d mol/hr', F);
fprintf('\n\tand feed entrance at Stage 4\n');
end

if (feedEntry == 5) 
%% Entry at Stage 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n');
x(1) = 0.95; %x1 = x(1)
y(1) = 0.95; y1 = y(1);
t(1) = fzero(@(T) (y5(T) - 0.95), 360);
fprintf('T1 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(1), bpPen-t(1), Q2(V, L, t(1)));
% Stage 1 %
t(2) = fzero(@(T) (y5(T) - y(1)), 360);
x(2) = x5(t(2)); 
y(2) = q*x(2)/(q+1) + y(1)/(q+1); 
if(x(2) < 0.05), 
    fprintf('T2 = %f\t\tdeltaT = %f\n', t(2), bpHep-t(2));
end
% Stage 2 %
t(3) = fzero(@(T) (y5(T) - y(2)), 360);
x(3) = x5(t(3)); 
y(3) = (V*y(2) + (L)*x(3) - (L)*x(2))/V; 
if(x(3) < 0.05), 
    fprintf('T3 = %f\t\tdeltaT = %f\n', t(3), bpHep-t(3));
end
% Stage 3 %
t(4) = fzero(@(T) (y5(T) - y(3)), 360);
x(4) = x5(t(4)); 
y(4) = (V*y(3) + (L)*x(4) - (L)*x(3))/V; 
if(x(4) < 0.05), 
    fprintf('T4 = %f\t\tdeltaT = %f\n', t(4), bpHep-t(4));
end
% Stage 4 %
t(5) = fzero(@(T) (y5(T) - y(4)), 360);
x(5) = x5(t(5)); 
y(5) = (V*y(4) + (L)*x(5) - (L)*x(4))/V; 
if(x(5) < 0.05), 
    fprintf('T5 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(5), bpHep-t(5), Q((L+F), V, t(5))); 
end
% Stage 5: Feed Entry %
t(6) = fzero(@(T) (y5(T) - y(5)), 360);
x(6) = x5(t(6)); 
y(6) = (V*y(5) + (L+F)*x(6) - F*xP0 - L*x(5))/V; 
if(x(6) < 0.05), 
    fprintf('T6 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(6), bpHep-t(6), Q((L+F), V, t(6))); 
end
% Stage 6 %
t(7) = fzero(@(T) (y5(T) - y(6)), 360);
x(7) = x5(t(7)); 
y(7) = (V*y(6) + (L+F)*x(7) - (L+F)*x(6))/V;
if(x(7) < 0.05), 
    fprintf('T7 = %f\t\tdeltaT = %f\tQ = %f\n', ...
                                t(7), bpHep-t(7), Q((L+F), V, t(7))); 
end

stgNum = 7;
fprintf('Pentane Balances:\n\n');
for k = 1:stgNum
    fprintf('\t\tx: %f\ty: %f\n', x(k), y(k));
    if k < stgNum, fprintf('---Stage %d----------------------\n', k); end
end
fprintf('\n\twith feed of %d mol/hr', F);
fprintf('\n\tand feed entrance at Stage 5\n');
end


    end % end for loop for feed entry point
end % end for loop for all recycle ratios
