function [] = mccabe_thiele(x, y, F, zF, xD, xB, R, q, degree)

%McCabe-Thiele Graphical Method
%
% input:
%  x = liquid equilibrium data for volatile component (light key)
%  y = vapor equilibrium data for volatile component (heavy key)
%  F = feed flow [mole/time]
%  zF = volatile component feed composition [mole fraction]
%  xD = volatile component composition in distillate flow
%  xB = volatile component composition in bottom flow
%  R = L/V = reflux ratio
%  q = L_feed/F if partially vaporized (use q_fun.m)
%    = 1        if bubble-point liquid
%    = 0        if dew-point vapor
%    > 1        if subcooled liquid
%    < 1        if superheated vapor
%  degree = degree used to fit the equilibrium data (uses polyfit)

if nargin<8,error('at least 8 input arguments required'),end
if nargin==8, degree = 3; end
if nargin>9,error('at most 9 input arguments allowed'),end

% fit equilibrium data
p = polyfit(x, y, degree);
x_range = 0:0.01:0.99;
f = polyval(p, x_range); 
% plot(x,y,'r'), hold on; % original data plot
plot(x_range, f, 'g', 0:0.1:1, 0:0.1:1, 'k--'), grid on; 
xlabel('x'), ylabel('y'), title('');
hold on;

% solve for distillate & bottom flows
syms D B
S = solve(F == D + B, 0.4*F == 0.97*D + 0.02*B);
D = S.D; B = S.B;

%slope of q-line
q_slope = q / (q - 1);

% equilibirum function
eq = @(x) polyval(p, x); 
% intersection of equilibrium and q-line
q_point = fzero(@(x) (eq(x) - q_slope.*x + zF/(q-1)), 0.3);

% plot q-line
x_range = 0.42:0.01:0.48; % change range if needed
plot(x_range, q_slope.*x_range - zF/(q-1), 'b'); hold on;

%find intersection between q-line & rectifying
syms x
S = solve(q_slope*x - zF/(q-1) - (R/(R+1))*x - (1/(R+1))*xD);
x_inter = double(S);
y_inter = double(q_slope*x_inter - zF/(q-1));

% plot rectifying op-line
x_range = q_point:0.001:0.97;
plot(x_range, ((R/(R+1)).*x_range + (1/(R+1))*xD), 'r'); hold on;
% plot stripping op-line
x_range = 0:0.001:x_inter;
plot(x_range, (y_inter/x_inter).*x_range, 'r'); hold on;

% plot xD and xB
line([xD xD], [0 1]); line([xB xB], [0 1]); hold on;

% plot the steps of distillation column %
maxiter = 18; % max iterations
x_left = ones(1,maxiter); x_right = xD;
y_bottom = xD; y_top = y_bottom;
% rectifying op-line
rect = @(x) ((R/(R+1))*x + (1/(R+1))*xD);
% stripping op-line
strip =@(x) (y_inter/x_inter).*x;

i = 1; % index counter
x_left(1) = 1;

while (i < maxiter)

  % horizontal lines
  x_left(i) = fzero(@(x) (eq(x) - y_top), x_right);
  line([x_left(i) x_right], [y_bottom y_top]); hold on;

  % vertical lines, depending on side of q-line
  if (x_left > q_point), y_bottom = rect(x_left(i));
  else y_bottom = strip(x_left(i));
  end
  line([x_left(i) x_left(i)], [y_bottom y_top]);
  
  y_top = y_bottom;
  x_right = x_left(i);
  i = i + 1;
end
hold off;

% fprintf('xB = %1.3f\txD = %1.3f\n', x_left(16), xD);
% 
% % solve for distillate & bottom flows
% syms D B
% S = solve(13600 == D + B, 0.4*F == 0.97*D + x_left(maxiter-1)*B);
% D = S.D; B = S.B;
% fprintf('D = %1.1fkg/h\tB = %1.1f kg/h\n', double(D), double(B));
% fprintf('- Bottom composition (mol perc) - \nB: %1.3f\tT: %1.3f\n', ...
%         x_left(maxiter-1), (1 - x_left(maxiter-1)));
% 

end