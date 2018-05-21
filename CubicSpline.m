function [] = CubicSpline(A, B, n)

%% Define Runge Function
f = @(x) 1./(1 + 25*x.^2);
%% Define Runge Function Derivative
g = @(x) -(50*x)./(25*x.^2+1).^2;

%% Create window to plot Runge Function
%% and Cubic Spline Function
fig = figure(1);
clf();

%% Plot runge function points
x_runge = linspace(A,B,150);
y_runge = f(x_runge);

now = 1;

%% 2) Divide range to desired subintervals
data_points = linspace(A,B,n);
data_points_size = length(data_points);

%% 3) Iterasi per subinterval
%%        - plot error
for i = 1 : data_points_size - 1
  low = data_points(i);
  high = data_points(i+1);

  %% calculate spline function
  a = f(low)
  b = g(low);
  
  delta_x = high - low;
  y2 = (f(high) - f(low)) / delta_x;
  
  c = (y2 - g(low)) / delta_x;
  d = (g(high) + g(low) - 2*y2) / (delta_x^2);
  
  s = @(x) a + b*(x - low) + c*(x - low)^2 + d*((x - low)^2)*(x - high);
  
  y_approx = [];
  %% plot points in the subinterval
  while x_runge(now) <= high
    y_approx = horzcat(y_approx, [s(x_runge(now))])
    now++
    length(y_approx)
    if now > 150
      break;
    endif;
  endwhile;
    
endfor;

length(x_runge)
length(y_approx)
plot(x_runge, y_approx, 'r', 'linewidth', 4);

%% 4) Plot Runge function to graph
% plot(x_runge, y_runge, 'g', 'linewidth', 1);
title('Runge Function Cubic Spline Interpolation','fontsize',11);
legend('Runge Function', 'Cubic Spline', 'Data Points');
ylim([-0.25,1.25]);

