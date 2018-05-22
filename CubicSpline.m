function [] = CubicSpline(A, B, n)

%% Define Runge Function
f = @(x) 1./(1 + 25*x.^2);
%% Define Runge Function Derivative
g = @(x) -(50*x)./(25*x.^2+1).^2;
%% Define Error Function
e = @(fx, sx) abs((fx - sx)/fx);

%% Create window to plot Runge Function
%% and Cubic Spline Function
fig = figure(1);
clf();

%% Plot runge function points
x_runge = linspace(A,B,150);
y_runge = f(x_runge);


%% Set natural spline

now = 1;

%% 2) Divide range to desired subintervals
data_points = linspace(A,B,n);
f_data_points = f(data_points);
data_points_size = length(data_points);

offset = (B-A)/(n+1);
x_test = [];
y_f = [];
y_s = [];
emax = 0;
y_approx = [];
%% 3) Iterasi per subinterval
%%        - plot error
for i = 1 : data_points_size - 1
  low = data_points(i);
  high = data_points(i+1);

  %% calculate spline function
  s = calculateSpline(low, high, f, g);
  
  %% plot points in the subinterval
  while x_runge(now) <= high
    y_approx = horzcat(y_approx, [s(x_runge(now))]);
    now++;
    length(y_approx);
    if now > 150
      break;
    endif;
  endwhile;
  
  x_test(i) = low + offset; 
  f_y = f(x_test(i));
  s_y = s(x_test(i));
  y_f(i) = f_y;
  y_s(i) = s_y;
  err = e(f_y, s_y);
  if err > emax
    emax = err;
  end;
endfor;


%% Plot Runge function to graph
plot(x_runge, y_runge, 'g', 'linewidth', 3);
hold on;

%% Plot interpolation
plot(x_runge, y_approx, 'b', 'linewidth', 2);
hold on;

plot(data_points,f_data_points,'k.','markersize',15);
title('Runge Function Cubic Spline Interpolation','fontsize',11);
legend('Runge Function', 'Cubic Spline', 'Data Points');
ylim([-0.25,1.25]);

figure(2);
clf();
plot(x_test, y_s - y_f);
title('Error in Runge Cubic Spline Interpolation to |(f(x) - s(x))/f(x)|');
print -dpng ErrorRungeCubicSpline;
endfunction

function s = calculateSpline(XL, XR, f, f_drv)
  
  delta_x = XR - XL;
  y2 = (f(XR) - f(XL)) / delta_x;
  
  a = f(XL);
  b = f_drv(XL); 
  c = (y2 - f_drv(XL)) / delta_x;
  d = (f_drv(XR) + f_drv(XL) - 2*y2) / (delta_x^2);
  
  s = @(x) a + b*(x - XL) + c*(x - XL)^2 + d*((x - XL)^2)*(x - XR);
  
endfunction