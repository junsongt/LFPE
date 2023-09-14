% Examples


%% Example 1: y = sin(x)
close; clear; clc;
x = linspace(0,2*pi,100);
y = sin(x);
plot(x,y,'.')

%% Example 2: plot solution of y' = -t*y^2, y(0) = 1
% y(t) = 1/(t^2/2 + 1/y0)
close; clear; clc;
t = linspace(0,5,100);
y = 1./(t.^2+1);
plot(t,y,'.')


%% Example 3: plot solution of y' = -t*y^2, y(0) = y0
% y(t) = 1/(t^2/2 + 1/y0)
close; clear; clc;
t = linspace(0,2,100);
for y0 = 1:0.1:10
    y = 1./(t.^2 + 1/ y0);
    plot(t,y), hold on
end


%% Example 4 Directional field plot
close; clear; clc;
f = @(t,y) -t.*y;
t = 0:0.2:4; y = -2:0.2:2; y0 = -2:0.2:2;
dirfield(f,t,y,y0);


