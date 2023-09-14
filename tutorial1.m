% MATH 215 MATLAB Tutorial 1
% September 15, 2022

%% Example 1: Sine wave
x = linspace(0,1,100);
y = sin(2*pi*x);
plot(x,y)

%% Example 2: Plot solutions of y' = -t*y^2, y(0) = 1
%  y(t) = 1/(t^2/2 + 1/y0)
t = linspace(0,5,100);
y = 1./(t.^2/2 + 1);
plot(t,y)

%% Example 3: Plot solutions of y' = -t*y^2, y(0) = y0
%  y(t) = 1/(t^2/2 + 1/y0)
clear; close; clc;
t = linspace(0,2,100);
for y0=1:10
    y = 1./(t.^2/2 + 1/y0);
    plot(t,y,'b'), hold on
end

%% Example 4: Direction field plot
close; clear; clc;
f = @(t,y) -t.*y;
t = 0:0.2:4; y = -2:0.2:2; y0 = -2:1:2;
dirfield(f,t,y,y0);

