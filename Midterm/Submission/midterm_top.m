%% Author Nicholas DiGregorio, 1220871392
%% midterm top level 
clear;
close all; 

%% Parameters

Nb       = 16;    %block size     
mu       = .0051; %step size
N_points = 1008;  %number of points
eflag    = 0;     %1 for NEE and 0 for MSE
s        = Nb;    %shifting parameter
d_flag   = 1;     %1 for white noise input, 0 for colored noise input


%% Unknown systems
%b = [1 1 1 1];
b=[0.776  2.397  1.966  1.859  1.171  0.123  0.525  -0.994  0.588  -1.177  -0.102  1.471 3.161  4.329  2.023  2.666];
cb = [1];
ca = [1 0.64];
try
    bpad = zeros(1, Nb -length(b));
catch
    disp("b is already correct size")
end
b = [b, bpad];
a = [1];

%% Create the desired signal 
rng(3) % seed the random number generator to produce the same numbers
x = randn(N_points, 1); %white noise


d1 = filter(b, a, x); %nominal output

cd = filter(cb, ca, x); %colored noise
d2 = filter(b, a, cd); %colored output

if d_flag == 1
    d = d1; %desired signal for white noise input
else %d_flag == 0
    d = d2;
end

%% BLMS
blms(Nb, mu, x, d, N_points, b, eflag, s)


%% Fast LMS 
flms(Nb, mu, x, d, N_points, b, eflag, s)