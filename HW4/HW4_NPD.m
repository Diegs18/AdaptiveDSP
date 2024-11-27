%% Author Nicholas DiGregorio, 1220871392
clear;
close all; 

N_points = 1000; 
L = 4; 
coef = L+1;
M = 1; 
mu = .05;
eflag = 1;

bhat = zeros(L+1, 1);
ahat = zeros(M+1, 1);
ahat(1) = 1;
chat = zeros(L+M+1, 1);%[bhat; ahat];
e  = zeros(N_points-coef, 1);
NEE = zeros(N_points-coef, 1);
%% Create the desired signal 
b = [1 -1.3435 0.9025];
a = [1 -.45 .55];

rng(3) % seed the random number generator to produce the same numbers
x = randn(N_points, 1); %white noise


d = filter(b, a, x); %nominal output


%% System ID
%i = coef;
icount = 1; 
for i = coef:N_points
    %grab the x vector
    xv = x(i:-1:i-L);
    
    %grab the d vector
    dv = (-d(i-1:-1:i-M));
    %create the u vector by concatentating xv and dv
    uv = [xv; dv];
    
    %calculate the error
    e(i) = d(i) - uv'*chat;

    %calculate the gradient
    grad = 2*e(i)*uv;

    %calculate the weight update
    chat = chat + mu * grad;

    bhat = chat(1:L+1)';
    %ahat = [1; chat(L+2:end)]'; %use this for part a
    ahat = [1; chat(L+2:end); 0; 0; 0]'; %use this for part b
    h1 = fft(b/a);
    h2 = fft(bhat/ahat);
    NEE(i) = 10.*log10(sum(abs(h1-h2).^2)/sum(abs(h1).^2)+1e-16);
    icount = icount + 1;
end
 
bhat
ahat
title_str = ['IIR EEM:   Mu=', num2str(mu), '   L=', num2str(L), '   M=', num2str(M)];
if eflag == 1
    plote = ( NEE);
    plot(plote(:));
    title(title_str);
    ylabel('NEE (dB)')
    xlabel('Iteration')
else
    plote = 10*log10(e.*e +0.0000001);
    plot(plote(:));
    title(title_str);
    ylabel('MSE (dB)')
    xlabel('Iteration')
end
