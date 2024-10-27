%% Variable creation
clear;
N_points = 6000;
eflag = 1;
%algorithm controls
mu = .00051;
Nb = 8;
L  = Nb-1;
coef = L+1;
s = Nb; % shifting parameter size


%% Create the desired signal 
b = [1 1 1 1 1 1 1 1];
a = [1];
rng(3) % seed the random number generator to produce the same numbers
x = randn(N_points, 1);


d = filter(b, a, x);

%% Algorithm
xb   = zeros(Nb, coef);
yb   = zeros(coef, 1); 
eb   = zeros(coef, 1);
e    = zeros(N_points,1);
E    = zeros(floor((N_points-Nb)/s),1);
bhat = zeros(coef,1);


tic
icount = 0;
for i = 1: s : N_points- (2*coef) % i is the block index
    icount = icount + 1; 
    %create the input block 
    for j = 1:Nb
        xv = flip(x(i+j-1 : i+j-1+L)); %the ith block and the jth sample to the jth+L sample
        xb(j, :) = xv';
    end

    %calculate output
    yb = xb * bhat; % Nbx(L+1) * (L+1)x1

    %calculate the error
    db = d(i : (i-1)+Nb); %get the current block
    eb = db - yb;
    e(i : (i-1)+Nb) = eb; %save the error
    
    %calc gradient
    grad = (2/Nb) * xb' * eb;
    
    %update weights
    bhat = bhat + mu * grad;

    if eflag == 1
        B = fft(b);
        Bhat = fft(bhat);
        num = nee_num(B, Bhat, Nb); %freq of B - freq of Bhat
        denom = nee_denom(B, Nb);
        E(icount) = trapz(num)/trapz(denom);
    end

end
elapsed_time = toc;
disp(['Elapsed time: ', num2str(elapsed_time), ' seconds']);

if eflag == 1
    plote = 10*log10(E);
    plot(plote(:));
    ylabel('NEE (dB)')
    xlabel('Iteration')
else
    plote = 10*log10(e.*e +0.0000001);
    plot(plote(:));
    ylabel('MSE (dB)')
    xlabel('Iteration')
end
