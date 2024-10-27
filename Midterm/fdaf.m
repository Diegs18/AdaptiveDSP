%% Variable creation
clear;
N_points = 6000;
eflag = 1; 
%algorithm controls
mu = .001;
Nb = 8;
s = Nb; % shifting parameter size


%% Create the desired signal 
b = [1 1 1 1 1 1];
bpad = zeros(1, Nb -length(b));
b = [b, bpad];
a = [1];
rng(3) % seed the random number generator to produce the same numbers
x = randn(N_points, 1);


d = filter(b, a, x);

%% Algorithm
 
Ek   = zeros(Nb, 1);
e    = zeros(N_points,1);
E    = zeros(floor((N_points-Nb)/s),1);
Bhat = zeros(Nb,1);


tic
icount = 0;
for i = 1: s : N_points- Nb % i is the block index
    icount = icount + 1; 
    
    %get the input block 
    xv = x(i : (i-1)+Nb); %the ith block and the jth sample to the jth+L sample
    db = d(i : (i-1)+Nb); %get the current block
    
    %do the ffts
    Xk = fft(xv);    
    Dk = fft(db);

    %calculate output
    Xkd = diag(Xk);
    Yk = Xkd * Bhat; % Nbx(L+1) * (L+1)x1

    %calculate the error
    Ek = Dk - Yk;
    e(i : ((i-1)+Nb)) = ifft(Ek);
    
    %calc gradient + update weights
    Bhat = Bhat + 2*mu*Xkd'*Ek ;


    if eflag == 1
        B = fft(b);
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