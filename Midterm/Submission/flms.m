%% Author Nicholas DiGregorio, 1220871392
function bhat = flms(Nb, mu, x, d, N_points, b, eflag, s)

%Nb = block size
%mu = step size
%x  = input signal
%d  = desired signal 
%N_points = number of points in the signal 



%% Variable creation



%eflag = 0; 
%algorithm controls
%mu = .001;
%Nb = 8;
%s = Nb; % shifting parameter size




%% Algorithm
 
Ek   = zeros(2*Nb, 1);
e    = zeros(N_points,1);
E    = zeros(floor((N_points-Nb)/s),1);

Bhat = zeros(2*Nb,1); 
zpad = zeros(1, Nb); 
b = [b, zpad];

tic
icount = 0;
for i = 1: s : N_points- 2*Nb % i is the block index
    icount = icount + 1; 
    
    %get the input block 
    xv = x(i : (i-1)+Nb*2); %the ith block and the jth sample to the jth+L sample
    db = d(Nb + (i-1) : Nb + (i-1) + (2*Nb-1)); %get the current block, which to start is the Nbth d sample
    
    %do the ffts
    Xk = fft(xv);    
    Dk = fft(db);

    %calculate output
    Xkd = diag(Xk);
    Yk = Xkd * Bhat; % Nbx(L+1) * (L+1)x1
    
    %calculate the error
    yb = ifft(Yk);
    yn = yb(Nb+1: Nb*2); %grabbing the last N points
    dn = d(Nb+1 + (i-1) : Nb+1 + (i-1) + (Nb-1)); %get the Nbth output from desired signal
    en = dn - yn; %generating the error in time domain
    e(Nb + (i-1) : Nb + (i-1) + (Nb-1)) = en;
    eb = [zpad'; en]; %pad the error signal before ffting
    Ek = fft(eb); 


    %calc gradient + update weights
    Bhat = Bhat + 2*mu*Xkd'*Ek ;


    if eflag == 1
        B = fft(b);
        num = nee_num(B, Bhat, Nb); %freq of B - freq of Bhat
        denom = nee_denom(B, Nb);
        E(icount) = trapz(num)/trapz(denom);
    end

end
b = ifft(Bhat);
bhat = b(1:Nb);
elapsed_time = toc;
disp(['FLMS time: ', num2str(elapsed_time), ' seconds for Nb=', num2str(Nb)]);

figure;
title_str = ['FLMS: Nb=', num2str(Nb), ' Mu=', num2str(mu)];
if eflag == 1
    plote = 10*log10(E);
    plot(plote(:));
    title(title_str);
    ylabel('NEE (dB)');
    xlabel('Iteration');
else
    plote = 10*log10(e.*e +0.0000001);
    plot(plote(:));
    title(title_str);
    ylabel('MSE (dB)');
    xlabel('Iteration');
end