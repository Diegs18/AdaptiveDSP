%% Author Nicholas DiGregorio, 1220871392
function bhat = blms(Nb, mu, x, d, N_points, b, eflag, s)

%% Variable creation

%N_points = 6000;
%eflag = 1;
%algorithm controls
%mu = .001;
%Nb = 16;
L  = Nb-1;%Nb-1;
coef = L+1;
%s = Nb/2; % shifting parameter size



%% Algorithm
xb   = zeros(Nb, coef);
xbb   = zeros(Nb, coef);
yb   = zeros(coef, 1); 
eb   = zeros(coef, 1);
e    = zeros(N_points,1);
E    = zeros(floor((N_points-Nb)/s),1);
bhat = zeros(coef,1);


tic
icount = 0;
for i = 1: s : N_points- (3*coef) % i is the block index
    icount = icount + 1; 
    %disp(icount)
    %create the input block 
    for j = 1:Nb
        xv = (x(i+j-1 : i+j-1+L)); %the ith block and the jth sample to the jth+L sample
        xbb(j, :) = xv';
        xb(j, :) = flip(xv');
    end
    xbt = xb';
    %calculate output
    yb = xb * bhat; % Nbx(L+1) * (L+1)x1

    %calculate the error
    db = d(coef + i-1 : Nb + coef + i-2); %get the current block %this is where I had all my issues, since I am starting the convolution at the Nbth sample I need to get the Nbth desired sample
    eb = db - yb;
    e(coef + i-1 : Nb + coef + i-2) = eb; %save the error
       
    %calc gradient
    grad =  xb' * eb;
    
    %update weights
    bhat = bhat + (2) * mu * grad;

    if eflag == 1
        B = fft(b);
        Bhat = fft(bhat);
        num = nee_num(B, Bhat, Nb); %freq of B - freq of Bhat
        denom = nee_denom(B, Nb);
        E(icount) = trapz(num)/trapz(denom);
    end

end
elapsed_time = toc;
disp(['BLMS time: ', num2str(elapsed_time), ' seconds for Nb=', num2str(Nb)]);

figure;
title_str = ['BLMS: Nb=', num2str(Nb), ' Mu=', num2str(mu)];
if eflag == 1
    plote = 10*log10(E);
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
