fs = 425e3;
fs1 = 59e3;
fp1 = 54e3;
fp2 = 104e3;
fs2 = 99e3;

fc1 = (fs1+fp1)/2;
fc2 = (fp2+fs2)/2;

Wc1 = (fc1*2*pi)/fs;
Wc2  = (fc2*2*pi)/fs;

A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-8) / (2.285*0.0235*pi));         

bs = ideal_lp(pi,N_min) - ideal_lp(Wc2,N_min) + ideal_lp(Wc1,N_min);   %Ideal bandstop impulse response of length N_min

window = (kaiser(N_min,beta))'; %Kaiser Window of length N_min

FIR_BS = bs .* window;
fvtool(FIR_BS);         %frequency response

[H,f] = freqz(FIR_BS,1,1024, fs);   %magnitude response
plot(f,abs(H))
grid

function h = ideal_lp(wc,M)

    alpha = (M-1)/2;
    n = 0:M-1;
    m = n - alpha + eps;
    h = sin(wc*m) ./ (pi*m);

end
