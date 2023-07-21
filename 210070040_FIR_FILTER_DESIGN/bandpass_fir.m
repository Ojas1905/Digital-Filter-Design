fs = 600e3;
fs1 = 54e3;
fp1 = 59e3;
fp2 = 134e3;
fs2 = 139e3;
d = 0.15;

fc1 = (fs1+fp1)/2;
fc2 = (fp2+fs2)/2;

Wc1 = (fc1*2*pi)/fs;
Wc2  = (fc2*2*pi)/fs;

%Kaiser paramters
A = -20*log10(0.15);
if(A < 21)
    beta = 0;
elseif(A <51)
    beta = 0.5842*(A-21)^0.4 + 0.07886*(A-21);
else
    beta = 0.1102*(A-8.7);
end

N_min = ceil((A-8) / (2.285*0.01667*pi));           

bp = ideal_lp(Wc2,N_min) - ideal_lp(Wc1,N_min);    %Ideal bandpass impulse response of length N_min

window = (kaiser(N_min,beta))';     %Kaiser Window of length N_min

FIR_BP = bp .* window;
fvtool(FIR_BP);         %frequency response

[H,f] = freqz(FIR_BP,1,1024, fs);   %magnitude response
plot(f,abs(H))
grid

function h = ideal_lp(wc,M)

    alpha = (M-1)/2;
    n = 0:M-1;
    m = n - alpha +eps;
    h = sin(wc*m) ./ (pi*m);

end