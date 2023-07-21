w1 = 59;
w2 = 134;
w3 = 54;
w4 = 139;
fs = 600;
d = 0.15;

wp1 = tan(w1*pi/fs);
wp2 = tan(w2*pi/fs);
ws1 = tan(w3*pi/fs);
ws2 = tan(w4*pi/fs);

w0 = sqrt(wp1*wp2);
B = wp2-wp1;

wls1 = calculatewl(ws1, w0, B);
wlp1 = calculatewl(wp1, w0, B);
wlp2 = calculatewl(wp2, w0, B);
wls2 = calculatewl(ws2, w0, B);

wls = min(abs(wls1), abs(wls2));

ep = sqrt((2*d - d*d)/(1 - 2*d + d*d));
k = 1/wls;
k_ = sqrt(1-(k*k));
k1 = ep/sqrt((1/d^2)-1);
k1_ = sqrt(1-(k1*k1));

[K, K_] = ellipke(k);
[K1, K1_] = ellipke(k1);

N = ceil((K1_*K)/(K_*K1));

k = ellipdeg(N, k1);
wls_new = 1/k;

L = floor(N/2);
r = mod(N,2);

v0 = (-1i/(N))*asne((1i/ep), k1);

u = 1/N;

p1 = 2*pi*1i*cde((u-1i*v0), k);
p0 = 2*pi*1i*sne((1i*v0), k);
p1_ = conj(p1);

zeta = cde(u, k);
zi = 2*pi*wlp2*1i./(k*zeta);
z_ = conj(zi);

%poles = [p0 p1];

[Nr, Dr] = zp2tf([zi z_]', [p0 p1 p1_], 0.85);

syms s z;
H_analog_lpf(s) = poly2sym(Nr, s)/poly2sym(Dr,s); 
sL = (w0*w0 + s*s)/(B*s);                                % Converting function
H_analog_bpf(s) = H_analog_lpf(sL);              % Analog bandstop
H_discrete_bpf(z) = H_analog_bpf((1-(z^-1))/(1+(z^-1)));          % Discrete (bilinear transformation

% Coefficients of analog transfer function
[Nr_analog, Dr_analog] = numden(H_analog_bpf(s)); 
Nr_analog1 = sym2poly(Nr_analog);
Dr_analog1 = sym2poly(Dr_analog);
Dr_analog2 = Dr_analog1/Dr_analog1(1);
Nr_analog2 = Nr_analog1/Dr_analog1(1);

% Coefficients of discrete transfer function
[Nr_discrete, Dr_discrete] = numden(H_discrete_bpf(z));
Nr_discrete1 = sym2poly(Nr_discrete);
Dr_discrete1 = sym2poly(Dr_discrete);
Dr_discrete2 = Dr_discrete1/Dr_discrete1(1);
Nr_discrete2 = Nr_discrete1/Dr_discrete1(1);
fvtool(Nr_discrete2, Dr_discrete2)              % Plotting the frequency response and phase response

[H,f] = freqz(Nr_discrete2, Dr_discrete2,1024*1024, 600e3);    % Magnitude response
plot(f,abs(H))
%phasez(Nr_discrete2, Dr_discrete2,1024*1024, 600e3)
grid

function wl = calculatewl(w, w0, B)
    wl = (w^2 - w0^2)/(B*w) ;
end