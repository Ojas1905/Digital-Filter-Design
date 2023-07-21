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

D1 = (1/((1-d)^2))-1;
D2 = (1/(d^2))-1;

N = ceil(acosh(sqrt(D2/D1))/acosh(wls));

p1 = 1i*cosh((acosh(-100i/sqrt(3841)))/7);
p2 = 1i*cosh((acosh(-100i/sqrt(3841))-2i*pi)/7);
p3 = 1i*cosh((1/7)*(acosh(-100i/sqrt(3841))+2i*pi));
p4 = 1i*cosh((1/7)*(acosh(-100i/sqrt(3841))-4i*pi));
p5 = 1i*cosh((1/7)*(acosh(-100i/sqrt(3841))+4i*pi));
p6 = 1i*cosh((1/7)*(acosh(-100i/sqrt(3841))-6i*pi));
p7 = 1i*cosh((1/7)*(acosh(-100i/sqrt(3841))+6i*pi));
p8 = 1i*cosh((1/7)*(acosh(100i/sqrt(3841))));
p9 = 1i*cosh((1/7)*(acosh(100i/sqrt(3841))-2i*pi));
p10 = 1i*cosh((1/7)*(acosh(100i/sqrt(3841))+2i*pi));
p11 = 1i*cosh((1/7)*(acosh(100i/sqrt(3841))-4i*pi));
p12 = 1i*cosh((1/7)*(acosh(100i/sqrt(3841))+4i*pi));
p13 = 1i*cosh((1/7)*(acosh(100i/sqrt(3841))-6i*pi));
p14 = 1i*cosh((1/7)*(acosh(100i/sqrt(3841))+6i*pi));

poles = [p3, p5, p7, p8, p10, p12, p14];

[Nr, Dr] = zp2tf([], poles, (-1*p3*p5*p7*p8*p10*p12*p14));

syms s z;
H_analog_lpf(s) = poly2sym(Nr, s)/poly2sym(Dr,s);       % Analog lowpass
%H_analog_lpf(s) = (-1*p3*p5*p7*p8*p10*p12*p14)/((s-p3)*(s-p5)*(s-p7)*(s-p8)*(s-p10)*(s-p12)*(s-p14));
sL = (w0*w0 + s*s)/(B*s);                                % Converting function
H_analog_bpf(s) = (-1*p3*p5*p7*p8*p10*p12*p14)/((sL-p3)*(sL-p5)*(sL-p7)*(sL-p8)*(sL-p10)*(sL-p12)*(sL-p14));                % Analog bandstop
H_discrete_bpf(z) = H_analog_bpf((1-(z^-1))/(1+(z^-1)));          % Discrete (bilinear transformation

% Coefficients of analog transfer function
[Nr_analog, Dr_analog] = numden(H_analog_bpf(s)); 
Nr_analog1 = sym2poly(expand(Nr_analog));
Dr_analog1 = sym2poly(expand(Dr_analog));
Dr_analog2 = Dr_analog1/Dr_analog1(1);
Nr_analog2 = Nr_analog1/Dr_analog1(1);

% Coefficients of discrete transfer function
[Nr_discrete, Dr_discrete] = numden(H_discrete_bpf(z));
Nr_discrete1 = sym2poly(expand(Nr_discrete));
Dr_discrete1 = sym2poly(expand(Dr_discrete));
Dr_discrete2 = Dr_discrete1/Dr_discrete1(1);
Nr_discrete2 = Nr_discrete1/Dr_discrete1(1);
%fvtool(Nr_discrete2, Dr_discrete2)              % Plotting the frequency response and phase response

[H,f] = freqz(Nr_discrete2, Dr_discrete2,1024*1024, 600e3);    % Magnitude response
plot(f,abs(H))
phasez(Nr_discrete2, Dr_discrete2,1024*1024, 600e3)
grid


function wl = calculatewl(w, w0, B)
    wl = (w^2 - w0^2) / (B*w);
end


