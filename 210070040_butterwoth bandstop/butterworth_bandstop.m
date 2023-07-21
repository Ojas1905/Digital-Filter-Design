% Analog lowpass filter specs
wc = 1.037;
N = 14;

% Analog bandstop filter specs after transformation tan(w/2)
ws1 = 0.451;
ws2 = 0.88;
wp1 = 0.414;
wp2 = 0.939;

w0 = 0.6234;   % sqrt(w1p1*wp2)
B = 0.525;     % wp1-wp2

% calculating the poles:
x = zeros(1,28);
y = zeros(1,28);

for j=1:28
    m=wc*cos((pi+2*j*pi)/28);
    n=wc*sin((pi+2*j*pi)/28); 
    y(j) = m;
    x(j) = n;
end

re_p = x(14:27);
im_p = y(14:27);

% left half complex plane poles:
poles = [];

for j = 1:14
    poles(j) = re_p(j) + i*im_p(j);
end

% Transfer Function:
[Nr, Dr] = zp2tf([], [poles], wc^N);

% Defining our system
syms s z;
H_analog_lpf(s) = poly2sym(Nr, s)/poly2sym(Dr,s);       % Analog lowpass
sL = (B*s)/(w0^2 + s^2);                                % Converting function
H_analog_bsf(s) = H_analog_lpf(sL);                     % Analog bandstop
H_discrete_bsf(z) = H_analog_bsf((z-1)/(z+1));          % Discrete (bilinear transformation

% Coefficients of analog transfer function
[Nr_analog, Dr_analog] = numden(H_analog_bsf(s)); 
Nr_analog1 = sym2poly(Nr_analog);
Dr_analog1 = sym2poly(Dr_analog);
Dr_analog2 = Dr_analog1/Dr_analog1(1);
Nr_analog2 = Nr_analog1/Dr_analog1(1);

% Coefficients of discrete transfer function
[Nr_discrete, Dr_discrete] = numden(H_discrete_bsf(z));
Nr_discrete1 = sym2poly(Nr_discrete);
Dr_discrete1 = sym2poly(Dr_discrete);
Dr_discrete2 = Dr_discrete1/Dr_discrete1(1);
Nr_discrete2 = Nr_discrete1/Dr_discrete1(1);
fvtool(Nr_discrete2, Dr_discrete2)              % Plotting the frequency response and phase response


[H,f] = freqz(Nr_discrete2, Dr_discrete2, 425e3);    % Magnitude response
plot(f,abs(H))
grid


