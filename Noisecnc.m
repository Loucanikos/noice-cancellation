clear; close all; clc;
U = 2000;
P = 50; Ts = 1/P; %Ts is the standard sampling time which is 0.02 secs
T = [0:U-1]*Ts;
f1 = 1;
s = sin(2*pi*f1*T); 
x = randn(1,U); %Generating white noise

G = filter([1 0.5],1,x);  %generating noise G(T) which is correlated to x(T)
D = s + G;  %generating noisy sinewave
figure
subplot(411), plot(T,s), title('Signal') 
subplot(412), plot(T,G), title('Background Noise') 
subplot(413), plot(T,D), title('Background Noise and signal combined')        

%LMS algorithm
L = 2;
w = [-4 zeros(1,L-1)]';
xx = [zeros(1,L-1) x]; 
mu = 0.3;
for i = 1:U
y(i) = w(:,i)'*xx(i+L-1:-1:i)'; 
e(i) = D(i) - y(i);
w(:,i+1) = w(:,i) + mu*e(i)*xx(i+L-1:-1:i)'; 
end
subplot(414), plot(T,e), title('Signal Using Noise Cancellation')

% MSE performance curve
z = w(1,:);
npts = 40;
scale = 10; 
variance = var(D);
r_dx = xcorr(D,x,0,'unbiased'); 
R_xx = var(x);
for w11 = -npts:npts
w1 = w11/scale;
J(w11+npts+1) = variance - 2*r_dx*w1 + R_xx*w1.^2 ;
end
figure, plot(-npts/scale:1/scale:npts/scale,J), xlabel('w'), ylabel('MSE'), 
title('(MSE)Mean Square Error Performance Curve')
for p = 1:10:length(z)
J1(p) = variance - 2*r_dx*z(p) + R_xx*z(p).^2 ; 
hold on, plot(z(p),J1(p),'or') , hold off 
pause(0.1)
end

% Using a filter to filter out noise.
L = 11;
k = -(L-1)/2:(L-1)/2;
fc = f1; 
Fc = fc/fs;
h = 2*Fc*sinc(2*k*Fc); 
s_filter = filter(h,1,D);

% Using the optimum filter from LMS
wopt = w(:,end);
v_est = filter(wopt,1,x); 
s_est = D - v_est;
figure
subplot(311), plot(T,x), title('Noisy signal')
subplot(312), plot(T,s_filter), title('Estimated signal after filtering using a lowpass filter')
subplot(313), plot(T,s_est), title('Estimated signal with noise canceller')

clear all; close all; clc

d = audioread('noisy_speech.wav')'; 
x = audioread('noise.wav')';
N = length(d);
L = 2;
w = [1 0]';
xx = [zeros(1,L-1) x]; 
mu = 0.001;
for i = 1:N
y(i) = w(:,i)'*xx(i+L-1:-1:i)'; 
e(i) = d(i) - y(i);
w(:,i+1) = w(:,i) + mu*e(i)*xx(i+L-1:-1:i)';
end

% Using a lowpass filter to filter out noise. 
fs = 8000; Ts = 1/fs;
t = [0:N-1]*Ts; 
L = 11;
k = -(L-1)/2:(L-1)/2; 
fc = 1000;
Fc = fc/fs;
h = 2*Fc*sinc(2*k*Fc); 
s_filter = filter(h,1,d);

% Using the optimum filter from LMS 
wopt = w(:,end);
est_noise = filter(wopt,1,x); 
s_est = d - est_noise;
figure
subplot(311), plot(t,d), title('Noisy signal')
subplot(312), plot(t,s_filter), title('Estimated signal after filtering using a lowpass filter')
subplot(313), plot(t,s_est), title('Estimated signal with noise canceller')
