figure(1)
t = -5.11:0.01:5.12;
g = exp(-t.^2/2);
plot(t,g,'k')
xlabel('time (s)')
figure(2)
omega = -50.5:0.195:50.5;
g_hat_t = sqrt(2*pi)*exp(-(2*pi*omega).^2/2);
plot(omega,g_hat_t,'k')
xlabel('frequency (Hz)')
figure(3)
g_t(1:513) = g(512:1024);
g_t(514:1024) = g(1:511);
g_hat = fft(g_t)*0.01;
ghat(1:511) = g_hat(514:1024);
ghat(512:1024) = g_hat(1:513);
ghat = ghat;
sample = length(omega);
c = length(ghat);
start = c/2-sample/2;
g_hat_n = real(ghat(start:start+sample-1));
plot(omega,g_hat_t,'k', omega,g_hat_n,'r--')
xlabel('frequency (Hz)')
xlim([-56 56])
legend('theoretical', 'numerical')
ylim([0 3])
maxi = max(imag(ghat));
display(maxi)
