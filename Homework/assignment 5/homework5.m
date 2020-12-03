figure(1)
A = 4*pi*10^-6;
Cm = 1;
gCl = 0.3;
tau = Cm/gCl*10^-3;
dt = 1/8*10^-3;
N = 8192;
T = N*dt;
df = 1/T;
F = N*df;
t = linspace(0,T,N);
f = linspace(0,F,N);
w = exp(-t/tau).*(t>=0);
w_hat = fft(w)*dt;
m = abs(w_hat)/A/Cm;
phi = angle(w_hat)*180/pi;
subplot(1,2,1)
plot(f,m,'k')
xlabel('\omega (Hz)')
ylabel('R_{in} (M\Omega)')
%xlim([0 100]) %to produce same figures as Fig 3.1
box off
subplot(1,2,2)
plot(f,phi,'k')
xlabel('\omega (Hz)')
ylabel('phase (deg)')
%xlim([0 100]) %to produce same figures as Fig 3.1
box off
figure (2)
rc = 0.24;
rs = 0.96;
c = 0.06;
kc = 1;
ks = kc*c;
fx = 0.01:0.01:2;
u_hat = kc*rc^2*pi*exp(-(rc*pi*fx).^2) - ks*rs^2*pi*exp(-(rs*pi*fx).^2);
scale = 50/max(u_hat);
dx = 0.006;
N = 2048;
x_end = dx*N/2;
x_start = -dx*(N/2-1);
x = x_start:dx:x_end;
u = kc*rc*sqrt(pi)*exp(-(x/rc).^2)-ks*rs*sqrt(pi)*exp(-(x/rs).^2);
u = scale*u;
u_t(1:1025) = u(1024:2048);
u_t(2048:-1:1026)= u(1023:-1:1);
subplot(1,2,2)
plot(x,u,'k')
box off
xlim([-7 7])
xticks([-5 0 5])
xlabel('Spatial angle (deg)')
yticks([0 50 100])
ylabel('Firing frequency(ref.spont.)(spks/s)')
df = 1/dx/N;
u_hat_t = fft(u_t)*dx;
fx_t= (0:1:1025)*df;
subplot(1,2,1)
u_hat = u_hat*scale;
u_hat_t = real(u_hat_t);
fx_g = [0.1 1 2];
ft = 3;
dt_g = 1;
T_g = 1000;
Nt = T_g/dt_g;
t_g = linspace(0,T_g,Nt);
for j= 1:3
  g = zeros(1000,2048);
  
  for i = 1:Nt
    g(i,1:2048) = cos(2*pi*fx_g(j)*x - 2*pi*ft*t_g(i)*10^-3);
  end
  c = g*u'*dx;
  
  cx(j) = max(c);
end

loglog(fx,u_hat,'k--',fx_t,u_hat_t(1:1026),'r--');
hold on
plot(fx_g,cx,'Marker','o','Color','k','LineStyle','none')
hold off
xlabel('Spatial frequency (c/deg)')
xlim([10^-2 10^1])
xticks([10^-2 10^0])
ylabel('Contrast sensitivity')
ylim([10^0 10^2])
yticks([10^0 10^1 10^2])
