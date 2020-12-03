figure(1)
alpha = 100; beta = 0.9;
n1 = 6; n2 = 9;
k = 120; %s^-1
kx = 5;
sx = 0.12;
dx = 0.004;
N = 128;
xstart = -(N/2-1)*dx;
xend = N/2*dx;
x = xstart:dx:xend;
gse = 1/(sqrt(2*pi)*sx)*exp(-x.^2/(2*sx^2)).*cos(2*pi*kx*x);
gso = 1/(sqrt(2*pi)*sx)*exp(-x.^2/(2*sx^2)).*sin(2*pi*kx*x);
plot(x,gse,'k',x,gso,'r')
xlabel('degree')
legend('g_{se}','g_{so}', 'location','best')
box off
figure(2)
dt = 2*10^-3; %s
Tstart = -(N/2-1)*dt;
Tend = N/2*dt;
t = Tstart:dt:Tend;
f1 = alpha*(k*t).^n1.*exp(-k*t).*(1/factorial(n1) - beta*(k*t).^2/factorial(n1+2));
f2 = alpha*(k*t).^n2.*exp(-k*t).*(1/factorial(n2) - beta*(k*t).^2/factorial(n2+2));
plot(t,f1,'k',t,f2,'r')
xlabel('time (s)')
legend('f_1','f_2','location','best')
box off
figure(3)
g = gse.*f1'+gso.*f2';
mesh(x,t,g)
ylabel('time (s)')
xlabel('degree')
figure(4)
nx = 4.5;
nt = 9;
T0 = 2;
T = (0:dt:T0)';
gr = cos(2*pi*(nx*repmat(x,[length(T) 1]) - nt*repmat(T,[1,N])));
gr_e_s = gr*gse';
con = conv(gr_e_s,f1);
gr_e_st = con(64:64+length(T)-1)*dt;
gl = cos(2*pi*(nx*repmat(x,[length(T) 1]) + nt*repmat(T,[1,N])));
gl_e_s = gl*gse';
con = conv(gl_e_s,f1);
gl_e_st = con(64:64+length(T)-1)*dt;
gr_o_s = gr*gso';
con = conv(gr_o_s,f2);
gr_o_st = con(64:64+length(T)-1)*dt;
gl_o_s = gl*gso';
con = conv(gl_o_s,f2);
gl_o_st = con(64:64+length(T)-1)*dt;
g_r = gr_e_st+gr_o_st;
g_l = gl_o_st+gl_o_st;
plot(T,g_r,'b',T,g_l,'r')
legend('right grating','left grating', 'location', 'best')
xlabel('time (s)')