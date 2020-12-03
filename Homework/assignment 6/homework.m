kt = 8;
st = 31*10^-3; %s
dt = 2 *10^-3; %s
kx = 4.2;
sx = 0.1;
N = 128;
Tend = dt*N/2;
Tstart = -dt*(N/2-1);
T = [Tstart:dt:Tend];
ft = 25;
fx = 10;
ft = linspace(-ft,ft,N);
fx = linspace(-fx,fx,N);
x = linspace(-0.25,0.25,N);
ge = zeros(N);
ge_non = zeros(N);
for j = 1:N
    for i = 1:N
        ge(j,i) = 1/(2*pi*st*sx)*exp(-T(i)^2/(2*st^2)-x(j)^2/(2*sx^2))*cos(2*pi*kt*T(i))*cos(2*pi*kx*x(j));
        ge_non(j,i) = 1/(2*pi*st*sx)*exp(-T(i)^2/(2*st^2)-x(j)^2/(2*sx^2))*cos(2*pi*(kx*x(j)-kt*T(i)));

    end
end
figure(1)
subplot(2,2,1)
meshc(x,T,ge)
colormap('gray')
xticks([-0.2 0 0.2])
xlim([-0.25 0.25])
yticks([-0.1 0 0.1])
ylim([-0.1 0.1])
ylabel('time (s)')
zticks([-50 0 50])
zlim([-55 55])
zlabel('response (arbitrary units)')
grid on
set(gca,'GridLineStyle','--')
subplot(2,2,2)
meshc(x,T,ge_non)
xticks([-0.2 0 0.2])
xlim([-0.25 0.25])
xlabel('space (deg)')
yticks([-0.1 0 0.1])
ylim([-0.1 0.1])
zticks([-50 0 50])
zlim([-55 55])
zlabel('response (arbitrary units)')
set(gca,'GridLineStyle','--')
ge_hat = zeros(N);
ge_non_hat = ge_hat;
for j = 1:N
    for i = 1:N
        c1 = exp((-sx^2*(2*pi)^2*(fx(j)+kx)^2)/2);
        c2 = exp((-sx^2*(2*pi)^2*(fx(j)-kx)^2)/2);
        c3 = exp((-sx^2*(2*pi)^2*(ft(i)+kt)^2)/2);
        c4 = exp((-sx^2*(2*pi)^2*(ft(i)-kt)^2)/2);
        ge_hat(j,i) = 1/4*(c1*c3+c1*c4+c2*c3+c2*c4);
        ge_non_hat(j,i) = 1/2*(c1*c4+c2*c3);
    end
end
subplot(2,2,3)
mesh(fx,ft,ge_hat)
xlim([-10 10])
xticks([-10 0 10])
ylim([-25 25])
yticks([-20 0 20])
ylabel('temporal frequency (c/s)') 
zlim([0 0.3])
set(gca,'GridLineStyle','--')
subplot(2,2,4)
mesh(fx,ft,ge_non_hat)
xlim([-10 10])
xticks([-10 0 10])
xlabel('spatial frequency (c/deg)') 
ylim([-25 25])
yticks([-20 0 20])
zlim([0 0.6])
zlabel('response (arbitrary units)')
set(gca,'GridLineStyle','--')

figure(2)
nx = 4;
nt = 8;
t = (0:dt:2)';
gse = 1/(sqrt(2*pi)*sx)*exp(-x.^2/(2*sx^2)).*cos(2*pi*nx*x);
gso = 1/(sqrt(2*pi)*sx)*exp(-x.^2/(2*sx^2)).*sin(2*pi*nx*x);
fte = 1/(sqrt(2*pi)*st)*exp(-T.^2/(2*st^2)).*cos(2*pi*nt*T);
fto = 1/(sqrt(2*pi)*st)*exp(-T.^2/(2*st^2)).*sin(2*pi*nt*T);
gr = cos(2*pi*(nx*repmat(x,[length(t) 1]) - nt*repmat(t,[1,N])));
gr_e_s = gr*gse';
con = conv(gr_e_s,fte);
gr_e_st = con(64:64+length(t)-1)*dt;
gl = cos(2*pi*(nx*repmat(x,[length(t) 1]) + nt*repmat(t,[1,N])));
gl_e_s = gl*gse';
con = conv(gl_e_s,fte);
gl_e_st = con(64:64+length(t)-1)*dt;
subplot(1,2,1)
plot(t, gr_e_st, 'k', t, gl_e_st, 'r--')
legend('right', 'left','location', 'best')
box off 
xticks([0:0.2:2])
xlabel('time (s)')
ylabel('response (arbitrary units)')
subplot(1,2,2)
gr_o_s = gr*gso';
con = conv(gr_o_s,fto);
gr_o_st = con(64:64+length(t)-1)*dt;
gl_o_s = gl*gso';
con = conv(gl_o_s,fto);
gl_o_st = con(64:64+length(t)-1)*dt;
plot(t,gr_e_st+gr_o_st,'k',t,gl_o_st+gl_o_st,'r--')
box off
xticks([0:0.2:2])
xlabel('time (s)')
figure(3)
X = 0:1:20;
p1 = binopdf(X, 100, 0.1);
p2 = poisspdf(X, 100*0.1);
p3 = binopdf(X, 12, 0.83);
p4 = normpdf([-5:0.1:5],0,1);
subplot(2,2,1)
plot(X,p1,'-xk')
xticks([0 10 20])
ylim([0 0.2])
yticks([0:0.05:0.2])
ylabel('probability')
title('binomial')
box off
subplot(2,2,2)
plot(X,p3,'-xk')
xticks([0 10 20])
ylim([0 0.4])
yticks([0:0.1:0.4])
ylabel('probability')
title('binomial')
box off
subplot(2,2,3)
plot(X,p2,'-xk')
xticks([0 10 20])
ylim([0 0.2])
yticks([0:0.05:0.2])
ylabel('probability')
title('poisson')
box off
subplot(2,2,4)
plot([-5:0.1:5],p4,'-k')
ylim([0 0.4])
yticks([0:0.1:0.4])
ylabel('probability density')
title('normal')
box off
