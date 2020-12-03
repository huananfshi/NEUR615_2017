g_syn_max = 1*10^-6; %mS
tau_alpha = 0.5; %ms
tk = [1 3]; %ms
v_syn = [70 70]; %mV
T = 10; %ms
a = 1e-4*[1 1 1];
l = [250 250 250]*10^-4; %cm
dx = 1*10^-4; %cm
dt = 0.05; %ms
Nx = l/dx;
Nt = T/dt;
A = 2*pi*a(3)*dx;
As = 400*pi*10^-8;%cm^2
rho = A/As;
Ra = 0.3; %kOhm cm
gCl = 1/15; %mS/cm^2
Cm = 1; %µF/cm^2
tau = Cm/gCl;
loc = [100*10^-4, 100*10^-4]; %cm
lammda2 = a/(2*Ra*gCl)/dx^2; %lambda^2/dx^2  
r = a/a(3);
H1 = [-2*lammda2(1)*ones(1,Nx(1)) -2*lammda2(2)*ones(1,Nx(2)) -2*lammda2(3)*ones(1,Nx(3)+1)];
H1(1) = -lammda2(1);
H1(Nx(1)+1) = -lammda2(2);
H1(Nx(1)+Nx(2)+1) =  -lammda2*r';
H1(end) = -rho*lammda2(3);
N = length(H1);

H2 = [lammda2(1)*ones(1,Nx(1)-1) 0 lammda2(2)*ones(1,Nx(2)) lammda2(3)*ones(1,Nx(3))];
H3 = [lammda2(1)*ones(1,Nx(1)-1) 0 lammda2(2)*ones(1,Nx(2)-1) r(2)*lammda2(2) lammda2(3)*ones(1,Nx(3))];
H3(end) = rho*H3(end);

H = spdiags( [[H3 0]' H1' [0 H2]'], -1:1, N, N);

H(Nx(1)+Nx(2)+1,Nx(1)) = r(1)*lammda2(1);
H(Nx(1),Nx(1)+Nx(2)+1) = lammda2(1);

I = speye(N);

B = (H-I)/tau;
B = I-dt*B/2;
c = round(N*loc/sum(l));
Ba = B(c,c);
x3 = 0:dx:l(3);
x1 = l(3):dx:l(3)+l(1);
x2 = l(3):dx:l(3)+l(2);
v = zeros(N,1);         
r1 = zeros(N,1);
v1 = zeros(2,Nt);
t = 0;
g_syn = g_syn_max*(dt/2)/A/Cm;
c_syn_1 = g_syn.*((t-tk)./tau_alpha).*exp(1-(t-tk)./tau_alpha).*(t>tk);
t = dt;
c_syn_2 = g_syn.*((t-tk)./tau_alpha).*exp(1-(t-tk)./tau_alpha).*(t>tk);

r = zeros(N,1);
r(c) = (v_syn.*sum(c_syn_1 + c_syn_2))';

for j=2:Nt

    B(c,c) = Ba + sum(c_syn_2);

    v = B\r;

    v1(:,j) = v(c);

    t = j*dt;
    c_syn_1 = c_syn_2;
    c_syn_2 = g_syn.*((t-tk)./tau_alpha).*exp(1-(t-tk)./tau_alpha).*(t>tk);
    r = 2*v - r;
    r(c) = r(c) + (v_syn.*sum(c_syn_1 + c_syn_2))';

end
t = linspace(0,T,Nt);
figure(1)
   plot(t,v1,'k')
   xlabel('t(ms)')
   ylabel('V(mV)')
figure(2)
[W,K,In]=figure910;
subplot(2,1,1)
plot(W(W~=W(1)),K(W~=W(1)),'k')
ylabel('$\bar g_N \cdot M(mS/cm^2)$','interpreter','latex')
ylim([0 100])
xlim([-100 50])
xticks([-100 -50 0 50])
yticks([0 50 100])
box off
subplot(2,1,2)
plot(W(W~=W(1)),In(W~=W(1)),'k')
ylabel('current(µA/cm^2')
xlabel('membrane potential (mV)')
ylim([-1500 0])
xlim([-100 50])
xticks([-100 -50 0 50])
yticks([-1500 -1000 -500 0])
box off
