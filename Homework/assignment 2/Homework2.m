v = -100:0.01:50;
figure(1)
alpha_m = 0.1*(46+v)./(1-exp(-(46+v)/10));
beta_m = 4*exp(-(v+71)/18);
tau_m = 1./(alpha_m + beta_m);
m_inf = alpha_m.* tau_m;
alpha_h = 0.07*exp(-(v+71)/20);
beta_h = 1./(exp(-(41+v)/10)+1);
tau_h = 1./(alpha_h + beta_h);
h_inf = alpha_h .* tau_h;
subplot(1,2,1)
xticks([-100 -75 -50 -25 0 25 50])
xlabel('V(mV)')
yyaxis left
plot(v,m_inf, 'r')
ylabel( 'm_\infty')
yticks([0 0.2 0.4 0.6 0.8 1])
yyaxis right
plot(v, tau_m, 'k')
ylabel( '\tau_m (ms)')
yticks([0.1 0.2 0.3 0.4 0.5 0.6])
ylim([0 0.6])
legend('\tau_m','m_\infty','location', 'best')
box off
subplot(1,2,2)
xticks([-100 -75 -50 -25 0 25 50])
xlabel('V(mV)')
yyaxis left
plot(v,h_inf, 'r')
ylabel( 'h_\infty')
yticks([0 0.2 0.4 0.6 0.8 1])
yyaxis right
plot(v, tau_h, 'k')
ylabel( '\tau_h (ms)')
yticks([0 2 4 6 8 10])
ylim([0 10])
legend('\tau_h','h_\infty','location', 'best')
box off
figure(2)
V_K = -77; %mV
V_Na = 56; %mV
V_Cl = -68; %mV
g_K_max = 36; %mS/cm^2
g_Na_max = 120; %mS/cm^2
g_Cl = 0.3; %mS/cm^2
alpha_n = 0.01*(61+v)./(1-exp(-(61+v)/10));
beta_n = exp(-(v+71)/80)/8;
tau_n = 1./(alpha_n + beta_n);
n_inf = alpha_n .* tau_n;
Iss = g_K_max * n_inf.^4.*(v-V_K) + g_Na_max* m_inf.^3.*h_inf.*(v-V_Na)+g_Cl*(v-V_Cl);
plot(v, Iss,'k')
xlabel('V(mV)')
ylabel('I_{ss} (µA/cm^2)')
xlim([-100 -60])
ylim([-10 10])
grid on
box off
figure(3)
Vr = -71; %mV
Cm = 1; %membrane capacitane, µF
r = 10*10^-4; %cell radius, cm
T = 20; % total recording time, ms
A=4*pi*r^2; %Surface area, cm*cm
I = 80*10^-6; %20ms stimulation current, µA
dt = 0.01;
alpha_m = @(v) 0.1*(46+v)./(1-exp(-(46+v)/10));
beta_m = @(v) 4*exp(-(v+71)/18);
tau_m = @(v) 1./(alpha_m + beta_m);
m_inf = @(v) alpha_m.* tau_m;
alpha_h = @(v) 0.07*exp(-(v+71)/20);
beta_h = @(v) 1./(exp(-(41+v)/10)+1);
tau_h = @(v) 1./(alpha_h + beta_h);
h_inf = @(v) alpha_h .* tau_h;
alpha_n = @(v)0.01*(61+v)./(1-exp(-(61+v)/10));
beta_n = @(v)exp(-(v+71)/80)/8;
tau_n = @(v)1./(alpha_n + beta_n);
n_inf = @(v)alpha_n .* tau_n;

N = ceil(1+T/dt); % steps
t = zeros(N,1); 
v = zeros(N,1);   
n = zeros(N,1); m = zeros(N,1); h = zeros(N,1);
I_Cl = zeros(N,1); I_K = zeros(N,1); I_Na = zeros(N,1);Istim = zeros(N,1);
t(1) = 0;
v(1) = Vr;
n(1) = alpha_n(Vr)/(alpha_n(Vr)+beta_n(Vr));  
m(1) = alpha_m(Vr)/(alpha_m(Vr)+beta_m(Vr));  
h(1) = alpha_h(Vr)/(alpha_h(Vr)+beta_h(Vr));
for j = 2:N
    t(j) = (j-1)*dt;
    Istim(j) = I*(t(j)-dt/2>2)*(t(j)-dt/2<4);
    n(j) = ((1/dt-(alpha_n(v(j-1))+beta_n(v(j-1)))/2)*n(j-1)+alpha_n(v(j-1)))/(1/dt+(alpha_n(v(j-1))+beta_n(v(j-1)))/2);
    m(j) = ((1/dt-(alpha_m(v(j-1))+beta_m(v(j-1)))/2)*m(j-1)+alpha_m(v(j-1)))/(1/dt+(alpha_m(v(j-1))+beta_m(v(j-1)))/2);
    h(j) = ((1/dt-(alpha_h(v(j-1))+beta_h(v(j-1)))/2)*h(j-1)+alpha_h(v(j-1)))/(1/dt+(alpha_h(v(j-1))+beta_h(v(j-1)))/2);
    v_half_top = 2*Cm*v(j-1)/dt+g_K_max*n(j)^4*V_K+g_Na_max*m(j)^3*h(j)*V_Na+g_Cl*V_Cl+Istim(j)/A;
    v_half_bottom = 2*Cm/dt+g_K_max*n(j)^4+g_Na_max*m(j)^3*h(j)+g_Cl;
    v(j) = 2*v_half_top/v_half_bottom - v(j-1);
    I_Cl(j) = g_Cl*(v(j)-V_Cl);
    I_K(j) = g_K_max*n(j)^4*(v(j)-V_K);
    I_Na(j) = g_Na_max*m(j)^3*h(j)*(v(j)-V_Na);
end
subplot(3,1,1)
plot(t,v,'k')
xticklabels([])
ylabel('V(mV)')
ylim([-100 100])
yticks([-100 0 100])
box off
subplot(3,1,2)
plot(t,n,'r', t, m, 'k--', t, h, 'k')
xticklabels([])
yticks([0 0.5 1])
legend('n', 'm', 'h', 'orientation', 'horizontal')
box off
subplot(3,1,3)
plot(t,I_Cl/1000,'r', t, I_K/1000, 'k--',t, I_Na/1000, 'k')
xlabel('t(ms)')
ylabel('I(mV/cm^2)')
yticks([-1 0 1])
legend('I_{Cl}', 'I_K', 'I_{Na}', 'orientation', 'horizontal')
box off
