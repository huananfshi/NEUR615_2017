VCl = -68; %Cl channel Nernst potential, mV
r = 10*10^-4; %cell radius, cm
gCl = 0.3; %Cl channel conductance mS/cm^2
Cm = 1; %membrane capacitane, µF
I = 10*10^-6; %20ms stimulation current, µA
tau = Cm/gCl; %membrane time constant, ms
A=4*pi*r^2; %Surface area, cm*cm
T = 40; % total recording time, ms
dt = 0.01; %ms
N = ceil(1+T/dt); % steps
v = zeros(N,1);  v(1)=VCl; %define v(t)
t = zeros(N,1); Istim = zeros(N,1);    % define t and Istim
Ic = zeros(N,1); ICl = zeros(N,1); % define Ic and ICl
for j = 2:N
    t(j)=j*dt;
    
    Istim(j) = (t(j)>2)*(t(j)<22)*I;%from 2ms to 22ms, apply 10pA stimulation current
    v(j) = (v(j-1) + dt*(VCl/tau + Istim(j)/(A*Cm)))/(1+dt/tau); 
    %backward Euler method, f = V'(t) = VCl/tau + Istim(j)/(A*Cm)
    ICl(j) = (v(j)-VCl)*gCl;  %Ohm's law 
    Ic(j) = (Istim(j)-A*ICl(j))/A;  %Kirchhoff?s Current Law

    
   
end
subplot(1,2,1);  %plot v - t
plot(t,v)
xlabel('t (ms)')
ylabel('V (mV)')
xlim([0 40]) 
xticks([0 5 10 15 20 25 30 35 40])
subplot(1,2,2); %plot Ic & ICl - t
plot(t,Ic)
hold on 
plot(t,ICl)
hold off
legend('Ic','ICl')
xlabel('t (ms)')
ylabel('Membrane Current Density (µA/cm^2)')
xlim([0 40])
xticks([0 5 10 15 20 25 30 35 40])
