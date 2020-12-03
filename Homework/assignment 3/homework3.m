l = 0.1; %cm, cable length
a = 1*10^-4; %cm, canle radius 
Cm = 1; %µF/cm^2
T = 10; %ms
I = 100*10^-6; %µA
gCl = 1/15; %mS/cm^2
Ra = 0.3; %kOhm cm
dx = 1*10^-4; %cm
dt = 0.05;
loc = [0.06 0.04]; %cm, stimulation location  
t1 = [1 2];%ms,   stimulation start-time 
t2 = [2 3];%ms,   stimulation end-time 
tau = Cm/gCl;
lambda = sqrt(a/(2*Ra*gCl));
figure(1)
N_x = l/dx;
A = 2*pi*a*dx;
x = dx/2:dx:l-dx/2;
N_t = T/dt;
v = zeros(N_x,1);
e = ones(N_x,1);
S = spdiags([e -2*e e], -1:1, N_x, N_x)/dx/dx;
S(1,1) = 1/dx/dx;
S(N_x,N_x) = 1/dx/dx;
B = (lambda^2*S-speye(N_x))/tau;
[L,U] = lu(speye(N_x)-(dt/2)*B);
I = I/Cm/A;
v1 = zeros(2,N_t);
subplot(1,2,1)
ts = zeros(N_t,1);
t = 0;
f0 = I*(t>t1).*(t<t2);
plot3(x,t*e,v)
hold on
t = dt;
f1 = I*(t>t1).*(t<t2);
c = round(loc/dx);
r = zeros(N_x,1);
r(c) = dt/2*(f0 + f1)';


for j=2:200
    
    v = U\(L\r);

    v1(:,j) = v(c);

    if mod(j,4) == 0
        plot3(x,t*e,v, 'k')
    end

    t = j*dt;
    ts(j)=j*dt;

    f0 = f1;
    f1 = I.*(t>t1).*(t<t2);
    r = 2*v - r;
    r(c) = r(c) + dt/2*(f0 + f1)';

end
xlabel('x(cm)')
ylabel('t(ms)')
zlabel('v(mV)')
hold off
xticks([0 0.02 0.04 0.06 0.08 0.1])
yticks([0 2 4 6 8 10])
zlim([0 10])
zticks([0 2 4 6 8 10])
subplot(1,2,2)
plot(ts,v1(1,:),'k',ts,v1(2,:),'r')
box off
xlabel('t(ms)')
ylabel('v(mV)')
xticks([0 1 2 3 4 5 6 7 8 9 10])
yticks([0 1 2 3 4 5 6 7 8 9])
ylim([0 9])

figure(2)
Rin =  1/(2*pi*a*lambda*gCl)*cosh(x/lambda).*cosh((l-x)/lambda)./sinh(l/lambda);
plot(x,Rin/1000,'k')
box off
xlabel('x_s(cm)')
ylabel('R_{in}(M\Omega)')
ylim([300 500])

figure(3)
I = 1*10^-3; %µA
v_inf = I*Ra*lambda*cosh((l-x)/lambda)./(pi*a^2*sinh(l/lambda)); %mV
plot(x,v_inf,'k')
xlabel('x(cm)')
ylabel('V(mV)')
