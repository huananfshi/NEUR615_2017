figure(1)
syn = struct('t1',1,'gsyn',1*1e-4,'loc',0.06,'vsyn',-20);
[~, v1_1] = mid1(syn);
syn = struct('t1',3,'gsyn',1*1e-4,'loc',0.04,'vsyn',70);
[t, v1_2] = mid1(syn);
v1(1,:) = v1_1;
v1(2,:) = v1_2;
prt(t,v1)
figure(2)
syn = struct('t1',3,'gsyn',1*1e-4,'loc',0.04,'vsyn',-20);
[~, v1_1] = mid1(syn);
syn = struct('t1',1,'gsyn',1*1e-4,'loc',0.06,'vsyn',70);
[t, v1_2] = mid1(syn);
v1(1,:) = v1_1;
v1(2,:) = v1_2;
prt(t,v1)

function [t,v1] = mid1(syn)
l = 0.1; %cm, cable length
a = 1*10^-4; %cm, canle radius 
Cm = 1; %µF/cm^2
T = 10; %ms
gCl = 1/15; %mS/cm^2
Ra = 0.3; %kOhm cm
dx = 1*10^-4; %cm
dt = 0.05;  
tau = Cm/gCl;
tau_alpha = 1/2; %ms
lambda = sqrt(a/(2*Ra*gCl));
N_x = l/dx;
A = 2*pi*a*dx;
x = dx/2:dx:l-dx/2;
N_t = T/dt;
v = zeros(N_x,1);
e = ones(N_x,1);
S = spdiags([-e 2*e -e], -1:1, N_x, N_x)/dx/dx;
S(1,1) = 1/dx/dx;
S(N_x,N_x) = 1/dx/dx;
B = dt/2*(lambda^2*S+speye(N_x))/tau;
B = speye(N_x)+B;
c = round(syn.loc/dx);
d = (c-1)*(N_x+1)+1;
Ba= B(d);
I = syn.gsyn*(dt/2)/Cm/A;
v1 = zeros(length(syn.loc),N_t);
t = 0;
c0 = I.*((t-syn.t1)./tau_alpha).*exp(1-(t-syn.t1)./tau_alpha).*(t>syn.t1);
t = dt;
c1 = I.*((t-syn.t1)./tau_alpha).*exp(1-(t-syn.t1)./tau_alpha).*(t>syn.t1);
r = zeros(N_x,1);
r(c) = syn.vsyn*(c0 + c1)';


for j=2:N_t
    
    B(d) = Ba+c1;
    v = B\r;
    v1(:,j) = v(c);

    t = j*dt;

    c0 = c1;
    c1 = I.*((t-syn.t1)./tau_alpha).*exp(1-(t-syn.t1)./tau_alpha).*(t>syn.t1);
    r = 2*v - r;
    r(c) = r(c) + syn.vsyn*(c0 + c1)';

end
t = linspace(0,T,N_t);

end

function prt(t,v1)
plot(t,v1(1,:),'k',t,v1(2,:),'r')
box off
xlabel('t (ms)')
ylabel('v (mV)')
xticks([0 1 2 3 4 5 6 7 8 9 10])
legend('inhibitory','excitatory','location', 'best')
peak1 = min(v1(1,:));
loc1 = find(v1(1,:) == peak1);
peak2 = max(v1(2,:));
loc2 = find(v1(2,:) == peak2);
ploc = [loc1 loc2];
tp = t(ploc);
X = sprintf('For inhibitory, the peak value is %f mV, the time of peak is at %f ms.',peak1,tp(1));
Y = sprintf('For excitatory, the peak value is %f mV, the time of peak is at %f ms.',peak2,tp(2));
display(X)
display(Y)
end