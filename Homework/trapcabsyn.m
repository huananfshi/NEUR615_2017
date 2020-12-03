%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%
%  trapcabsyn.m
%
% solve the passive fiber with alpha func synaptic input via the trapezoid rule
%
% usage:        trapcabsyn(cab,stim,pinc)
%
%               cab.rad = cable radius (cm)
%               cab.ell = cable length (cm)
%    	        cab.dx = space step (cm)
%               cab.dt = timestep (ms)
%               stim.t1 = stim start time (ms)
%               stim.tau = time constant of alpha function (ms)
%               stim.Gsyn = amplitude of synaptic conductance (mS)
%               stim.loc = location of current pulse (cm)
%               stim.Tfin = stopping time (ms)
%               pinc = number of time steps between plots
%
cab = struct('rad',1e-4,'ell',.1,'dx',1e-3,'dt',0.05);
%               stim = struct('t1',1,'tau',1,'Gsyn',1e-4,'loc',0.05,'Tfin',10)
%               pinc = 4 
%
% or, for multiple stimuli
%
stim = struct('t1',[1 3],'tau',[1 1]/2,'Gsyn',1e-6*[1 1],'loc',[0.06 0.04],'Tfin',10);
%


Cm = 1;		% micro F / cm^2
g_L = 1/15; %0.3;     		% mS / cm^2
R_2 = 0.3; % 0.034;		% k Ohm cm
dx = cab.dx;
dt = cab.dt;
Nx = cab.ell/dx;                % patch length
A = 2*pi*cab.rad*dx;            % patch surface area
x = dx/2:dx:cab.ell-dx/2;       % vector of patch midpoints

Nt = ceil(stim.Tfin/dt);

v_1 = zeros(Nx,1);		% initial conditions

e = ones(Nx,1);
S1 = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
S1(1,1) = 1/dx/dx;
S1(Nx,Nx) = 1/dx/dx;

tau = Cm/g_L;
lambda = sqrt(cab.rad/(2*R_2*g_L));
A = 2*pi*cab.rad*dx;

e1 = zeros(Nx,1);
eloc = round(Nx*stim.loc/cab.ell);
vsyn = [-20 70];
Iapp = stim.Gsyn*(dt/2)/A/Cm;

B1 = (dt/2)*(speye(Nx)+lambda^2*S1)/tau;
B1 = speye(Nx) + B1;
bloc = (eloc-1)*(Nx+1)+1;
dBe = B1(bloc);

t = 0;


c0_1 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);
t = dt;
c1_1 = Iapp.*(t./stim.tau).*exp(1-t./stim.tau).*(t>stim.t1);

r1 = zeros(Nx,1);
r1(eloc) = vsyn*(c0_1 + c1_1)';

vhot = zeros(length(eloc),Nt);

for j=2:200
    
    B1(bloc) = dBe + c1_1;

    v_1 = B1\r1; 

    vhot(:,j) = v_1(eloc);

 

    t = j*dt;
    c0_1 = c1_1;
    c1_1 = Iapp.*((t-stim.t1)./stim.tau).*exp(1-(t-stim.t1)./stim.tau).*(t>stim.t1);

    r1 = 2*v_1 - r1;
    r1(eloc) = r1(eloc) + vsyn*(c0_1 + c1_1)';

end


figure(2)
t = linspace(0,stim.Tfin,Nt);
plot(t,vhot(1,:),'k') 
hold on
plot(t,vhot(2,:),'r')
hold off
box off
xlabel('t (ms)','fontsize',14)
ylabel('v (mV)','fontsize',14)

