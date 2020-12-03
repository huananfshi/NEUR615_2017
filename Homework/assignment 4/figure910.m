function [W,K,In]=figure910()

Cm = 1;		
Ra = 0.3;
cab = struct('r',10^-4,'l',0.1);
stim = struct('t1',1,'t2',2,'dose',1,'xs',0.04,'T',12);
dx = 10^-3;
x = dx/2:dx:cab.l-dx/2;	

E = struct('K', -77, 'Na', 56, 'Cl', -68);
g = struct('K', 36, 'Na', 120, 'Cl', 1/15);
lsn = 1e-4; asn = 1e-5; Ash = 1e-8;
Rss = lsn*Ra/pi/asn^2;
gamma1 = 1/Rss/Ash;
gamma2 = 1/Rss/2/pi/cab.r;
Vsyn = 20;
dt = 0.01;
Nx = cab.l/dx;		
Nt = stim.T/dt+1;



KA = struct('kp',1.1,'kn',0.19,'gmax',200);
KN = struct('kp',0.072,'kn',0.0066,'gmax',100);

c = round(Nx*stim.xs/cab.l);

e = ones(Nx,1);
B = spdiags([-e 2*e -e], -1:1, Nx, Nx)/dx/dx;
B(1,1) = 1/dx/dx;
B(Nx,Nx) = 1/dx/dx;
B = (cab.r/2/Ra)*B;
dB = diag(B);



co =  gamma2*g.Cl/(g.Cl+gamma1)/dx;
options = optimset('Jacobian','on');
V = fsolve(@(V) Issp(V,E,g,B,co,c),-70*e,options);    

W = zeros(Nt,1); 
W(1) = (g.Cl*E.Cl + gamma1*V(c))/(g.Cl+gamma1);

n = an(V)./(an(V)+bn(V)); 
m = am(V)./(am(V)+bm(V)); 
h = ah(V)./(ah(V)+bh(V)); 

Ra = 0;
Rn = 0;

t = 0;
Ia = zeros(Nt,1);
In = zeros(Nt,1);
K = zeros(Nt,1);
for j=2:Nt
      t = (j-1)*dt;
      T = stim.dose.*(t>stim.t1)*(t<stim.t2);
      Ra = ( Ra + dt*KA.kp*T ) / (1 + dt*(KA.kp*T + KA.kn) );
      Rn = ( Rn + dt*KN.kp*T ) / (1 + dt*(KN.kp*T + KN.kn) );
      ga = KA.gmax*Ra;
      gn= KN.gmax*Rn;

      tmp = ga + gn*W(j-1)/(1+2*exp(-0.062*W(j-1))/3.57);

      Wnum = (Cm/dt)*W(j-1) + g.Cl*E.Cl + tmp*Vsyn + gamma1*V(c);
      Wden = (Cm/dt) + g.Cl + tmp + gamma1;
      W(j) = Wnum/Wden;

      a = an(V);  
      n = (n + dt*a) ./ (1 + dt*(a+bn(V)) ); n4 = n.^4;

      a = am(V);
      m = (m + dt*a) ./ (1 + dt*(a+bm(V)) );

      a = ah(V);  
      h = (h + dt*a) ./ (1 + dt*(a+bh(V)) ); m3h = m.^3.*h;

      d = g.Na.*m3h + g.K.*n4 + g.Cl;
      d(c) = d(c)  + gamma2/dx;

      f = g.Na.*m3h*E.Na + g.K.*n4*E.K + g.Cl*E.Cl;
      f(c) = f(c)  + gamma2*W(j)/dx;

      B(1:Nx+1:end) = dB + d + Cm/dt;         % update the diagonal

      V = B\(Cm*V/dt + f);

      
      K(j) = KN.gmax/(1+2*exp(-0.062*W(j))/3.57);
      In(j) = KN.gmax*(W(j)-Vsyn)/(1+2*exp(-0.062*W(j))/3.57);

     

end


return

function [val, jac] = Issp(V,E,g,B,co,eloc)
Nx = length(V);
a = am(V); b = bm(V);
m = a./(a+b);
dm = (dam(V).*b - a.*dbm(V))./(a+b).^2;

a = ah(V); b = bh(V);
h = a./(a+b);
dh = (dah(V).*b - a.*dbh(V))./(a+b).^2;

a = an(V); b = bn(V);
n = a./(a+b);
dn = (dan(V).*b - a.*dbn(V))./(a+b).^2;

m3h = m.^3.*h;
n4 = n.^4;

val = B*V + g.Na.*m3h.*(V-E.Na) + g.K.*n4.*(V-E.K) + g.Cl.*(V-E.Cl);
val(eloc) = val(eloc) - co*(E.Cl-V(eloc));

dj = g.Na.*((3*dm.*m.^2.*h + m.^3.*dh).*(V-E.Na) + m3h) + ...
     g.K.*(4*dn.*n.^3.*(V-E.K) + n4) + g.Cl;
jac = B + spdiags(dj,0,Nx,Nx);
jac(eloc,eloc) = jac(eloc,eloc) + co;

function val = an(v)
val = .01*(10-(v+71))./(exp(1-(v+71)/10)-1);

function val = dan(v)
tmp = exp(-(61+v)/10);
val = -( tmp.*(71+v) - 10 )./(tmp-1).^2/1000;

function val = bn(v)
val = .125*exp(-(v+71)/80);

function val = dbn(v)
val = -exp(-(v+71)/80)/640;

function val = am(v)
val = .1*(25-(v+71))./(exp(2.5-(v+71)/10)-1);

function val = dam(v)
tmp = exp(-(46+v)/10);
val = -( tmp.*(56+v) - 10 )./(tmp-1).^2/100;

function val = bm(v)
val = 4*exp(-(v+71)/18);

function val = dbm(v)
val = -(2/9)*exp(-(v+71)/18);

function val = ah(v)
val = 0.07*exp(-(v+71)/20);

function val = dah(v)
val = -(7/2000)*exp(-(v+71)/20);

function val = bh(v)
val = 1./(exp(3-(v+71)/10)+1);

function val = dbh(v)
tmp = exp(-(v+41)/10);
val = tmp./(tmp+1).^2/10;