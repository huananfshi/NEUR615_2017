figure(1)
k1 = [2 2 2 2 0.1 2 2 2];
k1n = [1 1 1 1 1 0.1 1 1];
k2 = [1.5 1.5 1.5 1.5 1.5 1.5 1.5 2];
T = 10;
dt = 0.01;
N = T/dt;
S = zeros(8,N);
E = S; ES = S; SP = S; 
S(:,1) = [8 2 8 2 8 8 8 8];
E(:,1) = [4 2 4 2 4 4 4 4];
ES(:,1) = [0 0 1 1 0 0 0 0];
SP(:,1) = [0 0 0.5 0 0 0 0.5 0];
t = linspace(0,T,N);
for j = 1:8
    for i = 2:N
    ds = k1n(j)*ES(j,i-1)-k1(j)*S(j,i-1)*E(j,i-1);
    de = k1n(j)*ES(j,i-1)+k2(j)*ES(j,i-1)-k1(j)*S(j,i-1)*E(j,i-1);
    des = k1(j)*S(j,i-1)*E(j,i-1)-k1n(j)*ES(j,i-1)-k2(j)*ES(j,i-1);
    dsp = k2(j)*ES(j,i-1);
    S(j,i) = S(j,i-1)+ds*dt;
    E(j,i) = E(j,i-1)+de*dt;
    ES(j,i) = ES(j,i-1)+des*dt;
    SP(j,i) = SP(j,i-1)+dsp*dt;
    end
end
for j = 1:8
subplot(2,4,j)
plot(t,S(j,:),'b-',t,E(j,:),'r-',t,ES(j,:),'k-',t,SP(j,:),'g-')
legend('S','E','ES','S^P','location','Best')
end
figure(2)
SP2 = zeros(8,N);
SP2(:,1) = [0 0 0.5 0 0 0 0.5 0];
Etotal = [4 2 5 3 4 4 4 4];
vmax = Etotal.*k2;
km = (k1n+k2)./k1;
for j = 1:8
   for i = 2:N
    dsp2 = vmax(j)*S(j,i-1)/(km(j)+S(j,i-1));
    SP2(j,i) = SP2(j,i-1)+dsp2*dt;
   end
end
for j = 1:8
subplot(2,4,j)
plot(t,S(j,:),'b-',t,SP2(j,:),'k-')
legend('S','S^P','location','Best')
end
figure(3)
Atotal = 1;
Ca = 0:0.01:3;
K = 30./(10+90*exp(-2*Ca));
P = 9./(10+20*exp(-20*Ca));
eta = K+P;
omega = Atotal*K./eta;
subplot(1,3,1)
plot(Ca,K,'k-',Ca,P,'r-')
legend('K(Ca)','P(Ca)')
xlabel('[Ca^{2+}]')
ylabel('Activity')
subplot(1,3,2)
plot(Ca,eta,'k-')
xlabel('[Ca^{2+}]')
ylabel('\eta')
subplot(1,3,3)
plot(Ca,omega,'k-')
xlabel('[Ca^{2+}]')
ylabel('\Omega')
figure(4)
Ap = zeros(2,N);
dAp = zeros(2,1);
for i = 2:N
dAp(1)=eta(Ca==0.5)*(omega(Ca==0.5)-Ap(1,i-1));
dAp(2)=eta(Ca==1.2)*(omega(Ca==1.2)-Ap(2,i-1));
Ap(1,i)=dAp(1)*dt+Ap(1,i-1);
Ap(2,i)=dAp(2)*dt+Ap(2,i-1);
end
F = 1/dt;
f = linspace(0,F,N);
plot(f,Ap(1,:),'k-',f,Ap(2,:),'r')
xlabel('frequency')
ylabel('LTP/LTD')
legend('[Ca^{2+}]=0.5','[Ca^{2+}]=1.2')