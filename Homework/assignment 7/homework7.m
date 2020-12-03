figure(1)
pr = 0.2;
N =10;
n = 100;
mu = 1;
s = 0.5;
pk = binopdf(0,N,pr);
v = 0:0.1:20;
nv = length(v);
pv = zeros(1,nv);
pv(1) = pk;
m_sn = N*pr;
m_v = N*pr*mu;
s_sn2 = N*pr*(1-pr);
s_v2 = s_sn2*mu^2+m_sn*s^2;
for i = 1:10
    pk = binopdf(i,N,pr);
    pv_k = normpdf(v,i*mu,i*sqrt(s));
    pv = pv+pk*pv_k;
end
pv(2:nv)= pv(2:nv)*(1-pv(1))/sum(pv(2:nv)*0.01);
count = zeros(1,nv);
count(1) = n*pv(1);
count(2:nv) = pv(2:nv)*(n-count(1))/(sum(pv(2:nv))*0.1);
bar(v,count)
xlabel('evoked e.p.p.s.(mV)');
ylabel('number of observations');
xlim([-0.5 15])

figure(2)
pr1 = 0.1:0.1:0.8;
np = length(pr1);
pk1 = binopdf(0,N,pr1);
pv1 = zeros(np,nv);
pv1(:,1) = pk1;

for i = 1:10
    pk1 = binopdf(i,N,pr1);
    pv1_k = normpdf(v,i*mu,i*sqrt(s));
    pv1 = pv1+pv1_k.*pk1';
end
pv1(:,2:nv)= pv1(:,2:nv).*(1-pv1(:,1))./sum(pv1(:,2:nv)*0.01);
count1 = zeros(np,nv);
count1(:,1) = n*pv1(:,1);
for j = 1:8
    count1(j,2:nv) = pv1(j,2:nv).*(n-count1(j,1))/(sum(pv1(j,2:nv))*0.1);
    subplot(2,4,j)
    bar(v,count1(j,:))
    xlabel('evoked e.p.p.s.(mV)');
    ylabel('number of observations');
    tt = sprintf('pr = %.1f',pr1(j));
    title(tt)
    xlim([-0.5 20])
end

figure(3)
m_sn = N*pr1;
m_v = N*pr1*mu;
s_sn2 = N*pr1.*(1-pr1);
s_v2 = s_sn2*mu^2+m_sn*s^2;
plot(m_v,s_v2,'*k')
hold on 
lambda = N*pr1;
m_vp = lambda*mu;
s_v2p = lambda*(mu^2+s^2);
plot(m_vp,s_v2p,'+k')
cv_minis = s/mu;
mv = 0:0.1:9;
sv2 = mu*(1+cv_minis^2)*mv-mv.^2./N;
plot(mv,sv2,'r-')
hold off
legend('Binomial','Possion','Estimate')
xlabel('Mean peak')
ylabel('Variance')
xlim([0 9])
ylim([0 11])
