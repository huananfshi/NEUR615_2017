figure(1)
m = [30 15];
s2 = 1.5*m;
s = sqrt(s2);
x = 0:0.01:60;
pp = normpdf(x,m(1),s(1));
pn = normpdf(x,m(2),s(2));
plot(x,pp,'k',x,pn,'r')
xlabel('spikes')
legend('preferred stimulus','non-preferred stimulus')
figure(2)
mu = (m(1)-m(2))/s(2);
sigma = s(1)/s(2);
rho_0 = 2*sigma^2*log(sigma);
l_val_vect = -3:0.05:3;
for j = 1:length(l_val_vect)
        l_val = l_val_vect(j);
        rho = rho_0 + 2*sigma^2*l_val;

        d_2 = sigma^2*(mu^2+rho)-rho;
    
        if ( d_2 > 0 )
            d = sqrt(sigma^2*(mu^2+rho)-rho);

            y_p = (-mu + d)/(sigma^2-1);
            y_m = (-mu - d)/(sigma^2-1);

            pfa(j) = normcdf(y_m,0,1) + (1-normcdf(y_p,0,1));
            pd(j) = normcdf(y_m,mu,sigma) + (1-normcdf(y_p,mu,sigma));
        end
end
inds1 = find(~isnan(pfa));
inds2 = find(~isnan(pd));
[p_fa, ix] = sort(pfa(inds1));
p_fa = [0 p_fa 1];
p_d = pd(inds1(ix));
p_d = [0 p_d 1];
plot(p_fa,p_d,'k');
xlabel('P_{FA}')
ylabel('P_D')
xlim([0 1])
ylim([0 1])
xticks([0 0.5 1])
yticks([0 0.5 1])
figure(3)
error = 0.3*p_fa + 0.7*(1-p_d);
plot(p_fa,error,'k')
xlabel('P_{FA}')
ylabel('error rate')
xticks([0 0.5 1])
yticks([0:0.2:0.8])
min_error = min(error);
min_pfa = pfa(find(error == min_error));
x = sprintf('The minumum error rate is %f. and the corresponding fals-alarm rate is %f.',min_error,min_pfa);
display(x)
area = trapz(p_fa,p_d);
x = sprintf('The area under the ROC curve is %f.',area);
display(x)