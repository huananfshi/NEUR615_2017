dt = 0.01; %ms
T = 35; %ms
g_syn_max_1 = 0.2;  %synaptic maximum conductance mS/cm^2, excitatory AMAP only
tau_alpha_1 = 2; %ms
t1_1 = 5; %ms 
V_syn_1 = 0; % mV
g_syn_max_2 = [0.2 0.2]; %synaptic maximum conductance mS/cm^2, paired excitatory AMAP and inhibitory GABA
tau_alpha_2 = [2 2]; %ms
t1_2 = [5 4]; %ms 
V_syn_2 = [0 -68]; % mV
[t, g_syn_1, v_1] = trapezoid(dt, T, g_syn_max_1, tau_alpha_1, t1_1, V_syn_1); % call trapezoid
[t, g_syn_2, v_2] = trapezoid(dt, T, g_syn_max_2, tau_alpha_2, t1_2, V_syn_2);
subplot(1,2,1)
plot(t,g_syn_1,'k')
hold on
plot(t,g_syn_2(:,1),'r--', t, g_syn_2(:,2), 'r')
hold off
box off
legend('excitatory solo','excitatory paired','inhibitory paired','location','best')
xlabel('t(ms)')
ylabel('g_{syn}(µS/cm^2)')
xlim([0 35])
xticks([0 5 10 15 20 25 30 35])
subplot(1,2,2)
plot(t,v_1,'k')
hold on 
plot(t,v_2,'r')
hold off
box off
legend('excitatory only','inhibitory and excitatory','location','best')
xlabel('t(ms)')
ylabel('V(mV)')
xlim([0 35])
xticks([0 5 10 15 20 25 30 35])
