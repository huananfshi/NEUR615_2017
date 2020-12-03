%learning process, equal probability
figure(1)
x1 = [1 0.1; 0.2 0.9];
w0 = [0.1;0.1];
p = [1/2 1/2];
tau = 5;
eta = 0.1;
y = w0'*x1;
theta = p*y'.^2;
w = w0;
plot([0 x1(1,1)],[0 x1(2,1)],'-r')
hold on 
plot([0 x1(1,2)],[0 x1(2,2)],'-b')
plot([0 w(1)],[0 w(2)],'k')
xlabel('x_1')
xlim([-1 3])
ylabel('x_2')
ylim([-0.5 1])
for i = 1:500
    dtheta = 1/tau*(p*y'.^2-theta);
    theta = dtheta+theta;
    dw = eta*x1(:,1)*p(1)*y(1)*(y(1)-theta)+eta*x1(:,2)*p(2)*y(2)*(y(2)-theta);
    w = dw+w;
    y = w'*x1;
    if (mod(i,20) == 0)
    plot([0 w(1)],[0 w(2)],'k')
    pause(0.2)
    end
end
hold off
box off
legend('x_1','x_2','w')
title('equal probability')
%Fixed point, equal probability
figure(2)
w1 = bcm(x1, p);
plotbcm(x1,w1)
title('equal probability')
%Leaning process, different probabilities
figure(3)
p = [1/4 3/4];
y = w0'*x1;
theta = p*y'.^2;
w = w0;
plot([0 x1(1,1)],[0 x1(2,1)],'-r')
hold on 
plot([0 x1(1,2)],[0 x1(2,2)],'-b')
plot([0 w(1)],[0 w(2)],'k')
xlabel('x_1')
xlim([-1 2])
ylabel('x_2')
ylim([-0.5 2])
for i = 1:500
    dtheta = 1/tau*(p*y'.^2-theta);
    theta = dtheta+theta;
    dw = eta*x1(:,1)*p(1)*y(1)*(y(1)-theta)+eta*x1(:,2)*p(2)*y(2)*(y(2)-theta);
    w = dw+w;
    y = w'*x1;
    if (mod(i,20) == 0)
    plot([0 w(1)],[0 w(2)],'k')
    pause(0.2)
    end
end
hold off
box off
legend('x_1','x_2','w')
title('different probabilities')
%Fixed point, different probabilities
figure(4)
w2 = bcm(x1,p);
plotbcm(x1,w2)
title('different probabilities')
%Learning process, 4 inputs, equal probability 
x2 = [1 0.1 0.9 0.2; 0.2 0.9 0.3 1; 1 0.2 1.1 0.3; 0.2 1 0.1 1.1];
p = [1/4 1/4 1/4 1/4];
w = [0.1 0.2 0.3 0.4]';
y = w'*x2;
theta = p*y'.^2;
for i = 1:5000000
    dtheta = 1/tau*(p*y'.^2-theta);
    theta = dtheta+theta;
    dw = x2(:,1)*p(1)*y(1)*(y(1)-theta)+x2(:,2)*p(2)*y(2)*(y(2)-theta);
    dw = dw +x2(:,3)*p(2)*y(3)*(y(3)-theta)+x2(:,4)*p(4)*y(4)*(y(4)-theta);
    w = eta*dw+w;
    y = w'*x2;
end

%fixed point, 4 inputs, equal probability
x2 = [1 0.1 0.9 0.2; 0.2 0.9 0.3 1; 1 0.2 1.1 0.3; 0.2 1 0.1 1.1];
p3 = [1/4 1/4 1/4 1/4];
y1 = zeros(4,16);
y1(:,1) = [0 0 0 0]';
y1(:,2) = [1/p3(1) 0 0 0]';
y1(:,3) = [0 1/p3(2) 0 0]';
y1(:,4) = [0 0 1/p3(3) 0]';
y1(:,5) = [0 0 0 1/p3(4)]';
y1(:,6) = [1/(p3(1)+p3(2)) 1/(p3(1)+p3(2)) 0 0]';
y1(:,7) = [1/(p3(1)+p3(3)) 0 1/(p3(1)+p3(3)) 0]';
y1(:,8) = [1/(p3(1)+p3(4)) 0 0 1/(p3(1)+p3(4))]';
y1(:,9) = [0 1/(p3(3)+p3(2)) 1/(p3(3)+p3(2)) 0]';
y1(:,10) = [0 1/(p3(4)+p3(2)) 0 1/(p3(4)+p3(2))]';
y1(:,11) = [0 0 1/(p3(3)+p3(4)) 1/(p3(3)+p3(4))]';
y1(:,12) = [1/(p3(1)+p3(2)+p3(3)) 1/(p3(1)+p3(2)+p3(3)) 1/(p3(1)+p3(2)+p3(3)) 0]';
y1(:,13) = [1/(p3(1)+p3(2)+p3(4)) 1/(p3(1)+p3(2)+p3(4)) 0 1/(p3(1)+p3(2)+p3(4))]';
y1(:,14) = [1/(p3(1)+p3(3)+p3(4)) 0 1/(p3(1)+p3(3)+p3(4)) 1/(p3(1)+p3(3)+p3(4))]';
y1(:,15) = [0 1/(p3(2)+p3(3)+p3(4)) 1/(p3(2)+p3(3)+p3(4)) 1/(p3(2)+p3(3)+p3(4))]';
y1(:,16) = [1/sum(p3) 1/sum(p3) 1/sum(p3) 1/sum(p3)]';
w1 = zeros(4,16);
for i = 1:16
    w1(:,i) = x2\y1(:,i);
end

function w = bcm(x, p)
y = zeros(2,4);
w = y;
y(:,1) = [0;0];
y(:,2) = [1/p(1);0];
y(:,3) = [0;1/p(2)];
y(:,4) = [1/(p(1)+p(2));1/(p(1)+p(2))];
for i = 1:4
    w(:,i) = x\y(:,i);
end
end

function plotbcm(x,w)
plot([0 x(1,1)],[0 x(2,1)],'-r')
hold on 
plot([0 x(1,2)],[0 x(2,2)],'-b')
for i = 1:4
plot([0 w(1,i)],[0 w(2,i)], '-k')
end
hold off
legend('x_1','x_2', 'w')
xlabel('x_1')
ylabel('x_2')
box off 
end 