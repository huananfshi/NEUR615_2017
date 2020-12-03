figure(1)
x1 = [0.253 0.973];
x2 = [0.055 1.003];
x3 = [-0.099 0.015];
x4 = [0.099 -0.015];
x5 = [-0.253 -0.973];
x6 = [-0.055 -1.003];
Q = zeros(2);
Q(1,1) = (x1(1)^2+x2(1)^2+x3(1)^2+x4(1)^2+x5(1)^2+x6(1)^2)/6;
Q(2,2) = (x1(2)^2+x2(2)^2+x3(2)^2+x4(2)^2+x5(2)^2+x6(2)^2)/6;
Q(1,2) = (x1(1)*x1(2)+x2(1)*x2(2)+x3(1)*x3(2)+x4(1)*x4(2)+x5(1)*x5(2)+x6(1)*x6(2))/6;
Q(2,1) = (x1(2)*x1(1)+x2(2)*x2(1)+x3(2)*x3(1)+x4(2)*x4(1)+x5(2)*x5(1)+x6(2)*x6(1))/6;
[V,D] = eig(Q);
z1 = V(:,1);
z2 = V(:,2);
l1 = D(1,1);
l2 = D(2,2);
plot(x1(1),x1(2),'*k');
hold on 
plot(x2(1),x2(2),'*k');
plot(x3(1),x3(2),'*k');
plot(x4(1),x4(2),'*k');
plot(x5(1),x5(2),'*k');
plot(x6(1),x6(2),'*k');
plot([0 V(1,1)],[0 V(2,1)],'r')
plot([0 V(1,2)],[0 V(2,2)],'b')
hold off
figure(2)
n = 200;
x1 = normrnd(0,sqrt(1),1,n);
x2 = normrnd(0,sqrt(3),1,n);
X = [x1;x2];
k = [cos(pi/6) sin(pi/6);-sin(pi/6) cos(pi/6)];
X = k*X;
    w = zeros(2,n);
    w(:,1) = rand(2,1);
    for i = 1:n-1
        a = X(:,i);
        b = w(:,i)'*a;
        w(:,i+1) = w(:,i) + (1/200)*b*a;
    end
for m = 1:n
    plot(X(1,m),X(2,m),'kx')
    hold on
    plot(w(1,m),w(2,m),'r+')
    pause(0.05)
end 
hold off 
xlabel('x_1')
ylabel('x_2')
legend('data point', 'synaptic weights', 'location', 'best')
title('Hebb Rule \theta = 30\circ')
Q1 = zeros(2);
for j = 1:n
    Q1(1,1) = Q1(1,1)+X(1,i)^2;
    Q1(2,2) = Q1(2,2)+X(2,i)^2;
    Q1(1,2) = Q1(1,2)+X(1,i)*X(2,i);
    Q1(2,1) = Q1(2,1)+X(2,i)*X(1,i);
end
Q1 = Q1/n;
[V1,D1] = eig(Q1);
z1_h = V1(:,1);
z2_h = V1(:,2);
l1_h = D1(1,1);
l2_h = D1(2,2);
figure(3)
w1 = zeros(2,n);
w1(:,1) = w(:,1);
for j = 1:n-1
    a = X(:,j);
    b = w1(:,j)'*a;
    w1(:,j+1) = w1(:,j) + (1/(100+i))*b*(a-b*w1(:,j));
end
for q = 1:n
    plot(X(1,q),X(2,q),'kx')
    hold on 
    plot(w1(1,q),w1(2,q),'r+')
    pause(0.05)
end
hold off 
xlabel('x_1')
ylabel('x_2')
legend('data point', 'synaptic weights', 'location', 'best')
title('Oja Neuron \theta = 30\circ')
y = sum(X.*w1);
lambda = sum(y.^2)/n;