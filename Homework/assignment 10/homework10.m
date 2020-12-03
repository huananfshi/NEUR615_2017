rand('state',sum(100*clock));
figure(1)
input = [0 0;0 1;1 0;1 1];
out_and = [0;0;0;1];
out_or = [0;1;1;1];
out_xor = [0;1;1;0];
w = [1 1];
p = rand(10,2);
subplot(1,3,1)
ppt(w,input,out_and,p,1)
title('Perceptron learning rule for AND gate')
box on
subplot(1,3,2)
ppt(w,input,out_or,p,1)
title('Perceptron learning rule for OR gate')
box on
subplot(1,3,3)
ppt(w,input,out_xor,p,2)
title('Perceptron learning rule for XOR gate')
box on 
figure(2)
bias = -1;
subplot(2,3,1)
x0 = rand(1,10000);
w = [1 1];
for i = 1:10000
    y0 = 0.5+2*x0(1,i);
    o0 = w(2)*x0(1,i)+bias*w(1);
    dw = (y0-o0)*x0(1,i);
    w(2) = w(2)+dw;
    w(1) = w(1)+bias*(y0-o0);
end
p = rand(1,20);
o = w(2)*p+bias*w(1);

plot(p,o,'k+')
hold on 
x  = linspace(0,1,20);
y = 0.5+2*x;
plot(x,y,'k')
hold off
title('Linear Network for linear function')
subplot(2,3,4)
w = [1 1];
for i = 1:10000
    y0 = 0.5+2*x0(1,i);
    o0 = w(2)*x0(1,i)+bias*w(1);
    o0 = 1/(1+exp(-o0));
    dw = (y0-o0)*x0(1,i);
    w(2) = w(2)+dw;
    w(1) = w(1)+bias*(y0-o0);
end
p = rand(1,20);
o = w(2)*p+bias*w(1);
o = 1./(1+exp(-o));
plot(p,o,'k+')
hold on 
x  = linspace(0,1,20);
y = 0.5+2*x;
plot(x,y,'k')
hold off
title('Sigmoidal Network for linear function')
subplot(2,3,2)
w = [1 1];
for i = 1:10000
    y0 = 1/(1+exp(-4*(x0(1,i)-0.2)));
    o0 = w(2)*x0(1,i)+bias*w(1);
    dw = (y0-o0)*x0(1,i);
    w(2) = w(2)+dw;
    w(1) = w(1)+bias*(y0-o0);
end
p = rand(1,20);
o = w(2)*p+bias*w(1);
plot(p,o,'k+')
hold on 
y = 1./(1+exp(-4*(x-0.2)));
plot(x,y,'k')
hold off
title('Linear Network for sigmoidal function')
subplot(2,3,5)
w = [1 1];
for j = 1:10000
    for i = 1:10000
    y0 = 1/(1+exp(-4*(x0(1,i)-0.2)));
    o0 = w(2)*x0(1,i)+bias*w(1);
    os = 1/(1+exp(-o0));
    dw = 2*0.7*(y0-o0)*x0(1,i)*os;
    w(2) = w(2)+dw;
    w(1) = w(1)+bias*(y0-o0);
    end
end
p = rand(1,20);
o = w(2)*p+bias*w(1);
o = 1./(1+exp(-o));
plot(p,o,'k+')
hold on 
y = 1./(1+exp(-4*(x-0.2)));
plot(x,y,'k')
hold off
title('Sigmoidal Network for sigmoidal function')
subplot(2,3,3)
w = [1 1];
for i = 1:10000
    y0 = (x0(1,i)-0.5)^2;
    o0 = w(2)*x0(1,i)+bias*w(1);
    dw = (y0-o0)*x0(1,i);
    w(2) = w(2)+dw;
    w(1) = w(1)+bias*(y0-o0);
end
p = rand(1,20);
o = w(2)*p+bias*w(1);
plot(p,o,'k+')
hold on 
y = (x-0.5).^2;
plot(x,y,'k')
hold off
title('Linear Network for quadratic function')
subplot(2,3,6)
w = [1 1];
for i = 1:10000
    y0 = (x0(1,i)-0.5)^2;
    o0 = w(2)*x0(1,i)+bias*w(1);
    o0 = 1/(1+exp(-o0));
    dw = (y0-o0)*x0(1,i);
   w(2) = w(2)+dw;
    w(1) = w(1)+bias*(y0-o0);
end
p = rand(1,20);
o = w(2)*p+bias(1);
o = 1./(1+exp(-o));
plot(p,o,'k+')
hold on 
y = (x-0.5).^2;
plot(x,y,'k')
hold off
title('Sigmoidal Network for quadratic function')
figure(3)
subplot(1,2,1)
eta = 0.7;
train_time = 10000;
w3 = -1 +2.*rand(3,3);
bias = [-1 -1 -1];
for i = 1:train_time
   out = zeros(4,1);
   N = length (input(:,1));
   for j = 1:N
      x1 = bias(1,1)*w3(1,1)+ input(j,1)*w3(1,2)+ input(j,2)*w3(1,3);
      x2 = bias(1,2)*w3(2,1)+ input(j,1)*w3(2,2)+ input(j,2)*w3(2,3);
      y2(1) = s(x1);
      y2(2) = s(x2);
      y3 = bias(1,3)*w3(3,1)+ y2(1)*w3(3,2)+ y2(2)*w3(3,3);
      out(j) = s(y3);
     
      delta1 = out(j)*(1-out(j))*(out_xor(j)-out(j));
      
      delta2_1 = y2(1)*(1-y2(1))*w3(3,2)*delta1;
      delta2_2 = y2(2)*(1-y2(2))*w3(3,3)*delta1;
      
      for k = 1:3
         if k == 1 
            w3(1,k) = w3(1,k) + eta*bias(1,1)*delta2_1;
            w3(2,k) = w3(2,k) + eta*bias(1,2)*delta2_2;
            w3(3,k) = w3(3,k) + eta*bias(1,3)*delta1;
         else 
            w3(1,k) = w3(1,k) + eta*input(j,1)*delta2_1;
            w3(2,k) = w3(2,k) + eta*input(j,2)*delta2_2;
            w3(3,k) = w3(3,k) + eta*y2(k-1)*delta1;
         end
      end
   end   
end
p = rand(10,2);
for q = 1:10
x1 = bias(1,1)*w3(1,1)+ p(q,1)*w3(1,2)+ p(q,2)*w3(1,3);
x2 = bias(1,2)*w3(2,1)+ p(q,1)*w3(2,2)+ p(q,2)*w3(2,3);
y2(1) = s(x1);
y2(2) = s(x2);
y3 = bias(1,3)*w3(3,1)+ y2(1)*w3(3,2)+ y2(2)*w3(3,3);
outq(q) = s(y3);
hold on
if round(outq(q))==1
   plot(p(q,1),p(q,2),'rx');
else if round(outq(q))==0
   plot(p(q,1),p(q,2),'b+');
    end
end
hold off
title('XOR gate learning test for multi-layer neurons')
box on 
end
subplot(1,2,2)
eta = 0.7;
train_time = 10000;
w3 = -1 +2.*rand(7,4);
x0 = linspace(0,1,200);
input_quad(:,1)=x0;
input_quad(:,2:3)=zeros(200,2);
bias = [-1 -1 -1 -1 -1 -1 -1];
y0 = (x0-0.5).^2;
for i = 1:train_time
   out = zeros(200,1);
   N = length (input_quad(:,1));
   for j = 1:N
      x1 = bias(1,1)*w3(1,1)+ input_quad(j,1)*w3(1,2)+ input_quad(j,2)*w3(1,3);
      x2 = bias(1,2)*w3(2,1)+ input_quad(j,1)*w3(2,2)+ input_quad(j,2)*w3(2,3);
      x3 = bias(1,3)*w3(3,1)+ input_quad(j,1)*w3(3,2)+ input_quad(j,2)*w3(3,3);
      y2(1) = s(x1);
      y2(2) = s(x2);
      y2(3) = s(x3);
      x4 = bias(1,4)*w3(4,1)+ y2(1)*w3(4,2)+ y2(2)*w3(4,3);
      y3(1) = s(x4);
      x5 = bias(1,5)*w3(5,1)+ y2(1)*w3(5,2)+ y2(3)*w3(5,3);
      y3(2) = s(x5);
      x6 = bias(1,6)*w3(6,1)+ y2(2)*w3(5,2)+ y2(3)*w3(6,3);
      y3(3) = s(x6);
      y4 = bias(1,7)*w3(7,1)+ y3(1)*w3(7,2)+ y3(2)*w3(7,3)+y3(3)*w3(7,4);
      out(j)=s(y4);
      delta1 = out(j)*(1-out(j))*(y0(j)-out(j));
      
      delta2_1 = y3(1)*(1-y3(1))*w3(7,2)*delta1;
      delta2_2 = y3(2)*(1-y3(2))*w3(7,3)*delta1;
      delta2_3 = y3(2)*(1-y3(2))*w3(7,4)*delta1;
      delta3_1 = y2(1)*(1-y2(1))*w3(4,2)*delta2_1...
          +y2(1)*(1-y2(1))*w3(5,2)*delta2_2...
          +y2(1)*(1-y2(1))*w3(6,2)*delta2_3;
      delta3_2 = y2(2)*(1-y2(2))*w3(4,3)*delta1...
          +y2(2)*(1-y2(2))*w3(5,3)*delta2_2...
          +y2(2)*(1-y2(2))*w3(6,3)*delta2_3;
      delta3_3 = y2(3)*(1-y2(3))*w3(4,4)*delta1...
          +y2(3)*(1-y2(3))*w3(5,4)*delta2_2...
          +y2(3)*(1-y2(3))*w3(6,4)*delta2_3;
      for k = 1:4
         if k == 1 
            w3(1,k) = w3(1,k) + eta*bias(1,1)*delta3_1;
            w3(2,k) = w3(2,k) + eta*bias(1,2)*delta3_2;
            w3(3,k) = w3(3,k) + eta*bias(1,3)*delta3_3;
            w3(4,k) = w3(4,k) + eta*bias(1,4)*delta2_1;
            w3(5,k) = w3(5,k) + eta*bias(1,5)*delta2_2;
            w3(6,k) = w3(6,k) + eta*bias(1,6)*delta2_3;
            w3(6,k) = w3(7,k) + eta*bias(1,7)*delta1;
         else 
            w3(1,k) = w3(1,k) + eta*input_quad(j,1)*delta3_1;
            w3(2,k) = w3(2,k) + eta*input_quad(j,2)*delta3_2;
            w3(3,k) = w3(3,k) + eta*input_quad(j,3)*delta3_3;
            w3(4,k) = w3(4,k) + eta*y2(k-1)*delta2_1;
            w3(5,k) = w3(5,k) + eta*y2(k-1)*delta2_2;
            w3(6,k) = w3(6,k) + eta*y2(k-1)*delta2_3;
            w3(7,k) = w3(7,k) + eta*y3(k-1)*delta1;
         end
      end
   end   
end
plot(x0,y0,'k',x0,out,'r--')
title('Two hidden layers for quadratic regression')
legend('Expected Results','Neural Network Training Results')
function y=s(x)
y = 1./(1+exp(-x));
end
function ppt(weight,input,output,p,n)
if n == 1
b = output'-weight*input'; 
b = sum(b)/2;
x = 0:0.1:2;
y = -weight(1)/weight(2)*x-b/weight(2);
plot(x,y,'k')
hold on
for j = 1:4
    a = hardlim(weight*input(j,:)'+b);
    if a ==1
        plot(input(j,1),input(j,2),'rx')
    else
        plot(input(j,1),input(j,2),'b+')
    end
end
for i = 1: 10
    a = hardlim(weight*p(i,:)'+b);
    if a == 1
        plot(p(i,1),p(i,2),'rx');
    else
        plot(p(i,1),p(i,2),'b+');
    end
    pause(0.1)
end
else
b = output'-weight*input'; 
b1 = sum(b)/4;
b2 = sum(b)-b1;
x = 0:0.1:2;
y1 = -weight(1)/weight(2)*x-b1/weight(2);
y2 = -weight(1)/weight(2)*x-b2/weight(2);

plot(x,y1,'k',x,y2,'k')
hold on
for j = 1:4
    a1 = hardlim(weight*input(j,:)'+b1);
    a2 = hardlim(weight*input(j,:)'+b2);

    if a1 == 1 && a2 == 0
        plot(input(j,1),input(j,2),'rx')
    else if a1 == 0 ||a2 == 1
        plot(input(j,1),input(j,2),'b+')
        end
    end
end
for i = 1: 10
    a1 = hardlim(weight*p(i,:)'+b1);
    a2 = hardlim(weight*p(i,:)'+b2);

    if a1 ==1&&a2==0
        plot(p(i,1),p(i,2),'rx');
    else if a1 ==0||a2==1
        plot(p(i,1),p(i,2),'b+');
        end
    end
    pause(0.1)
end
end
hold off
box off
xlim([0 2])
ylim([0 2])
xticks([0 1 2])
yticks([0 1 2])
end
