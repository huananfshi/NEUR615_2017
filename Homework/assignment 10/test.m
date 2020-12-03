eta = 0.7;
train_time = 10000;
w3 = -1 +2.*rand(7,4);
x0 = linspace(0,1,200);
input(:,1)=x0;
input(:,2:3)=zeros(200,2);
bias = [-1 -1 -1 -1 -1 -1 -1];
y0 = (x0-0.5).^2;
for i = 1:train_time
   out = zeros(200,1);
   N = length (input(:,1));
   for j = 1:N
      x1 = bias(1,1)*w3(1,1)+ input(j,1)*w3(1,2)+ input(j,2)*w3(1,3);
      x2 = bias(1,2)*w3(2,1)+ input(j,1)*w3(2,2)+ input(j,2)*w3(2,3);
      x3 = bias(1,3)*w3(3,1)+ input(j,1)*w3(3,2)+ input(j,2)*w3(3,3);
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
            w3(1,k) = w3(1,k) + eta*input(j,1)*delta3_1;
            w3(2,k) = w3(2,k) + eta*input(j,2)*delta3_2;
            w3(3,k) = w3(3,k) + eta*input(j,3)*delta3_3;
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
