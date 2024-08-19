
function [] = singen(w,n)
%1-5
%creating unit function
n1 = linspace(-2,10,n);
u = zeros(size(n1));
u(n1>=0) = 1;
y= u.*sin(w*n1);
figure;
plot(n1,u)
title('Generating Sinusoidal Signal')
xlabel('time')
ylabel('amp')
end