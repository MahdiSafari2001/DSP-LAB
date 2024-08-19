%Experiment 3-1
%section a
clc,close all;
w = linspace(0,pi,500);
f0 = 500;
fs = 10000;
w0 = 2*pi*f0/fs;
R = 0.99;
G = (1-R)*sqrt((1-2*R*cos(2*w0)+R^2));
H2 = (G^2)./((1-2*R*cos(w-w0)+R^2).*(1-2*R*cos(w-w0)+R^2));
plot(w,H2)
title('|H(jw)|^2')
xlabel('w')
ylabel('amp')
%section b
a1 = -2*R*cos(w0);
a2 = R^2;
% Initialize impulse response array
N =300;
h = zeros(1, N+1);

% Set h[0] (for n = 0)
h(1) = G;

% Calculate h[n] for n = 1 to N
for n = 1:N
    h(n+1) = -a1*h(n);
    if n > 1
        h(n+1) = h(n+1) - a2*h(n-1);
    end
end

% Plot impulse response
n = 0:N;
figure;
plot(0:N, h);
xlabel('n');
ylabel('h[n]');
title('Impulse Response');
h1 = (G/sin(w0))*((R.^n).*sin(w0*n+w0));
figure;
plot(0:N, h1);
xlabel('n');
ylabel('h[n]');
title('Impulse Response using equation 5-3');
%section c
n1 = 0:300;
sn = cos(w0*n1);
vn = randn(1,301);%creating noise
xn = sn + vn;
yn = zeros(1,301);
% Set y[0] (for n = 0)
yn(1) = G*xn(1);
% Calculate y[n] for n = 1 to 300
for n = 1:300
    yn(n+1) = -a1*yn(n) + G*xn(n);
    if n > 1
        yn(n+1) = yn(n+1) - a2*yn(n-1);
    end
end
figure;
plot(n1,yn,n1,sn);
xlabel('n');
ylabel('y[n]');
legend('recovered','original')
title('Response to noisy input vs Original signal');
%section d
yv = zeros(1,301);
% Set yv[0] (for n = 0)
yv(1) = G*vn(1);
% Calculate yv[n] for n = 1 to 300
for n = 1:300
    yv(n+1) = -a1*yv(n) + G*vn(n);
    if n > 1
        yv(n+1) = yv(n+1) - a2*yv(n-1);
    end
end
figure;
plot(n1,yv);
xlabel('n');
ylabel('y[n]');
title('Response to noisy input');
figure;
plot(n1,vn);
xlabel('n');
ylabel('v[n]');
title('noise');
%%
%section e
std_yv = std(yv);
std_vn = std(vn);
res1 = (std_yv/std_vn)^2;
res2 = (1+R^2)/((1+R)*(1+2*R*cos(w0)+R^2));
X = [' Result of 7-3 is ',num2str(res1),' and the result of 8-3 is ',num2str(res2)];
disp(X)
%%
%Experiment 3-2
%section a
numerator1 = [0.969531 -1.923772 0.969531];
denum1 = [1 -1.923772 0.939063];
numerator2 = [0.996088 -1.976468 0.996088];
denum2 = [1 -1.976468 0.992177];
figure;
freqz(numerator1,denum1)
title('H1(z)')
figure;
freqz(numerator2,denum2)
title('H2(z)')
n = 1:200;
x_in = n>=0;
y1 = filter(numerator1,denum1,x_in);
figure;
plot(n,y1)
title('filtering step input for H1')
y2 = filter(numerator2,denum2,x_in);
figure;
plot(n,y2)
title('filtering step input for H2')
%section b
ts = 0.0025;
sys1 = tf(numerator1,denum1, ts);
S1 = stepinfo(sys1,'SettlingTimeThreshold',0.01);
st1 = S1.SettlingTime;
sys2 = tf(numerator2,denum2, ts);
S2 = stepinfo(sys2,'SettlingTimeThreshold',0.01);
st2 = S2.SettlingTime
%section c
f1 = 4;
f2 = 8;
f3 = 12;
tt = 0:ts:6;
xt = zeros(size(tt));
xt = (tt>=0 & tt<2) .* cos(2*pi*f1*tt) + (tt>=2 & tt<4) .* cos(2*pi*f2*tt) + (tt>=4 & tt<6) .* cos(2*pi*f3*tt);
figure;
plot(tt,xt)
title('x(tn)')
xlabel('tn')
ylabel('x(tn)')
%filtering
yt = filter(numerator1,denum1,xt);
figure;
plot(tt,yt)
title('y(tn)for H1')
xlabel('tn')
ylabel('y(tn)')
%section e
%repeating last section scenario for H2(z)
yt1 = filter(numerator2,denum2,xt);
figure;
plot(tt,yt1)
title('y(tn) for H2')
xlabel('tn')
ylabel('y(tn)')
%section f
f = 0:0.01:20;
z = exp(i*2*pi*f*ts);
Hf1 = (0.969531 - 1.923772*(z.^(-1)) + 0.969531*(z.^(-2)))./(1 - 1.923772*(z.^(-1)) + 0.939063*(z.^(-2)));
Hf2 = (0.996088 - 1.976468*(z.^(-1)) + 0.996088*(z.^(-2)))./(1 - 1.976468*(z.^(-1)) + 0.992177*(z.^(-2)));
figure;
plot(f,abs(Hf1))
title('H1(f)')
xlabel('f')
ylabel('H1')
figure;
plot(f,abs(Hf2))
title('H2(f)')
xlabel('f')
ylabel('H2')
%%
%Experiment 3-2 last section : Repeating previous sections for a new H1 and H2
%section a
numerator1 = [0.030469 0 -0.030469];
denum1 = [1 -1.923772 0.939063];
numerator2 = [0.003912 0 -0.003912];
denum2 = [1 -1.976468 0.992177];
figure;
freqz(numerator1,denum1)
title('H1(z)')
figure;
freqz(numerator2,denum2)
title('H2(z)')
n = 1:200;
x_in = n>=0;
y1 = filter(numerator1,denum1,x_in);
figure;
plot(n,y1)
title('filtering step input for H1')
y2 = filter(numerator2,denum2,x_in);
figure;
plot(n,y2)
title('filtering step input for H2')
%section b
ts = 0.0025;
sys1 = tf(numerator1,denum1, ts);
S1 = stepinfo(sys1,'SettlingTimeThreshold',0.01);
st1 = S1.SettlingTime;
sys2 = tf(numerator2,denum2, ts);
S2 = stepinfo(sys2,'SettlingTimeThreshold',0.01);
st2 = S2.SettlingTime
%section c
f1 = 4;
f2 = 8;
f3 = 12;
tt = 0:ts:6;
xt = zeros(size(tt));
xt = (tt>=0 & tt<2) .* cos(2*pi*f1*tt) + (tt>=2 & tt<4) .* cos(2*pi*f2*tt) + (tt>=4 & tt<6) .* cos(2*pi*f3*tt);
figure;
plot(tt,xt)
title('x(tn)')
xlabel('tn')
ylabel('x(tn)')
%filtering
yt = filter(numerator1,denum1,xt);
figure;
plot(tt,yt)
title('y(tn)for H1')
xlabel('tn')
ylabel('y(tn)')
%section e
%repeating last section scenario for H2(z)
yt1 = filter(numerator2,denum2,xt);
figure;
plot(tt,yt1)
title('y(tn) for H2')
xlabel('tn')
ylabel('y(tn)')
%section f
f = 0:0.01:20;
z = exp(i*2*pi*f*ts);
Hf1 = (0.969531 - 1.923772*(z.^(-1)) + 0.969531*(z.^(-2)))./(1 - 1.923772*(z.^(-1)) + 0.939063*(z.^(-2)));
Hf2 = (0.996088 - 1.976468*(z.^(-1)) + 0.996088*(z.^(-2)))./(1 - 1.976468*(z.^(-1)) + 0.992177*(z.^(-2)));
figure;
plot(f,abs(Hf1))
title('H1(f)')
xlabel('f')
ylabel('H1')
figure;
plot(f,abs(Hf2))
title('H2(f)')
xlabel('f')
ylabel('H2')
%%
%Experiment 3-3
%section a
