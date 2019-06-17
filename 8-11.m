

clc, clear, close all

N = 1000;
x = zeros(1,N);
p = zeros(1,N);

K = 0.5;
x(1) = 0.5;
p(1) = 0;
T = 1;

Jn = 1;
for n = 1:N-1
    p(n+1) = p(n) + K * sin( x(n) );
    x(n+1) = x(n) + T * p(n+1);
    
    In = [1+K*cos(2*pi*x(n)) 1
          K*cos(2*pi*x(n))   1];
    Jn = In * Jn;
end

eig(Jn)

figure
plot(x,p,'o')
xlabel 'x'
ylabel 'p'
