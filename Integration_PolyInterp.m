% Ilan Halioua - 100472908

% format long
clear
clc

% x1 = 9 | x2 = 0 | x3 = 8

% t = 20 + 3*9 + 2*0 + 8 = 55


% METHOD 1: (DONE)
%   Direct integration of f_t (f_55)

fun = @(x) sin(55*x)./sqrt(x); % f55

n = zeros(1,12); % preallocates n
I1 = zeros(1,12); % preallocates I1

for k = 1:1:12
    n(k) = 2^k;
    I1(k) = integral(fun,1/(n(k)+1),1);
end
figure(1);
plot(n,I1)
title('I1-n');
subtitle('Integral of f55 from 1/n+1 to 1');
xlabel('n');
ylabel('Integral');

% METHOD 2:
%   Integration of f_55 using Natural Cubic Splines at equidistributed nodes

I2 = zeros(1,12); % preallocates I2

for k = 1:1:12
    % n(k) = 2^k; Already have it from method 1
    
    % (*SPLINES):

    r = zeros(1,13);
    y = zeros(1,13);
    for j = 1:1:n(k)+1 % n+1 nodes
        r(j) = (j)/(n(k)+1);
        y(j) = sin(55*r(j))/sqrt(r(j));
    end
    
    % SEL H m = u:
    
    % u: CHECK
    for i = 1:n(k)-1
        u(i) = ((y(i+2)-y(i+1))-(y(i+1)-y(i)))/(1/(n(k)+1)); 
    end 
    
    % H:
    H = zeros(n(k)-1);
    for i = 1:1:n(k)-1
       for j = 1:1:n(k)-1
           if i==j
               H(i,j) = (2/(n(k)+1))/3;
           elseif (j == i+1) || (j == i-1)
               H(i,j) = (1/(n(k)+1))/6;
           end
       end
    end
    
    % Solve SEL for m
    m = tridiagonal_matrix(H,u); % m = H\u;
    m = [0;m;0]; % natural spline
    
    % Calculating coefficients of cubic spline
    a0 = zeros(n(k),1);
    a1 = zeros(n(k),1);
    a2 = zeros(n(k),1);
    a3 = zeros(n(k),1);
    
    for j = 1:n(k) % CHECK
        a0(j) = y(j);
        a1(j) = (y(j+1)-y(j))/(1/(n(k)+1)) - ((1/(n(k)+1))*(2*m(j)+m(j+1)))/6;
        a2(j) = m(j)/2;
        a3(j) = (m(j+1)-m(j))/(6*(1/(n(k)+1)));
    end
    % Integral of S(x):
    I2(k) = 0;
    for j = 1:n(k)
        for d = 0:3
            if d == 0
                I2(k) = I2(k) + a0(j)/(((n(k)+1)^(d+1))*(d+1));
            elseif d == 1
                I2(k) = I2(k) + a1(j)/(((n(k)+1)^(d+1))*(d+1));
            elseif d == 2
                I2(k) = I2(k) + a2(j)/(((n(k)+1)^(d+1))*(d+1));
            elseif d == 3
                I2(k) = I2(k) + a3(j)/(((n(k)+1)^(d+1))*(d+1));
            end
        end
    end
end
figure(2);
plot(n,I2)
title('I2 - f55 integral');
subtitle('Natural Cubic Splines interpolation');
xlabel('n');
ylabel('Integral');

% METHOD 3:
%   Integration of f_55 using Hermite Interpolation

n = zeros(1,5);
I3 = zeros(1,5); % preallocates I3

for k = 1:1:5
    n(k) = 2^k;

    % (*FIND HermPoly):    
    r = zeros(1,n(k)+1);
    Q = zeros(2*n(k)+2,2*n(k)+2);
    for j = 0:1:n(k)
        r(j+1) = (j+1)/(n(k)+1); % rj
        Q(2*j+1,1) = sin(55*r(j+1))/sqrt(r(j+1)); % f(rj)
        Q(2*j+2,2) = -(sin(55*r(j+1))-110*r(j+1)*cos(55*r(j+1)))/(2*(r(j+1)^(3/2))); % f'(rj)
    end
    
    z = zeros(1,2*n(k)+2);
    for i = 0:n(k) % 0,...,n
        z(2*i+1) = r(i+1);
        z(2*i+2) = r(i+1);
        Q(2*i+2,1) = Q(2*i+1,1); % = f(r(i+1))
        if i ~= 0
            Q(2*i+1,2) = (Q(2*i+1,1)-Q(2*i,1))/(z(2*i+1)-z(2*i));
        end
    end
    for i = 2:2*n(k)+1
        for j = 2:i
            Q(i+1,j+1) = (Q(i+1,j)-Q(i,j))/(z(i+1)-z(i-j+1));
        end
    end

    syms x % (MORE PRECISE MTHD:)
    HP = Q((2*n(k)+1)+1,(2*n(k)+1)+1)*(x-z(2*n(k)+1));
    for i = 2:2*n(k)+1
        j = (2*n(k)+1)-i+1;
        HP = (HP+Q(j+1,j+1))*(x-z(j));
    end
    
    hermitePoly = HP + Q(1,1);
    I3(k) = int(hermitePoly,1/(n(k)+1),1);
end
figure(3);
plot(n,I3)
title('I3 - f55 integral');
subtitle('Hermite Interpolation');
xlabel('n');
ylabel('Integral');

% ERRORS

n = zeros(1,12); 
K = zeros(1,12); 
I2_err = zeros(1,12);
I3_err = zeros(1,5);

for k = 1:1:12
    K(k) = k;
    n(k) = 2^k;
    I2_err(k) = abs(I2(k)-I1(k));
end

for k = 1:1:5
    I3_err(k) = abs(I3(k)-I1(k));
end

% I1 vs I2 plot (n):
figure(4);
hold on

L(1) = plot(n,I1);
L(2) = plot(n,I2);
hold off
title('I1 vs I2 - n');
xlabel('n');
ylabel('Integral');
legend(L, {'I1 (f55) ~ Correct One', 'I2 (NCS)'},'Location','east')

% I1 vs I2 plot (K):
figure(5);
hold on

L(1) = plot(K,I1);
L(2) = plot(K,I2);
hold off
title('I1 vs I2 - k');
xlabel('k');
ylabel('Integral');
legend(L, {'I1 (f55) ~ Correct One', 'I2 (NCS)'},'Location','east') 

% I1 vs I2 vs I3 plot (n):
figure(6);
hold on

L(1) = plot(n,I1);
L(2) = plot(n,I2);
n3=[2,4,8,16,32];
L(3) = plot(n3,I3);
hold off
title('I1 vs I2 vs I3 - n');
xlabel('n');
ylabel('Integral');
legend(L, {'I1 (f55) ~ Correct One', 'I2 (NCS)', 'I3 (Hermite)'}) 

% I1 vs I2 vs I3 plot (K):
figure(7);
hold on

L(1) = plot(K,I1);
L(2) = plot(K,I2);
K3=[1,2,3,4,5];
L(3) = plot(K3,I3);
hold off
title('I1 vs I2 vs I3 - k');
xlabel('k');
ylabel('Integral');
legend(L, {'I1 (f55) ~ Correct One', 'I2 (NCS)', 'I3 (Hermite)'}) 

% I2 & I3 Error plot:

figure(8);
plot(n,I2_err)
title('I2 Error - n');
subtitle('Hermite Interpolation');
xlabel('n');
ylabel('Error');

figure(9);
plot(n3,I3_err)
title('I3 Error - n');
subtitle('Hermite Interpolation');
xlabel('n');
ylabel('Error');

% CONCLUSIONS & OBSERVATIONS

fprintf("\nFROM MY OBSERVATIONS:\n");
 
fprintf("\nAs n tends to 2^12, the errors of the other 2 methods\n"); 
fprintf("with respect to the direct Matlab integration one, are expected to tend to 0.\n");

fprintf("\nIf we zoom (figure 4) at n very large (n=4096 for example) we can see how I2 is really close to I1.\n");
fprintf("More precisely, around 0.000001 lower than I1 (See Fig.8 for more accurate value)\n");

fprintf("\nAlthough n=32 is the maximum we can analyse for I3 due to Matlab's symsengine limitations,\n");
fprintf("if I3 up to n=2^12 could be done, we could probaly see a similar approach to I1.\n");

fprintf("\nPlots show how at low values of n, both I2 and I3 have high errors with respect to I1.\n");
fprintf("Specifically at interval n ∈ (2,32) is where more deviation from the 'correct' integral can be seen.\n");

fprintf("At n = 8, a maximum error can be seen at both I2 and I3 of ~5 and ~0.88 respectively.\n");
fprintf("\nAt n = 32 the integral of I2 and I3 get very close to the one of I1 as one can see at figure 7,");
fprintf("\nso maybe one could expect I3 to get significantly closer to I1 from that point on.\n");

fprintf("\nIn conclusion, just from what can be compared (i.e. n ∈ [2,32]), errors calculated in this interval\n");
fprintf("show how Natural Cubic Splines carry out a much better job than Hermite interpolating polynomial for f55.\n");

% FUNCTIONS:

function x = tridiagonal_matrix(A,d)
    % A\d optimized to tridiagonal matrices by MATLAB FileExchange 
    % (Tamas Kis)
    
    % determines n
    n = length(d);
    
    % preallocates x
    x = zeros(n,1);
    
    % forward elimination
    for i = 2:n
        w = A(i,i-1)/A(i-1,i-1);
        A(i,i) = A(i,i)-w*A(i-1,i);
        d(i) = d(i)-w*d(i-1);
    end
    
    % backward substitution
    x(n) = d(n)/A(n,n);
    for i = (n-1):(-1):1
        x(i) = (d(i)-A(i,i+1)*x(i+1))/A(i,i);
    end
    
end
