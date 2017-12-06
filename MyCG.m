function [x,steps] = MyCG(A,b,x0,epsilon,iterMax)
% general Conjugate Gradient method to solve linear equations
%
%   % A: square matrix generated from B*B
%   % x: unknown vector
%   % espilon: gradient error, 1e-6 for default
%   % iterMax: max iteration times, 1000 for default
%
%Copyright 2017@Hanyu

r0 = b - A*x0;
p0 = r0;
if nargin < 3
    error('Error: no sufficient input!');
end
if nargin < 4
    epsilon = 1.0e-6;
end
if nargin < 5
    iterMax = 1000; % default iterations is1000
end
for steps=1:iterMax
    if abs(p0) < epsilon
        break;
    end
    a0 = r0'*r0/(p0'*A*p0); %alpha: search step length
    x1 = x0 + a0*p0; % iterative solution x
    r1 = r0 -a0*A*p0; %r: error, =x* - xi
    b0 = r1'*r1/(r0'*r0); %b: regularization coefficient 
    p1 = r1 + b0*p0; % p1&p0: conjugate gradient(search oriention)
    x0 = x1;
    r0 = r1;
    p0 = p1;
end
x = x0;
end