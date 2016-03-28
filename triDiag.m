function [T,B,P,P2,con] = triDiag(a,b,c,n)

% Input: a,b,c - Coefficient Values
%        n - n x n Tridiagonal Matrix
% Output: T - Tridiagonal Matrix
%         B - b Vector
%         P - Explicit Direct Solver
%         P2 - Backslash Operator Solver
%         con - Condition Number
%         also plots both solved values

T = b*diag(ones(n,1)) + c*diag(ones(n-1,1),1) + a*diag(ones(n-1,1),-1);
B = linspace(0,n-1,n);
con = cond(T);
B = B';
T2 = T;
T = inv(T);
P = T*B;
P2 = T2\B;
figure
subplot(2,1,1)
plot(B,P)
title('Explicity Solved')
xlabel('b vector values')
ylabel('x values')
subplot(2,1,2)
plot(B,P2,'g')
title('Backslash Operator Solved')
xlabel('b vector values')
ylabel('x values')
end