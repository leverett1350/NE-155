function [X,i] = jacobi(A,B,P,tol,maxi)

% Input  - A is an N x N nonsingular matrix
%        - B is an N x 1 matrix
%        - P is an N x 1 matrix; the initial guess
%	      - tol is the tolerance for P
%	      - maxi is the maximum number of iterations allowed
% Output - X is an N x 1 matrix: the jacobi approximation to
%	        the solution of AX = B
%         - i is the number of iterations


N = length(B);
iter = 0;
for k = 1:maxi
   for j = 1:N
      X(j) = (B(j)-A(j,[1:j-1,j+1:N])*P([1:j-1,j+1:N]))/A(j,j);
      iter = iter+1;
   end
   err = abs(norm(X'-P));
   relerr = err/(norm(X)+eps);
   P = X';
      if (err<tol)|(relerr<tol)
     break
   end
end
i = iter;
X = X';
