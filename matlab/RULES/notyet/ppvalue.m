function value = ppvalue(x, j, B, C)
%PPVALUE Calculate the value of the j'th derivative of a piece-wise polynomial
%       
%       value = ppvalue(x, j, B, C)
%         x - a matrix of points to be evaluated
%         j - the derivative (0 for the actual value)
%         B - the list of l+1 breakpoints in the function
%         C - the kxl matrix of right derivatives at the breakpoints
%             i.e. C(n,m) = D^{n-1}f(B(m))
%
%       Ppvalue implements the ppvalu function from de Boor: "A practical
%       Guide to Splines", 1978, Springer Verlag, pp. 89, and evaluates the
%       j'th derivative of a piece-wise polynomial at x.
%
%       Copyright:IBM Almaden Research Center & Jon Sporring, January 15, 1997

value = zeros(size(x,1), size(x,2));
i = interval(x, B);
exponent = [0:(size(C,1)-1-j)];

for n = 1:size(x,1)
  for m = 1:size(x,2)
    X = ((x(n,m)-B(i(n,m))).^exponent)./gamma(exponent+1);

    value(n,m) = X*C(j+1:size(C,1),i(n,m));
  end
end

