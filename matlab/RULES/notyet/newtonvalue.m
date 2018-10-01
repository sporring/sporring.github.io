function f = newtonvalue(D,t,x)
%NEWTONVALUE Calculates the value of a polynomial on Newton form
%       
%       f = newtonvalue(D,t,x)
%         f - The real value of the polynomial
%         D - A vector of coefficients on Newton form
%         t - A corresponding vector of zero points of the polynomie
%         x - a matrix of points to be evaluated
%
%       This function calculates the value of a polynomie in
%       Newton form as
%          f = D(1) + (x-t(1))*[D(2) + (x-t(2))*[...]]
%       with a given function at x(i).
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, January 15, 1997

f = D(max(size(D)))*ones(size(x));
for i = max(size(D))-1:-1:1
  f = D(i) + (x-t(i)).*f;
end
