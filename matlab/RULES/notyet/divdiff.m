function [D, t] = divdiff(x,f,c)
%DIVDIFF Calculates the divided differences vector for a polynomie in Newton form.
%       
%       [D,t] = divdiff(x,f,c)
%         D - A vector of coefficients obtain the divided differences table.
%         t - A vector of indexes into x that corresponds to the coefficients.
%         x - A vector of x-values
%         f - A table of function values and it's derivatives:
%              [f(x(i)), f'(x(i)), ... ]
%         c - A vector of the number of derivatives to take into account from
%             F at each x.
%
%       This function calculates the coefficient of a polynomie in
%       Newton form which has order c(i)+1 contact in each point with the
%       with a given function at x(i).
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, January 15, 1997


t = zeros(sum(c),1);
k=1;
for i = 1:max(size(c))
  for j = 1:c(i)
    t(k) = i;
    k = k + 1;
  end
end

D = f(t,1);
%disp([999.9999 x(t)'])
%disp([999.9999 D'])
for k = 1:size(D,1)-1
  for i = size(D,1)-k:-1:1
    dx = x(t(i+k))-x(t(i));
    if dx ~= 0
      D(i+k) = (D(i+k)-D(i+k-1))/dx;
    else
      D(i+k) = f(t(i),k+1);
    end
  end
%  disp([dx, D'])
end

t = x(t);
