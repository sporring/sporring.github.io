function y = chi(x,n)

c = x.^(n/2-1).*exp(-x/2)./(2^(n/2).*gamma(n/2));
y = zeros(size(c));
index = find(x > 0);
y(index) = c(index);
