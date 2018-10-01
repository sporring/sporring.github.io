function y = beta(x,p,q,a,b)

y = zeros(size(x));
y = (x-a).^(p-1).*(b-x).^(q-1)/(b-a)^(p+q-1).*gamma(p+q)/(gamma(p)*gamma(q));
