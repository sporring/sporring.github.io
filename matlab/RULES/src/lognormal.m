function y = lognormal(x,m,d)

y = zeros(size(x));
index = find(x > 0);
y(index) = exp(-((log(x(index))-m)/d).^2/2)./(d*x(index)*sqrt(2*pi));
