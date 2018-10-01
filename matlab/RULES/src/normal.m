function y = normal(x,m,d)

if d == 0
  y = inf*ones(size(x));
elseif d == inf
  y = zeros(size(x));
else
  y = exp(-((x-m)/d).^2/2)/(d*sqrt(2*pi));
end
