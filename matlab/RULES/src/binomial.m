function y = binomial(x,n,p)

y=zeros(size(x));
for i = 1:size(x,1)
  for j = 1:size(x,2)
    y(i,j) = prod(n:-1:n-x(i,j)+1)/prod(1:x(i,j))*p^x(i,j)*(1-p)^(n-x(i,j));
  end
end
