function y = binomial(x,m)

y=m.^x*exp(-m);
for i = 1:size(x,1)
  for j = 1:size(x,2)
    y(i,j) = y(i,j)/prod(1:x(i,j));
  end
end
