function y = sign2(x)
%SIGN2  Alternative sign calculation.
%
%       y = sign2(x)
%         y - the sign matrix or vector of x
%         x - the input matrix/vector
%
%       Sign2 returns the sign of X, where sign2(0) = 1.
%       This is an extension to the sign function, which
%       evaluates sign(0) = 0.
%
%                                          Jon Sporring, January 1, 1994

for i = 1:size(x,1)
  for j = 1:size(x,2)
    if x(i,j) >= 0
      x(i,j) = +1;
    else
      x(i,j) = -1;
    end
  end
end

y = x;
