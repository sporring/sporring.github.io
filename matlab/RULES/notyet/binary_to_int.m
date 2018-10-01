function y = binary_to_int(x)
%BINARY_TO_INT Convert a binary number into an integer
%
%       y = binary_to_int(x)
%         y - is an integer
%         x - is a binary number represented as a string of 0's and 1's 
%
%       binary_to_int returns the positive integer vector of
%       the binary matrix x.  x must be a horizontal (column) vector.
%
%                                          Jon Sporring, January 1, 1994

i = zeros(size(x,1),1);
for j = 0:(size(x,2)-1)
  i = i + x(:,size(x,2)-j)*2^j;
end

y = i;
