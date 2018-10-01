function y = int_to_binary(x)
%INT_TO_BINARY Convert an integer to binary form.
%
%       y = int_to_binary(x)
%         y - is the binary number represented as string of 0's and 1's
%         x - is an integer
%
%       Int_to_binary returns the binary matrix of the positive
%       integer vector X.  The binary string is a horizontal vectors
%
%                                          Jon Sporring, January 1, 1994

b = zeros(size(x,2), ceil(log2(max(x))));

j = 1;
while any(x > 0)
  b(:,j) = rem(x,2)';
  x = fix(x./2);
  j = j + 1;
end

y = fliplr(b);
