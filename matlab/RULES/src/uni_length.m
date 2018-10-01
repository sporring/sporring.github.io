function L = uni_length(A)
%UNI_LENGTH The ideal codelenth of Rissanens Universal Prior of Integers.
%
%       L = uni_length(A)
%         L - the resultant matrix
%         A - a matrix of integers > 0
%
%       Uni_length calculates the log^* function on each element in A, where
%       log^*(n) = c+log2(n)+log2(log2(n))+... over all positive terms.
%       This is only defined for positive integers.
%
%       Copyright: Jon Sporring, January 1, 1995

L = log2(2.865)*ones(size(A));
t = log2(A);
while any(any(t > 0))
  t = t.*(t>0);
  L = L + t;
  t = log2(t+1*(t==0));
end
