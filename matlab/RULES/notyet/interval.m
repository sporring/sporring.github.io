function i = interval(x, B)
%INTERVAL Find which interval x belongs to in B
%       
%       i = interval(x, B)
%         i - the matrix of indexes of the left interval boundary in B
%         x - a real valued matrix
%         B - a non-decreasing list of interval borders.
%
%       For each x(n,m) the index of the left border of the intervals in B
%       is returned.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, January 13, 1997

i = zeros(size(x,1),size(x,2));
for n = 1:size(x,1)
  % This inner loop can be implemented faster, but then again, if B is large
  % a more sophisticate search should be implemented.
  for m = 1:size(x,2)
    [y,i(n,m)] = min(abs(B(1:size(B,2)-1)-x(n,m)));
  end
  i(n,:) = i(n,:) - ((B(1,i(n,:))>x(n,:)) & (i(n,:)>1));
end
