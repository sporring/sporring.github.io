function c = extract(L,V,method)
%EXTRACT Calculates a function as a list of points at V in L.
%
%       c = extract(L,V,method)
%         c - the resulting measured function
%         L - the matrix to be used as basis
%         V - a list of points in L
%         method - the name of the interpolation method
%
%       Extract returns a list of function values measured at the points
%       specified in V from L using the interpolation method specified
%       with `method'.  `method' should be one of the names allowed by
%       the program Interpolate.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 5, 1997

c = zeros(1,size(V,2));
for i = 1:size(c,2)
  c(i) = interpolate(L,V(:,i),method);
end
