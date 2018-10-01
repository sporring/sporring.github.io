function E = entropy(P);
%ENTROPY The information theoretic entropy
%
%       E = entropy(P)
%         E - the mean entropy in bits (log2 used)
%         P - a 1D or 2D distribution
%
%       This function calculates the entropy of a 1D or 2D distribution P,
%       i.e. every point of P must greater than or equal to 0 and the P
%       should be normalized (sum(sum(P)) = 1).  The log2 unit is used.
%       The complexity of the algorithm is O(n*m) where n and m are the
%       dimensions of P.
%
%       To calculate the entropy of a uniform distribution
%       (as specified) use,
%         E = entropy(ones(10,10)/100);
%
%       Copyright: Jon Sporring, January 1, 1996
%
  
P = P+(P<=0);
E = -sum(sum(P.*log2(P)));
