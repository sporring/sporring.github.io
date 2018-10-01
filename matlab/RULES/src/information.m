function E = information(P,q);
%INFORMATION The information theoretic generalized entropy
%
%       E = information(P,q)
%         E - the entropy of order q (in bits, i.e. log2 is used)
%         P - a 1D or 2D distribution
%         q - the order, if q is a matrix then a matrix is returned
%             where each point is the corresponding generalized
%             entropy.
%
%       This function calculates the information of order q of a 1D or
%       2D distribution P, i.e. for every point 0 <= P(x) <= 1,
%       and sum(sum(P)) = 1. The log2 unit is used.  This measure is
%       not defined for q = 1, but it can be shown, that the limit for
%       q going to 1 is Shannon's Entropy.  Thus a first order
%       approximation is used close to 1.  At infinum the logarithm of
%       the max- or min-norm is returned
%
%       The complexity of the algorithm is O(n*m) where n and m are the
%       dimensions of P.
%
%       To calculate the entropy of a uniform distribution
%       (as specified) use,
%         E = information(ones(10,10)/100,0);
%
%       Copyright: Jon Sporring, June 13, 1996

E = zeros(size(q));
for i = 1:size(q,1)
  for j = 1:size(q,2)
    if q(i,j) == Inf
      % Essentially the logarithm of the max-norm
      E(i,j) = -log2(max(max(P)));
    else 
      if q(i,j) == -Inf
	% Essentially the logarithm of the min-norm
	E(i,j) = -log2(min(min(P+(P<=0))));
      else
        Pq = P.^q(i,j);
        X = sum(sum(Pq));
        if abs(q(i,j)-1) < 0.5
          % E(P,q) -> entropy(P) for q -> 1
          E(i,j) = -sum(sum(Pq.*log2(P+(P<=0))))/X;
        else 
	  E(i,j) = log2(X)/(1-q(i,j));
	end
      end
    end
  end
end

