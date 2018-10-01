function [K, J, E] = kmean(I,N,T);
%KMEAN  Calculate the matrix minimizing the squared distance using fewer values.
%
%       [K, J, E] = kmean(I,N,T)
%         K - is the set of (kluster-) values found to be minimized
%         J - the new matrix or vector using only values from K
%         E - the is the L-2 distance between J and I.
%         I - the original matrix or vector
%         N - the size of K (the maximum number of new values)
%         T - the number of different starting points that should
%             be investigated
%
%       Kmean is a simple algorithm for transforming an image (matrices or
%       vectors) into one with fewer intensity values.  The present
%       algorithm is the simplest but also slowest possible.  K is initiated
%       to be a random vector, and the squared distance is sought to be
%       minimized from the modified image to the original using only the
%       values from K. Thereafter K is modified as the center of the
%       belonging kluster.  Kluster-centers with no members are
%       reinitialized randomly.
%
%       To reduce the number of intensities to 4 (investigating only one
%       starting point) in matrix I use,
%         [K, J, E] = kmean(I, 4, 1);
%
%                                          Jon Sporring, January 1, 1994


if (N <= 0) | (T <= 0)
  error('N and T must be larger than zero');
else
  MAX = max(max(I));
  MIN = min(min(I));
  E_min = 256;
  E = zeros(T,1);

  for t = 1:T
    K = (MAX-MIN)*rand(N,1)+MIN;
    J = ones(size(I));

    changes = 1;
    while(changes == 1)
      changes = 0;
      % For all points calculate the nearest cluster
      E(t) = 0;
      for i = 1:size(I,1)
        for j = 1:size(I,2)
          [d, n] = min((K-I(i,j)).^2);
          J(i,j) = K(n);
          E(t) = E(t) + d;
        end
      end
      E(t) = sqrt(E(t)/(size(I,1)*size(I,2)));
%      disp(E(t))

      % For all points calculated the new clusters as the center of mass
      K_new = zeros(size(K));
      for n = 1:size(K,1)
        S = (J == K(n));
        if sum(sum(S)) == 0
          K_new(n) = (MAX-MIN)*rand+MIN;
        else
          K_new(n) = sum(sum(I.*S))/sum(sum(S));
        end
      end 
      if K_new ~= K
        K = K_new;
        changes = 1;
      end
    end
%    disp(E(t));
    if E(t) < E_min
      E_min = E(t);
      K_min = K;
    end
  end
%  disp(sprintf('min(E) = %6.2f, max(E) = %6.2f, mean(E) = %6.2f, std(E) = %6.2f', min(E), max(E), mean(E), std(E)));
  K = sort(K_min);
  E = E_min;
end
