function I = RandomIsingImage(alpha,north,west,south,east,T,w,N)
% RandomIsingImage - generate a random image drawn from the Ising model
%
%  I = RandomIsingImage(alpha,north,west,south,east,T,w,N)
%
%    I           - The resulting image
%    alpha,
%    north,
%    west,
%    south,
%    east,
%    T           - Parameters to the Ising model: 
%                  p(I(x,y) | neigbors)  
%                     = (1/Z) exp(-(1/T)*[I(x,y)*alpha
%                                 +north*I(x+1,y)+west*I(x,y-1)
%                                 +south*I(x-1,y)+east*I(x,y+1)]
%                  where Z is the normalizing constant.
%    w           - The number of sweeps
%    N           - The size of the resulting image (NxN)
%
% The Ising model is a Markov Random Field for binary images.  In this
% implementation, we do not alter pixels at the one-pixel border around
% the image. 
% 
%                                    Copyright, Jon Sporring, DIKU, 2003.
  
  I = rand(N,N)>0.5;
  % In the following, we assume that indeces are inside a 1-pixel border
  [X,Y] = meshgrid([2:N-1],[2:N-1]);
  ind = sub2ind([N,N],X,Y);
  ind = ind(:);
  for i = 1:w
    % To ensure that we visit every pixel once per sweep in a random
    % order, we use randperm on the indices.
    rand_ind = randperm(length(ind));
    for n = 1:length(rand_ind)
      [x,y] = ind2sub([N,N],ind(rand_ind(n)));
      
      % potential(0,...)=0.
      p0 = 1;
      p1 = exp(-potential(I,1,x,y,alpha,north,west,north,west)/T);
      Z = 1+p1;
      p0 = p0/Z;
      p1 = p1/Z;
      
      % Choose white (1) with probability p1.
      I(x,y) = rand < p1;
    end
  end

function U = potential(I,v,x,y,alpha,north,west,south,east)
% Potential - Calculate the Gibbs potential for an Ising model
  
  U = v*(alpha+north*I(x+1,y)+west*I(x,y-1)+south*I(x-1,y)+east*I(x,y+1));
  
