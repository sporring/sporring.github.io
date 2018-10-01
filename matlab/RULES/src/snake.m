function x = snake(x,alpha,beta,f,gamma,w)
% SNAKE an implementation of Kass's snake
%
%   x = snake(x,alpha,beta,f,gamma,w)
%
%   x     - the new and old curve as an nx2 matrix
%   alpha - the membrane weight
%   beta  - the thin-plate weight
%   f     - the image-energy image, x must be inside the borders of f
%   gamma - the inverse step-size for implicit solution
%   w     - the number of iterations to perform
%
% Snake is an implementation of Kass, Witkin, and Terzopoulos, "Snakes:
% Active Contour Models", International Journal of Computer Vision,
% pp. 321-331, 1988, with my corrections.
%
%                         Copyright, Jon Sporring, DIKU, December 2003.

  [fc,fr] = gradient(f);
 
  % Calculate the circular derivative matrix
  t = [beta,-alpha-4*beta,6*beta+2*alpha,-alpha-4*beta,beta];

  % If length(t) is uneven, then the filter is thought centered at
  % position 1, if length(t) enen, then the filter is thought centered at
  % 0.5.
  n = floor(length(t)/2);
  A = toeplitz([t(n+1:-1:1),zeros(1,size(x,1)-length(t)),t(end:-1:end-n+1)],...
               [t(n+1:end),zeros(1,size(x,1)-length(t)),t(1:n)]);
  
  % In the implicit, the system matrix is
  B = inv(gamma*eye(size(A))+A);
  
  for i = 1:w
    ind = sub2ind(size(f),round(x(:,1)),round(x(:,2)));
    
    %% Semi-Implicit Euler:
    x = [B*(gamma*x(:,1)+fr(ind)),B*(gamma*x(:,2)+fc(ind))];
  end
  
