function J = dilation(I, t, step);
%DILATION Morphological Gray-value dilation
%
%       J = dilation(I, iterations, stepsize);
%         J - the resulting double matrix
%         I - the double matrix original matrix
%         t - an optional diffusion time
%             (0 <= t, default 1)
%         stepsize - an optional time step size 
%             (0 < stepsize <= sqrt(1/2), default sqrt(1/2))
%
%       Copyright: Jon Sporring, October 24, 2000

  if nargin < 3
    step = 1/sqrt(2);
  elseif step <= 0
    warning('stepsize must be larger than 0');
    step = 1/sqrt(2);
  elseif step > 1/sqrt(2)
    warning(sprintf('stepsize must be smaller than %.2f',1/sqrt(2)));
    step = 1/sqrt(2);
  end
  if nargin < 2
    t = 1;
  elseif t < 0
    warning('t must be larger than or equal 0');
    t = 1;
  end
  if nargin < 1
    error('Usage: J = dilation(I, t, stepsize);');
  else
    J = double(I);
    iter = ceil(t/step);
    for i = 1:iter
      Jt = sqrt(min(J-J([2,1:end-1],:),0).^2 ...
		+ max(J([2:end,end-1],:)-J,0).^2 ...
		+ min(J-J(:,[2,1:end-1]),0).^2 ...
		+ max(J(:,[2:end,end-1])-J,0).^2);
      J = J+step*Jt;
    end
  end
