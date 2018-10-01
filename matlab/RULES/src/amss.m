function J = amss(I, t, step);
%AMSS Affine Morphological Scale Space
%
%       J = amss(I, t, stepsize);
%         J - the resulting double matrix
%         I - the double matrix original matrix
%         t - an optional diffusion time
%             (0 <= t, default 1)
%         stepsize - an optional time step size 
%             (0 < stepsize <= 0.1, default 0.1)
%
%       Copyright: Joachim Weickert & Jon Sporring, October 24, 2000

  if nargin < 3
    step = 0.1;
  elseif step <= 0
    warning('stepsize must be larger than 0');
    step = 0.1;
  end
  if nargin < 2
    t = 1;
  elseif iter < 0
    warning('t must be larger than or equal 0');
    t = 1;
  end
  if nargin < 1
    error('Usage: J = amss(I, t, stepsize);');
  else
    J = double(I);
    iter = ceil(t/step);
    for i = 1:iter
      dx = (J([2:end,end-1],:)-J([2,1:end-1],:))/2;
      dy = (J(:,[2:end,end-1])-J(:,[2,1:end-1]))/2;
      dxx = J([2:end,end-1],:)-2*J+J([2,1:end-1],:);
      dyy = J(:,[2:end,end-1])-2*J+J(:,[2,1:end-1]);
      dxdy = dx.*dy;
      dxy = (dxdy < 0).*(J([2:end,end-1],[2:end,end-1])-J(:,[2:end,end-1])-J([2:end,end-1],:)+J ...
			 +J([2,1:end-1],[2,1:end-1])-J(:,[2,1:end-1])-J([2,1:end-1],:)+J)/2 ...
	    +(dxdy >= 0).*(-J([2,1:end-1],[2:end,end-1])+J(:,[2:end,end-1])+J([2:end,end-1],:)-J ...
			   -J([2:end,end-1],[2,1:end-1])+J(:,[2,1:end-1])+J([2,1:end-1],:)-J)/2;

      tmp = dy.^2.*dxx+dx.^2.*dyy-2*dx.*dy.*dxy;
      Jt = sign(tmp).*(abs(tmp).^(1/3));

      J = J+step*Jt;
    end
  end
