function J = mcm(I, t, step);
%MCM Mean Curvature Motion
%
%       J = mcm(I, t, stepsize);
%         J - the resulting double matrix
%         I - the double matrix original matrix
%         t - an optional diffusion time
%             (0 <= t, default 1)
%         stepsize - an optional time step size 
%             (0 < stepsize <= 0.25, default 0.25)
%
%       Copyright: Joachim Weickert & Jon Sporring, October 24, 2000

  if nargin < 3
    step = 0.25;
  elseif step <= 0
    warning('stepsize must be larger than 0');
    step = 0.25;
  end
  if nargin < 2
    t = 1;
  elseif t < 0
    warning('t must be larger than or equal 0');
    t = 1;
  end
  if nargin < 1
    error('Usage: J = mcm(I, t, stepsize);');
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
      
      Jt = (dx.^2.*dyy+dy.^2.*dxx-2*dx.*dy.*dxy)./(eps+dx.^2+dy.^2); 
      J = J+step*Jt;
    end
  end
