function J = ind(I,t,step,lambda,sigma);
%IND Isotropic Nonlinear Diffusion
%
%       J = ind(I,t,stepsize,lambda,sigma);
%         J - the resulting double matrix
%         I - the double matrix original matrix
%         t - an optional difusion time
%             (0 <= iterations, default 1)
%         stepsize - an optional step size
%             (0 < stepsize < 0.25, default 0.25)
%         lambda - the optional edge enhancement parameter
%             (0 <= lambda, default 7)
%         sigma - the optional regularization parameter
%             (0 <= sigma, default 1)
%
%       Copyright: Joachim Weickert & Jon Sporring, October 20, 2000.

  if nargin < 5
    sigma = 1;
  elseif sigma < 0
    warning('sigma must be larger than or equal 0');
    sigma = 1;
  end
  if nargin < 4
    lambda = 7;
  elseif lambda < 0
    warning('lambda must be larger than or equal 0');
    lambda = 7;
  end
  if nargin < 3
    step = 0.25;
  elseif step <= 0
    warning('stepsize must be larger than 0');
    step = 0.25;
  elseif step > 0.25
    warning('stepsize must be smaller than 0.25');
    step = 0.25;
  end
  if nargin < 2
    t = 1;
  elseif t < 0
    warning('t must be larger than or equal 0');
    t = 1;
  end
  if nargin < 1
    error('Usage: J = ind(I,t,stepsize,lambda,sigma);');
  else
    J = double(I);

    iter = ceil(t/step);
    % Convert from standard deviation to time:
    sigma = sigma^2/2;
    for i = 1:iter
      T = ild(J, sigma,sigma/ceil(sigma),'implicit');
      grad2 = ((T([2:end,end-1],:)-T([2,1:end-1],:)).^2 ...
	       +(T(:,[2:end,end-1])-T(:,[2,1:end-1])).^2)/4;
      dc = 1-exp(-3.31488./(eps+grad2.^4/lambda^8));
      Jt = ((dc([2:end,end-1],:)+dc).*(J([2:end,end-1],:)-J)...
	    +(dc([2,1:end-1],:)+dc).*(J([2,1:end-1],:)-J))...
	   +((dc(:,[2:end,end-1])+dc).*(J(:,[2:end,end-1])-J)...
	     +(dc(:,[2,1:end-1])+dc).*(J(:,[2,1:end-1])-J)); 
      J = J+step*Jt/2;
    end
  end
