function J = eed(I,t,step,lambda,sigma);
%EED Edge-Enhancing Anisotropic Diffusion
%
%       J = eed(I,t,stepsize,lambda,sigma);
%         J - the resulting double matrix
%         I - the double matrix original matrix
%         t - an optional diffusion time
%             (0 <= iterations, default 1)
%         stepsize - an optional step size
%             (0 < stepsize <= 0.25, default 0.25)
%         lambda - the optional edge enhancement parameter
%             (0 <= lambda, default 3)
%         sigma - the optional regularization parameter
%             (0 <= sigma, default 2.5)
%
%       Copyright: Joachim Weickert & Jon Sporring, October 20, 2000.

  if nargin < 5
    sigma = 2.5;
  elseif sigma < 0
    warning('sigma must be larger than or equal 0');
    sigma = 2.5;
  end
  if nargin < 4
    lambda = 3;
  elseif lambda < 0
    warning('lambda must be larger than or equal 0');
    lambda = 3;
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
    error('Usage: J = eed(I,t,stepsize,lambda,sigma);');
  else
    J = double(I);
    iter = ceil(t/step);
    % Convert from standard deviation to time:
    sigma = sigma^2/2;
    for i = 1:iter
      T = ild(J, sigma,sigma/ceil(sigma),'implicit');
      
      % Calculate the Diffusion Tensor
      Vx = (T([2:end,end-1],:)-T([2,1:end-1],:))/2;
      Vy = (T(:,[2:end,end-1])-T(:,[2,1:end-1]))/2;
      grad = sqrt(Vx.^2+Vy.^2);
      c = Vx./(eps+grad);
      s = Vy./(eps+grad);
      lam1 = 1-exp(-3.31488./(eps+(grad/lambda).^8));
      lam2 = 1;
      dxx = c.^2.*lam1+s.^2*lam2;
      dyy = lam1+lam2-dxx;
      dxy = c.*s.*(lam1-lam2);

      % Calculate the increment
      Jt = step/2 *(dxx([2:end,end-1],:)             +dxx) .*(J([2:end,end-1],:)             -J) ...
	   +step/2*(dxx([2,1:end-1],:)               +dxx) .*(J([2,1:end-1],:)               -J) ...
	   +step/2*(dyy(:,[2:end,end-1])             +dyy) .*(J(:,[2:end,end-1])             -J) ...
	   +step/2*(dyy(:,[2,1:end-1])               +dyy) .*(J(:,[2,1:end-1])               -J) ...
	   +step/4*(dxy([2:end,end-1],[2:end,end-1]) +dxy) .*(J([2:end,end-1],[2:end,end-1]) -J) ...
	   +step/4*(dxy([2,1:end-1],[2,1:end-1])     +dxy) .*(J([2,1:end-1],[2,1:end-1])     -J) ...
	   -step/4*(dxy([2,1:end-1],[2:end,end-1])   +dxy) .*(J([2,1:end-1],[2:end,end-1])   -J) ...
	   -step/4*(dxy([2:end,end-1],[2,1:end-1])   +dxy) .*(J([2:end,end-1],[2,1:end-1])   -J);

      % Perform step
      J = J+Jt;
    end
  end
