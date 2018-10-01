function J = ced(I,t,step,lambda,sigma,rho,alpha);
%CED Coherence-Enhancing Anisotropic Diffusion
%
%       J = ced(I,t,stepsize,lambda,sigma,rho,alpha);
%         J - the resulting double matrix
%         I - the double matrix original matrix
%         t - an optional diffusion time
%             (0 <= iterations, default 1)
%         stepsize - an optional step size
%             (0 < stepsize < 0.25, default 0.25)
%         lambda - the optional contrast parameter
%             (0 <= lambda, default 1)
%         sigma - the optional noise parameter
%             (0 <= sigma, default 0.5)
%         rho - the optional integration parameter
%             (0 <= rho, default 4)
%         alpha - the optional linear difference fration parameter
%             (0 <= alpha, default 0.0001)
%
%       Copyright: Joachim Weickert & Jon Sporring, October 20, 2000

  if nargin < 7
    alpha = 0.0001;
  elseif alpha < 0
    warning('alpha must be larger than or equal 0');
    alpha = 0.0001;
  end
  if nargin < 6
    rho = 4;
  elseif rho < 0
    warning('rho must be larger than or equal 0');
    rho = 4;
  end
  if nargin < 5
    sigma = 4;
  elseif sigma < 0
    warning('sigma must be larger than or equal 0');
    sigma = 4;
  end
  if nargin < 4
    lambda = 1;
  elseif lambda < 0
    warning('lambda must be larger than or equal 0');
    lambda = 1;
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
    warning('t must be larger than or equal or equal 0');
    t = 1;
  end
  if nargin < 1
    error('Usage: J = ced(I,t,stepsize,lambda,sigma,rho,alpha);');
  else
    J = double(I);
    precision = 2;
    iter = ceil(t/step);
    % Convert from standard deviation to time:
    sigma = sigma^2/2;
    rho = rho^2/2;
    for i = 1:iter
      T = ild(J, sigma,sigma/ceil(sigma),'implicit');
      
      % Calculate the Structure Tensor
      Vx = (T([2:end,end-1],:)-T([2,1:end-1],:))/2;
      Vy = (T(:,[2:end,end-1])-T(:,[2,1:end-1]))/2;
      dxx = Vx.^2;
      dxy = Vx.*Vy;
      dyy = Vy.^2;
      dxx = ild(dxx, rho,rho/ceil(rho),'implicit');
      dxy = ild(dxy, rho,rho/ceil(rho),'implicit');
      dyy = ild(dyy, rho,rho/ceil(rho),'implicit');

      % Principal axis transformation. 
      % mu1 and mu2 are eigenvalues and v1 = (c,s)
      tmp = sqrt((dxx -dyy).^2 +4*dxy.^2);
      mu1 = (dxx +dyy +tmp)/2; 
      mu2 = (dxx +dyy -tmp)/2; 
      c = 2*dxy;
      s = dyy -dxx +tmp;
      norm = sqrt(c.^2+s.^2);
      c = c./(eps+norm);
      s = s./(eps+norm);

      % Set diffusion strength along eigenvectors
      coherence = (mu1-mu2).^2;
      lam1 = alpha;
      lam2 = alpha +(coherence>eps).*(1-alpha).*exp(-lambda./coherence);
      
      % Principal axis backtransformation of a symmetric (2*2)-matrix. 
      % A = U * diag(lam1, lam2) * U_transpose with U = (v1 | v2)     
      % v1 = (c, s) is first eigenvector
      dxx = c.^2.*lam1+s.^2.*lam2;
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
