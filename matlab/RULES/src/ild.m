function J = ild(I,t,step,method)
% ILD: Isotropic Linear Diffusion using various implementations.
%
%       J = ild(I,t,stepsize,method);
%         J - the resulting double matrix
%         I - the double matrix original matrix
%         t - an optional diffusion time
%             (0 <= iterations, default 1)
%         stepsize - an optional step size
%             (0 < stepsize < 0.25, default 0.1)
%         method - an optional string 
%             ('explicit' or 'implicit', default 'explicit')
%  
%       The Isotropic Linear Diffusion is implemented using the either the
%       explicit or the implicit solution to the heat equation with
%       homogeneous Neumann boundary conditions.  The computational
%       complexity of the algorithm is linearly proportional to the number
%       of pixels and iterations.  The implicit method is approximately
%       10 times slower than the explicit, but the implicit method allows
%       for abitrary high step sizes.
%
%       Corrected loop in the explicit scheme.  Jon Sporring January 4,
%       2001.
%       Corrected default stepsize. Jon Sporring, November 22, 2006.
%  
%       Copyright: Jon Sporring, October 19, 2000
  if nargin < 4
    method = 'explicit';
  end
  if nargin < 3
    step = 0.1;
  elseif step <= 0
    warning('stepsize must be larger than 0');
    step = 0.1;
  elseif (step > 0.25) & (lower(method(1)) == 'e')
    warning('stepsize must be smaller than 0.25 for ''explicit'' method');
    step = 0.1;
  end
  if nargin < 2
    t = 1;
  elseif t < 0
    warning('t must be larger than or equal 0');
    t = 1;
  end
  if nargin < 1
    error('Usage: J = ild(I,t,stepsize,method);');
  else
    J = double(I);
    clear I;

    iter = ceil(t/step);

    switch lower(method(1))
     case 'e'
      switch ndims(J)
       case 1
	for i = 1:iter
	  Jt = J([2:end,end-1])-2*J+J([2,1:end-1]);
	  J = J+step*Jt; 
	end
       case 2
	for i = 1:iter
	  Jt = J([2:end,end-1],:)-2*J+J([2,1:end-1],:) ...
	       +J(:,[2:end,end-1])-2*J+J(:,[2,1:end-1]);
	  J = J+step*Jt; 
	end
       case 3
	for i = 1:iter
	  Jt = J([2:end,end-1],:,:)-2*J+J([2,1:end-1],:,:) ...
	       +J(:,[2:end,end-1],:)-2*J+J(:,[2,1:end-1],:) ...
	       +J(:,:,[2:end,end-1])-2*J+J(:,:,[2,1:end-1]);
	  J = J+step*Jt; 
	end
       otherwise
	error('ild explicit has only been implemented for 0 < ndims < 4');
      end
     case 'i'
      % Joachim suggest:
      %   [L,U] = lu(speye(size(A)) - step*A);
      %   invL = full(inv(L));
      %   invU = full(inv(U));
      % However, since I don't know how to save space by LU decomposition in
      % Matlab, we might as well just use the inverse of A, which Matlab
      % is fast to compute. By the way, full representation is much
      % faster than sparse. 
      
      switch ndims(J)
       case 1
	A = sparse([-2 2 zeros(1,size(J,1)-2);
		    toeplitz([1,zeros(1,size(J,1)-3)],[1 -2 1 zeros(1,size(J,1)-3)]);
		    zeros(1,size(J,1)-2) 2 -2]);
	invA = full(inv(speye(size(A)) - step*A));
	for j = 1:iter
	  J = invA*J;
	end
       case 2
	A = sparse([-2 2 zeros(1,size(J,1)-2);
		    toeplitz([1,zeros(1,size(J,1)-3)],[1 -2 1 zeros(1,size(J,1)-3)]);
		    zeros(1,size(J,1)-2) 2 -2]);
	invA1 = full(inv(speye(size(A)) - step*A));
	A = sparse([-2 2 zeros(1,size(J,2)-2);
		    toeplitz([1,zeros(1,size(J,2)-3)],[1 -2 1 zeros(1,size(J,2)-3)]);
		    zeros(1,size(J,2)-2) 2 -2]);
	invA2 = full(inv(speye(size(A)) - step*A));
	for j = 1:iter
	  J = invA1*J;
	  % J(:,i) = invU*(invL*J(:,i));
	  J = (invA2*J')';
	  % J(i,:) = (invU*(invL*J(i,:)'))';
	end
       case 3
	A = sparse([-2 2 zeros(1,size(J,1)-2);
		    toeplitz([1,zeros(1,size(J,1)-3)],[1 -2 1 zeros(1,size(J,1)-3)]);
		    zeros(1,size(J,1)-2) 2 -2]);
	invA1 = full(inv(speye(size(A)) - step*A));
	A = sparse([-2 2 zeros(1,size(J,2)-2);
		    toeplitz([1,zeros(1,size(J,2)-3)],[1 -2 1 zeros(1,size(J,2)-3)]);
		    zeros(1,size(J,2)-2) 2 -2]);
	invA2 = full(inv(speye(size(A)) - step*A));
	A = sparse([-2 2 zeros(1,size(J,3)-2);
		    toeplitz([1,zeros(1,size(J,3)-3)],[1 -2 1 zeros(1,size(J,3)-3)]);
		    zeros(1,size(J,3)-2) 2 -2]);
	invA3 = full(inv(speye(size(A)) - step*A));
	for j = 1:iter
	  for j = 1:size(J,3)
	    J(:,:,j) = invA1*J(:,:,j);
	    J(:,:,j) = (invA2*J(:,:,j)')';
	  end
	  for i = 1:size(J,1)
	    for j = 1:size(J,2)
	      J(i,j,:) = permute(invA3*permute(J(i,j,:),[3,2,1]),[3,2,1]);
	    end
	  end
	end
       otherwise
	error('ild implicit has only been implemented for 0 < ndims < 4');
      end
     otherwise
      warning('method must be either ''explicit'' or ''implicit''');
    end
  end
