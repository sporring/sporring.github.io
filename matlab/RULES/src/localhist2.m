function H = localhist2(I,sigma,beta,alpha,N)
%LOCALHIST2  Local soft histogram of an image.
%
%       H = localhist2(I,sigma,beta,alpha,N)
%         I - the original an image
%         sigma - the standard deviation of the Gaussian to smooth the
%             image with. (default 1)
%         beta - the standard deviation of the Gaussian to soften the
%             isophotes with. (default 1)
%         alpha - the standard deviation of the Gaussian window.
%              (default 1)
%         N - the maximum intensity value. (default max(max(I))).
%
%       This function returns the local histogram in all points of an image
%       for a specific setting of sigma, beta, and alpha.  Hence, the output
%       H, is of dimension ndims(I)+1, and the histogram will be sampled at
%       intensity values (0..N).  REMEMBER, the user most likely will want
%       to rescale intensity values to lie in the interval!
%
%       Example:
%         I = 16*rand(64,64);
%         H = localhist2(I,2,1,1,16);
%         plot(squeeze(H(32,32,:)));
%
%       Copyright: Jon Sporring, December 1, 1999

  % Checking parameters and setting default;
  if nargin < 1
    error('Image I must be supplied! Usage: H = localhist2(I,sigma,beta,alpha,N)');
  end
  if nargin < 2
    sigma = 1;
  end
  if nargin < 3
    beta = 1;
  end
  if nargin < 4
    alpha = 1;
  end
  if nargin < 5
    if isreal(I)
      N = ceil(max(max(I)));
    else
      warning(['You will get faster computation if you specify N, when Image is ' ... 
	       'given as its Fourier Transform.']);
      N = ceil(max(max(real(ifft2(I)))));
    end
  end
    
  % Smoothing image
  if sigma > 0
    if isreal(I)
      L = real(ifft2(scale2(fft2(I),sigma,0,0)));
    else
      L = real(ifft2(scale2(I,sigma,0,0)));
    end
  else
    if isreal(I)
      L = I;
    else
      L = real(ifft2(I));
    end
  end
  
  % Calculating soft isophote and sum under window
  H = zeros([size(I),N]);
  for i = 1:N
    % disp(i);
    if beta == 0
      J = (L==i);
    else
      J = 1/sqrt(2*pi*beta^2)*exp(-(L-i).^2/(2*beta^2));
    end
    
    if alpha > 0
      H(:,:,i) = real(ifft2(scale2(fft2(J),alpha,0,0)));
      % Clamp values to ensure positivity
      H(:,:,i) = (H(:,:,i)>0).*H(:,:,i);
      % figure(1); plot(H(1,:,i+1)); pause
    else
      H(:,:,i) = J;
    end
  end
  
  % Ensure that the histogram sums to 1
  H = H./repmat(sum(H,3),[1,1,size(H,3)]);
