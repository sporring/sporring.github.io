function FIs = scalen(FI,s,d);
%SCALEN  Gaussian Scale-Space using the Fourier Domain
%
%       FL = scalen(I,s,d)
%         I - the Fourier transform of a matrix or an image
%         s - a scalar or list of standard deviations of the Gaussian (default 1)
%         d - a scalar or list of derivative orders (default 0)
%
%       This is an implementation of the Gaussian scale space on matrices.
%       The convolution is implemented in the Fourier domain and for that
%       reason the number of rows and columns of the matrix must be powers
%       of 2.
%       If either s or d are scalars, then the value is used along all
%       dimensions. Fractional valued d's are possible, but be warned the
%       result will probably be complex.
%       The complexity of this algorithm is O(n) where n is the total number
%       of elements of I.
%
%       To calculate an 2D image of scale (variance) 2^2 use,
%         L = real(ifftn(scalen(fftn(I),[2,2],[0,0])));
%       To derive an 3D image once in the direction of rows at scale 1^2 do,
%         L = real(ifftn(scalen(fftn(I),[1,1,1],[1,0,0])));
%
%       Copyright: Jon Sporring, January 1, 1996

if nargin < 3
  d = zeros(1,ndims(FI));
end
if nargin < 2
  s = ones(1,ndims(FI));
end
if nargin < 1
  error('No input!  Use: FL = scalen(FI,s,d);');
end

if (any(s < 0) | any(d < 0))
  error('Negative scales or derivative powers are not allowed.');
else 
  if length(s) == 1
    s = s*ones(1,ndims(FI));
  end
  if length(d) == 1
    d = d*ones(1,ndims(FI));
  end
  if (length(s) ~= ndims(FI)) | (length(d) ~= ndims(FI))
    error('s and d must be scalars or vectors of length ndims(FI).')
  else
    if all(s == 0)
      if all(d == 0)
	FIs = FI;
      else
	FIs = zeros(size(FI));
      end
    else
      FG = zeros(size(FI));
      DG = ones(size(FG));
      if any(rem(size(FI),2) ~= 0)
	error('Each image sides must be divisible by 2');
      else
	% Calculate the Fourier transform of a gaussian fct., FG, and the
	% differentiation matrix, DG.
	for i = 1:ndims(FG)
	  tiling = size(FG);
	  tiling(i) = 1;
	  reshaping = ones(1,ndims(FG));
	  reshaping(i) = size(FG,i);
	  x = [0:size(FG,i)/2-1,-size(FG,i)/2:-1]/size(FG,i);
	  FG = FG+repmat(reshape((2*pi*s(i)*x).^2,reshaping),tiling);
	  DG = DG.*repmat(reshape(x.^d(i),reshaping),tiling);
	end
	FG = exp(-FG/2);
	DG = DG.*(sqrt(-1)*2*pi)^sum(d);
	
	FIs = FI.*FG.*DG;
      end
    end
  end
end

