function Is = scale2(FI,s,dr,dc);
%SCALE2  Gaussian Scale-Space using the Fourier Domain
%
%       FL = scale2(I,s,dr,dc)
%         I - the Fourier transform of a matrix or an image
%         s - the standard deviation of the Gaussian (must be larger than 0)
%             (default 1)
%         dr - the derivative order in the direction of the rows
%              (default 0)
%         dc - the derivative order in the direction of the columns
%              (default 0)
%
%       This is an implementation of the Gaussian scale space on matrices.
%       The convolution is implemented in the Fourier domain and for that
%       reason the number of rows and columns of the matrix must be powers
%       of 2.
%       Fractional valued dr and dc are possible, but be warned the result
%       will probably be complex.
%       The complexity of this algorithm is O(n) where n is the total number
%       of elements of I.
%
%       To calculate an image of scale (variance) 2^2 use,
%         L = real(ifft2(scale2(fft2(I),2,0,0)));
%       To derive an image once in the direction of rows at scale 1^2 do,
%         L = real(ifft2(scale2(fft2(I),1,1,0)));
%
%       Copyright: Jon Sporring, January 1, 1996

if nargin < 4
  dc = 0;
end
if nargin < 3
  dr = 0;
end
if nargin < 2
  s = 1;
end
if nargin < 1
  error('No input!  Use: Is = scale(FI,s,dr,dc);');
end

if any([s,dr,dc] < 0)
  error('s, dr, and dc must be larger than zero');
else 
  if s == 0
    if (dr == 0) & (dc == 0)
      Is = FI;
    else
      Is = zeros(size(FI));
    end
  else
    rows = size(FI,1);
    cols = size(FI,2);

    FG = zeros(rows,cols);
    DG = ones(rows,cols);
    if s == Inf
      FG(1,1) = 1;
      Is = FI.*FG;
    else
      % Calculate the Fourier transform of a gaussian fct. and the
      % Differentiation matrix.
      if (rows > 1) & (cols > 1) 
	if any(rem(size(FI),2) ~= 0)
	  error('Each image sides must be divisible by 2');
	else
	  % 2 dimensional image
	  x = [0:rows/2]'/rows;
	  y = [0:cols/2]/cols;
	  
	  FG(1:rows/2+1, 1:cols/2+1) = exp(-(repmat(((s*2*pi)*x).^2/2,[1,cols/2+1])+repmat(((s*2*pi)*y).^2/2,[rows/2+1,1])));  
	  FG(rows/2+1:rows, 1:cols/2+1) = flipud(FG(2:rows/2+1, 1:cols/2+1));
	  FG(:, cols/2+1:cols) = fliplr(FG(:, 2:cols/2+1));

	  if any([dr,dc] > 0)
	    DG(1:rows/2+1, 1:cols/2+1) = repmat(x.^dr,[1,cols/2+1]).*repmat(y.^dc,[rows/2+1,1])*(sqrt(-1)*2*pi)^(dr+dc);
	    DG(rows/2+1:rows, 1:cols/2+1) = (-1)^dr*flipud(DG(2:rows/2+1, 1:cols/2+1));
	    DG(:,cols/2+1:cols) = (-1)^dc*fliplr(DG(:, 2:cols/2+1));
	  end
	end
      else
	if rem(length(FI),2) ~= 0
	  error('Number of vector elements must be divisible by 2');
	else
	  % 1 dimensional image
	  [val,ind] = max([rows,cols]);
	  x = [0:val/2]'/val;
	  FG(1:val/2+1) = exp(-x.^2*(s*2*pi)^2/2);
	  FG(val/2+1:val) = flipdim(FG(2:val/2+1),ind);

	  DG = (j*2*pi*x').^dr;
	  DG(val/2+1:val) = flipdim(DG(2:val/2+1),ind);
	end
      end

      Is = FI.*FG.*DG;
    end
  end
end
