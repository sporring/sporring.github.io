function J = power2expand(I);
%POWER2EXPAND Expand image to nearest higher power of 2 using imresize
%
%       J = power2expand(I)
%         J - the resized matrix or vector
%         I - any matrix or vector
%
%       This does a imresize of an image to the nearest higher power of 2
%       for e.g. FFT Fourier Transform algorithms.  This function requires
%       the image toolbox!
%
%       Copyright: Jon Sporring, September 3, 1996

m = 2^floor(log2(size(I,1)));
n = 2^floor(log2(size(I,2)));
if m < size(I,1) | n < size(I,2)
  if m < size(I,1)
    m = 2*m;
  end
  if n < size(I,2)
    n = 2*n;
  end

  J = imresize(I,[m,n]);
else
  J = I;
end
