function J = fbm_image(p1,p2,p3)
%FBM_IMAGE calculates a random fractal Brownian Motion Image
%
%       Fbm_image as two calling possibilities:
%
%       J = fbm_image(rows,cols,alpha);
%         J - The Fractal Brownian Motion image
%         rows, cols - the dimensions of the image
%         alpha - the coefficient for the 1/f^alpha power law
%
%       or
%
%       J = fbm_image(I,alpha);
%         J - The Fractal Brownian Motion image
%         I - an array of size [rows/2+1,cols]
%         alpha - the coefficient for the 1/f^alpha power law
%
%       Fbm_image calculates a real random image, which has a power spectrum
%       according to  1/f^alpha and.  An example is given below:
%
%         J = fbm_image(128,128,2);
%         colormap(gray(256));
%         imagesc(J);
%
%       The phase matrix may also be specified by the user:
%
%         rows = 128; cols = 128;
%         T = rand(rows/2+1,cols)*2*pi;
%         for i = linspace(1,3,10);
%           J = fbm_image(T,i);
%           contour(J,[0,0]);    
%           drawnow;
%         end
%
%       Copyright: Jon Sporring, January 13, 1997

if nargin == 3
  rows = p1;
  cols = p2;
  alpha = p3;
  T = rand(rows/2+1,cols)*2*pi;
else
  if nargin == 2
    T = p1;
    rows = 2*(size(p1,1)-1);
    cols = size(p1,2);
    alpha = p2;
  else
    error('Wrong number of parameters');
  end
end

J = zeros(rows,cols);
A = zeros(size(J));
A(1:rows/2+1,:) = T;

J(1:rows/2+1, 1:cols/2+1) = (([rows/2:-1:0]'*ones(1,cols/2+1)).^2 + (ones(rows/2+1,1)*[cols/2:-1:0]).^2).^(alpha/2);
J(rows/2+1,rows/2+1) = 1;
J(1:rows/2+1, 1:cols/2+1) = 1./J(1:rows/2+1,1:cols/2+1);
J(rows/2+1,rows/2+1) = 0;
J(1:rows/2+1, cols/2+2:cols) = fliplr(J(1:rows/2+1, 2:cols/2));
J(rows/2+2:rows, 1:cols) = flipud(J(2:rows/2, 1:cols));
A(rows/2+2:rows, :) = -flipud([A(2:rows/2,1),fliplr(A(2:rows/2, 2:cols))]);
A(rows/2+1, cols/2+2:cols) = -fliplr(A(rows/2+1, 2:cols/2));
J = J.*exp(i*A);
J = real(ifft2(fftshift(J)));
