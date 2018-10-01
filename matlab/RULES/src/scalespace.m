function Is = scalespace(FI, s, dr, dc)
%SCALESPACE calculates smoothed and possibly derived versions of FI
%
%       Is = scalespace(FI, s, dr, dc)
%         Is - The scale space of the images in spatial coordinates
%         FI - the Fourier transformed image
%         s  - a vector of scales
%         dr - the derivative order in the direction of the rows
%         dc - the derivative order in the direction of the columns
%
%       Scalespace calculates the derivative at all indicated scales.  The
%       result is given as an ndims+1 array, where the scale parameter is
%       the last coordinate.  Below is an example showing the evolution of
%       the mid isophote:
%
%         y = rand(128,128);
%         Y = fft2(y);
%         tau = [-Inf,linspace(-2,8,63)];
%         ys = scalespace(Y,sqrt(2*exp(tau)),0,0);
%         m = mean(mean(y));
%         for i = 1:size(ys,3)
%           contour(ys(:,:,i), [m,m]);
%           drawnow;
%         end
%
%       Copyright: Jon Sporring, January 13, 1997


Is = zeros([size(FI),max(size(s))]);

for i = 1:max(size(s))
    Is(:,:,i) = real(ifft2(scale(FI,s(i),dr,dc)));
end

if size(FI,1) == 1
   Is = squeeze(Is)';
else
  if size(FI,2) == 1;
    Is = squeeze(Is);
  end
end
