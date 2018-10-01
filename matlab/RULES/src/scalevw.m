function res = scalevw(I,s,dv,dw)

% SCALEVW Gaussian Scale-Space using the function scale2
%
%         Is = scalevw(I,s,dr,dc)
%           I - the Fourier transform of a matrix or an image
%           s - the standard deviation of the Gaussian (must be larger than 0)
%           dv - the derivative order in the perpendicular direction to the gradient
%	    dw - the derivative order in the direction of the gradient
%
%         To derive an image twice in the direction of gradient at scale 2^2 do,
%         Lww = scalevw(I,2,0,2);
%         
%         Copyright: Ole Fogh Olsen, November 1999
%
%         See also scale2



fI = fft2(I);
[m n] = size(I);
order = dv+dw;

if order == 0, 
  res = real(ifft2(scale2(fI,s,0,0)));
else
  Ix = real(ifft2(scale2(fI,s,1,0)));
  Iy = real(ifft2(scale2(fI,s,0,1)));
  gradientSquared = Ix .^2 + Iy .^2;

  gradientMagOrder = gradientSquared .^ (order/2);

  derivativeTensor = zeros(m,n,order);
  for i=0:order,
    derivativeTensor(:,:,i+1) = real(ifft2(scale2(fI,s,order-i,i)));
  end

  for j=1:dw,
    for i=1:order+1-j,
      derivativeTensor(:,:,i) = (derivativeTensor(:,:,i) .* Ix + derivativeTensor(:,:,i+1) .* Iy);
    end
  end

  for j=1:dv,
    for i=1:order+1-dw-j,
      derivativeTensor(:,:,i) = derivativeTensor(:,:,i) .* Iy - derivativeTensor(:,:,i+1) .* Ix;
    end
  end

  res = derivativeTensor(:,:,1) ./ gradientMagOrder;
end


    

