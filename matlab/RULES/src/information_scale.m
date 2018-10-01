function L = information_scale(I,s,q);
%INFORMATION_SCALE Calculate the information as function of order and scale
%
%       L = information_scale(I,s,q)
%         L - The information as function of order and scale
%         I - The image
%         s - The list of scales (2*variance)
%         q - The list of information order
%
%       This function evaluates the generalized entropy or information for
%       an image.  The rows corresponds to the scales and the columns to the
%       information order. Be warned: This is a slow function!
%
%       Copyright: Jon Sporring, January 1, 1996

L = zeros(size(s,2),size(q,2));
FI = fft2(extend(I,2,2));
for i = 1:size(s,2)
  I2 = real(ifft2(scale(FI,sqrt(s(i)/2),0,0)));
  I2 = I2(1:size(I,1),1:size(I,2));
  L(i,:) = information(I2/sum(sum(I2)),q);
end
