function [G,H] = curvature(I,s)
%CURVATURE Finds the Gaussian and the Mean curvature of the image surface
%
%       [G,H] = eliptical(I,s)
%         G - an image of Gaussian curvatures
%         H - an image of Mean curvatures
%         I - the original image
%         s - the scale (std) of evaluation.
%
%       curvature evaluates the Hessian of an image in every point and
%       assigns the Gaussian curvature as the determinant and the Mean
%       curvature as the trace.
%
%       Copyright: Jon Sporring, November 7, 1997

FI = fft2(I);
Lxx = real(ifft2(scale(FI,2,0,s)));
Lxy = real(ifft2(scale(FI,1,1,s)));
Lyy = real(ifft2(scale(FI,0,2,s)));

G = Lxx.*Lyy-Lxy.*Lxy;
H = Lxx+Lyy;
