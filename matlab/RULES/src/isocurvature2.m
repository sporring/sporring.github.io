function C = isocurvature2(FI,s);
%ISOCURVATURE2 Calculate the isophote curvature.
%       
%       C = isocurvature2(FI,s)
%         C - the curvature image
%         FI - the Fourier Transform of the original image
%         s - the scale at which the curvature is calculated
%
%       Isocurvature calculates the isophote curvature using image
%       using Scale-Space image derivatives including order 2, and
%       the scale should be set accordingly.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, January 13, 1997

C = isocurvature(real(ifft2(scale(FI,s,1,0))),real(ifft2(scale(FI,s,0,1))),real(ifft2(scale(FI,s,2,0))),real(ifft2(scale(FI,s,1,1))),real(ifft2(scale(FI,s,0,2))));
