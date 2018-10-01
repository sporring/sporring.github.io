function [F,A] = spectrum(I,t,q);
%SPECTRUM Calculate multifractal spectrum as function of exponents and scale
%
%       [F,A] = spectrum(I,t,q)
%         F - The multifractal spectrum
%         A - The corresponding scaling exponents
%         I - The image
%         t - The list of scales used (2*variance)
%         q - The list of information order used
%
%       This function evaluates the multifractal spectrum from the 
%       image.
%
%       Copyright: Jon Sporring, January 1, 1996

L = information_scale(I,t,q);
[F,A] = inf2spect(L,sqrt(12)/size(I,1),q);
