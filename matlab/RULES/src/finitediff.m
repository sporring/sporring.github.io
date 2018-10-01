function c = finitediff(d,p,h,t)
% FINITEDIFF - calculate coefficients for a linear approximation of a derivative
%
% c = finitediff(d,p,h,t)
%
%   c - a d+p vector of coefficients
%   d - the order of differentiation
%   p - the least order of approximation (rest term)
%   h - the spacing, typically 1.
%   t - the type, -1 = backwards, 0 = central, 1 = forward scheme
%
% The function solves a linear system of equations to find the coefficients
% for a finite difference approximation.
% 
% Example, the coeffients may be considered a kernel for convlution if
% flipped, i.e., for a square function x, the approximate central
% derivatives to first order are 
%
%   x = [-5:5].^2;
%   K1 = finitediff(1,1,1,0)
%   K2 = finitediff(2,1,1,0)
%   K3 = finitediff(3,1,1,0)
%   conv(x,fliplr(K1),'valid')
%   conv(x,fliplr(K2),'valid')
%   conv(x,fliplr(K3),'valid')
%
% October 27, 2014: Output transposed, documentation improved and a minor
% modification in the calculation.  
%
% Copyright 27/11-2003, Jon Sporring, DIKU
  
  if t < 0
    % Backward
    i = 1-d-p:0;
  elseif t>0 
    % Forward
    i = 0:d+p-1;
  else
    % Central
    if rem(d+p,2) == 0
      p = p+1;
    end
    i = -(d+p-1)/2:(d+p-1)/2;
  end
  
  n = 0:d+p-1;
  s = zeros(length(i)); 
  a = zeros(size(n))'; 
  for t = 1:length(n)
    s(t,:) = i.^(n(t));
    a(t) = (n(t)==d);
  end
  c = ((prod(1:d)/(h^d))*(s\a))';
