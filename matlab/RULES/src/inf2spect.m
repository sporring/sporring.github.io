function [F,A] = inf2spect(L,l,q);
%INF2SPECT Calculate multifractal spectrum as function of exponents and scale
%
%       [F,A] = inf2spect(L,l,q)
%         F - The multifractal spectrum
%         A - The corresponding scaling exponents
%         L - The generalized entropy (from e.g. information_scale)
%         l - The list of lengths used in the measurement usually (const. 1/N)
%         q - The list of information order used
%
%       This function evaluates the multifractal spectrum from the 
%       generalized entropy
%
%       Copyright: Jon Sporring, January 1, 1996

% The generalized dimensions
D = -L./log(l);

% The scaling exponents
T = ((q'-1)*ones(1,size(L,1)))'.*D;
A = (T(:,3:size(q,2))-T(:,1:size(q,2)-2))/(q(3)-q(1)); 
dq = q(2:size(q,2)-1);

% The spectrum
F = (dq'*ones(1,size(L,1)))'.*A - T(:,2:size(T,2)-1);
