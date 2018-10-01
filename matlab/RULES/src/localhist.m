function [H,vs] = localhist(I,row,col,sigma,beta,alpha,vs)
%LOCALHIST  Local soft histogram in a point of an image.
%
%       [H,vs] = localhist(I,row,col,sigma,beta,alpha,vs)
%         H - the resulting histogram
%         vs - the intensity values H is sampled in. If not given,
%             then default values will be chosen
%         I - the original an 2d image
%         row,col - The coordinates where the histogram is calculated.
%         sigma - the standard deviation of the Gaussian to smooth the
%             image with.
%         beta - the standard deviation of the Gaussian to soften the
%             isophotes with.
%         alpha - the standard deviation of the Gaussian window
%             within which the histogram is calcualted.
%
%       This function returns the local, normalized histogram in a point of an
%       image for a specific setting of sigma, beta, and alpha.  Hence,
%       the output H is one dimensional.  The number of bins in the
%       histogram is set equal to ceil((max(I(:))-min(I(:)))/max(eps,beta)),
%       and H is normalized such that sum(H*beta)=1.
%
%       Example:
%         I = 16*rand(64,64);
%         [H,vs] = localhist(I,32,32,2,1,16);
%         plot(H);
%
%       Copyright: Jon Sporring, December 1, 1999
%
%       History:
%       2009/08/13, sporring: Corrected histogram offset error
%       2011/05/15, sporring: Normalized histogram, set adaptive
%         intensity sampling values, added intensity values to be
%         returned, and rewritten histogram sampling procedure

L = real(ifft2(scale(fft2(I),sigma,0,0)));
%beta=max(0.01,beta);

minV = min(min(I));
maxV = max(max(I));
if nargin < 7
  vs = linspace(minV,maxV,ceil((maxV-minV)/beta));
end
H = zeros(size(vs));
if alpha == 0
  [val,ind] = min((vs-L(row,col)).^2);
  H(ind) = 1;
elseif alpha == inf
  for i = 1:length(vs);
    J = exp(-(L-vs(i)).^2/(2*beta^2));
    H(i) = sum(sum(J));
  end
else
  rows = (([1:size(I,1)]-row)'*ones(1,size(I,2)));
  cols = ones(size(I,1),1)*([1:size(I,2)]-col);
  A = exp(-(rows.^2+cols.^2)/(eps+2*alpha^2));
  A = A/sum(sum(A));

  for i = 1:length(vs);
    J = exp(-(L-vs(i)).^2/(2*beta^2));
    H(i) = sum(sum(J.*A));
  end
end
H=H/sum(H)/beta;
