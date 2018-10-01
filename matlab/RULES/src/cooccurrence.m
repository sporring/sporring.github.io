function J = cooccurrence(I,N)
% COOCCURRENCE  Calculate cooccurrence matrix for a gray valued image.
%
%       J = cooccurrence(I,N)
%         I - A gray valued image
%         N - An nx2 matrix of spatial vector displacements in row column
%             notation.
%
%       Cooccurrences is a statistical method for analysing texture images.
%       The present implementation can evaluate the n-dimensional
%       cooccurrence image for a 2 dimensional image.  The spatial
%       dependency is user supplied.  Be warned, a gray valued image has 256
%       gray values, hence the resulting cooccurence image will be 256^n.
%
%       Example:
%         % load image and calculate cooccurrence matrix
%         I = dikuread('/usr/local/image/archive/DIKU-Robot/eye-hand2.im');
%         J = cooccurrence(I,[0,1;1,0]);
%         % display cooccurrences between (i,j+1) and (i+1,j)
%         figure;
%         colormap(gray(256));
%         imagesc(squeeze(sum(J,1)));
%         % display cooccurrences between (i,j) and (i+1,j)
%         figure;
%         colormap(gray(256));
%         imagesc(squeeze(sum(J,2)));
%         % display cooccurrences between (i,j) and (i,j+1)
%         figure;
%         colormap(gray(256));
%         imagesc(squeeze(sum(J,3)));


% J should have been sparse, but unfortunately this only works for 2D
% matrices (with Matlab 5.2).

if (min(min(I)) < 1) | (max(max(I)) > 256)
  warning('Image intensities must be between 1 and 256.  Rescaling!');
  Imin = min(min(I));
  Imax = max(max(I));
  I = (I-Imin)*255/(Imax-Imin)+1;
end
maxI = 256;
J = zeros(maxI*(ones(1,size(N,1)+1)));
Jsize = size(J);
for i=1:size(I,1)-max(N(:,1))
  for j=1:size(I,2)-max(N(:,2))
    v = I(i,j);
    for k=1:size(N,1)
      v=v+prod(Jsize(1:k))*(I(i+N(k,1),j+N(k,2))-1);
    end
    J(v)=J(v)+1;
  end
end

% That was the general code for this specialized example
% for i=1:size(I,1)-max(N(1,:))
%   for j=1:size(I,2)-max(N(2,:))
%     for k=1:size(N,2);
%       J(I(i,j),I(i+N(1,k),j+N(2,k)))=J(I(i,j),I(i+N(1,k),j+N(2,k)))+1;
%     end
%   end
% end


