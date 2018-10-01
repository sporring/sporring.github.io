function L = local_information(I,s,q);
%LOCAL_INFORMATION Calculate the local entropy map of an image
%
%       L = local_entropy(I,s,q)
%         L - The local information map
%         I - The image
%         s - The window scale (2*variance)
%         q - The information order
%
%       This function evaluates the local information in each point.  Locallity
%       is defined as a multiplication with a gauss function of variance
%       s/2.  Each locallity is normalized before the information contents
%       is evaluated. Be warned: This is a slow function!
%
%       Copyright: Jon Sporring, January 1, 1996

% Local neighbourhood function
rows = floor(3*sqrt(s/2));
if rows > size(I,1)
  error(sprintf('The window function is too large!\nUse s < %6.2f', 2*(size(I,1)/3)^2));
end
cols = rows;
G = exp(-(([-rows:rows]'*ones(1,2*cols+1)).^2 + (ones(2*rows+1,1)*[-cols:cols]).^2)/s);
G = G/sum(sum(G));

% Calculate the local generalized entropy
L = zeros(size(I));
I3 = extend(I,3,3);
I3 = fliplr(flipud(I3(size(I,1)-rows-1+[1:size(I,1)+2*rows+1],size(I,2)-cols-1+[1:size(I,2)+2*cols+1])));
for x = 1:size(L,1)
  for y = 1:size(L,2)
    I4 = I3(x+[0:size(G,1)-1],y+[0:size(G,2)-1]).*G;
    L(x,y) = information(I4/sum(sum(I4)),q);
  end
end
