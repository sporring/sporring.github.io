function V = isophote(I, X)
%ISOPHOTE return a list of vectors indicating the isophote curve beginning at X
%
%       V = isophote(I, X)
%         V - A resulting list of vector displacements
%         I - An image containing the isophote
%         X - A 2-dimensional vector indicating the starting point.
%
%       Isophote makes a simple search of the image at scale s for it's
%       isophotes using sub-pixel accuracy up til a linear model of the
%       image surface.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, January 23, 1997

C = contourc(I, [min(min(I))-1,I(X(1),X(2))]);

a = 1;
Vval = abs(size(I,1)+1+size(I,2)+1);
while a < size(C,2)
  b = C(2,a)+a;
  % CONTOUR USES X,Y INSTEAD OF ROW,COL
  [minval,minpos] = min(sum(abs(C(:,a+1:b)-[X(2);X(1)]*ones(1,b-a))));
  if minval < Vval
    V = C(:,a+1:b);
    Vval = minval;
  end
%  disp([minval, minpos, C(:,a+minpos)',sum(abs(C(:,a+minpos)-[X(2);X(1)]))]);
%  disp([[a;b],C(:,a:a+2)]);
  a=b+1;
end
V = flipud(V);

% REMOVE DUBLETS, THAT IS POINTS WHICH ARE SO CLOSE THAT THEY ARE ROUND OFF ERRORS
i = 2;
while i <= size(V,2)
  if max(abs(V(:,i-1)-V(:,i))) < 0.1
    V = [V(:,1:i-1), V(:,i+1:size(V,2))];
  else
    i = i+1;
  end
end
