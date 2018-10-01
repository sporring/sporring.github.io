function C = mycontourc(I, v, Ix, Iy)
%MYCOUNTOURC A modification of contourc.
%
%        C = mycontourc(I, v)
%          C - A resulting list of lists of vector displacements
%          I - An image containing the isophote
%          v - An isophote value.
%          Ix - An optional x-gradient image.
%          Iy - An optional y-gradient image.
%
%       Mycontour modifies the output of contourc to remove too small
%       vector changes.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, February 7, 1997

C = contourc(I, [v,v]);

a = list_first(C);
while a > 0
  V = flipud(list_get(C,a));

  if exist('Ix') & exist('Iy') & (size(V,2) > 1)
    dx = [interpolate(Ix,V(:,1),'bilinear'),interpolate(Iy,V(:,1),'bilinear')];
    if det([dx; (V(:,2)-V(:,1))']) < 0
      V = fliplr(V);
    end
  end
  % REMOVE DUBLETS, THAT IS POINTS WHICH ARE SO CLOSE THAT THEY ARE ROUND OFF ERRORS
  V = [V, V(:,1)];
  dV = sum((V(:,2:size(V,2))-V(:,1:size(V,2)-1)).^2);
  index=1;
  for j = 2:size(V,2)
    if sum(dV(index(size(index,2)):j-1)) > 0.25
      index=[index,j];
    end
  end
  index(size(index,2)) = size(V,2);

  C = list_delete(C,a);
  C = list_insert(C,a,V(:,index),v);

  a = list_next(C,a);
end
