function zc = zerocrossing(Y)
%ZEROCROSSING evaluate the fractional index of the zero crossings of a vector
%
%       zc = zerocrossing(Y);
%         zc - a list of fractional indices where zero crossings occur in Y
%         Y - a 1D real vector.
%
%       Zerocrossing assumes a the x-component to be the indices of Y.

zc = [];
if min(size(Y)) > 1
  error('Zerrocrossing has only been implemented for 1D functions.');
else
  if size(Y,1) > 1
    Y = Y';
  end
  ind = find((sign(Y(1:max(size(Y))-1)) ~= sign(Y(2:max(size(Y))))));
  zc = zeros(1,size(ind,2));
  a = Y(ind+1)-Y(ind);
  b = Y(ind)-a.*ind;
  zc = -b./a;
  zc = [zc,find(Y==0)];
end
