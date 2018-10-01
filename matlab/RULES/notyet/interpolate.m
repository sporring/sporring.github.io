function f = interpolate(F, X, method)
%INTERPOLATE Does an interpolation of F values
%
%        f = interpolate(F, X, method)
%         f - the interpolated value from F at X
%         F - the matrix to be interpolated
%         X - the interpolation points using column vectors
%         method - the interpolation method being either one of the
%             strings: 'nearest', 'bilinear', or 'bicubic'.
%
%       Interpolate does an interpolation of the values taken from F
%       at X.  'bicubic' has not been implemented yet.
%
%       Copyright: IBM Almaden Research Center & Jon Sporring, January 24, 1997


f=NaN*ones(1,size(X,2));
if(size(X,1) < 2)
  X = [X; ones(1,size(X,2))];
end

if strcmp(lower(method),'nearest') 
  X = round(X);
  for i = 1:size(X,2)
    f(i) = F(round(X(1,i)), round(X(2,i)));
  end
else 
  if strcmp(lower(method),'bilinear') 
    F2 = zeros(size(F)+1);
    F2(1:size(F,1),1:size(F,2)) = F;
    A = floor(X);
    B = X-A;
    for i = 1:size(X,2)
      f(i) = [1-B(1,i),B(1,i)]*F2(A(1,i)+[0:1],A(2,i)+[0:1])*[1-B(2,i);B(2,i)];
    end
  else 
    if strcmp(lower(method),'bicubic') 
      error(['The method ''' method ''' has not been implemented yet']);
    else
      error(['The method ''' method ''' is not recognized']);
    end
  end
end
