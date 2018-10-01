function J = extend(I,M,N);
%EXTEND Periodic flip-extension of matrices
%       
%       J = extend(I,M,N)
%         J - the extended matrix
%         I - the original matrix
%         M - number of periods in row direction
%         N - number of periods in column direction
%
%       Extend makes a MxN inverted periodic extension of the matrix
%       I with origional as 1,1.  This is especially usefull when
%       implemention scale space with mirroring boundaries.
%
%       Copyright: Jon Sporring, September 3, 1996

Ilr = fliplr(I);
Iud = flipud(I);
Ilrud = fliplr(Iud);

J = zeros(M*size(I,1),N*size(I,2));
for n = 1:N
  for m = 1:M
    r = (m-1)*size(I,1)+[1:size(I,1)];
    c = (n-1)*size(I,2)+[1:size(I,2)];
    
    if rem(n,2) & rem(m,2)
      J(r,c) = I;
    else
      if rem(n,2)
	J(r,c) = Iud;
      else
	if rem(m,2)
	  J(r,c) = Ilr;
	else
	  J(r,c) = Ilrud;
	end
      end
    end
  end
end
