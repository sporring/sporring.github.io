function borders = boundaries(f);

% borders = boundaries(f);

[m n] = size(f);
 
RowDown = (f(1:m-1,:) ~= f(2:m,:));
RowDown(m,:) = zeros(1,n);
ColumnRight = (f(:,1:n-1) ~= f(:,2:n));
ColumnRight(:,n) = zeros(m,1);

borders =  RowDown | ColumnRight | RowDown([m 1:m-1],:) | ColumnRight(:,[n 1:n-1]);

borders = 255 * borders;
