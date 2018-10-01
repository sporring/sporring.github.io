function [basins, minima] = basins(f);

% BASINS   Catchment basins (influenze zones) on Morse functions, Pixelbased
%   [bas, minima] = basins(f)
%    - f  a  matrix/image
%   Returns basins and minima
%   - basins
%     All pixels in each basin is labelled the same integer number
%     All basins have distinct numbers
%   - minima
%     The minima, the lowest value locally, is labelled 255, all other pixels
%     is labelled 0
%
% Note it is assumed that each pixelvalue do not have neighbours with the 
% same value (Morse function)
%
%			Ole Fogh Olsen, February 1999

global  bestDiff bestDirection diff  % Used by update_indici

[m n] = size(f);
Size = [m n];
mn = m*n;

bestDiff = zeros(Size);
bestDirection = 1:mn;

%%% Find the lowest neighbouring pixel and put its index in dirRow and dirCol

metric = 1/sqrt(2);


% check one row down
diff = zeros(Size);
diff(1:m-1,:) = f(1:m-1, :) - f(2:m, :);
update_indici(1);

% check row above
diff(2:m, :) = -diff(1:m-1,:);
diff(1,:) = 0;
update_indici(-1);

% check row down and column right
diff(1:m-1,1:n-1) = metric * (f(1:m-1, 1:n-1) - f(2:m, 2:n));
diff(m,:) = 0;
diff(:,n) = 0;
update_indici(1+m);

% check row above and left column
diff(2:m, 2:n) = -diff(1:m-1,1:n-1);
diff(1,:) = 0;
diff(:,1) = 0;
update_indici(-m-1);

% check column right
diff(:, 1:n-1) = (f(:, 1:n-1) - f(:, 2:n));
diff(:,n) = 0;
update_indici(m);

% check left column
diff(:, 2:n) = -diff(:, 1:n-1);
diff(:,1 ) = 0;
update_indici(-m);

% check column right and above
diff(2:m, 1:n-1) = metric * (f(2:m, 1:n-1) - f(1:m-1, 2:n));
diff(1,:) = 0;
diff(:,n) = 0;
update_indici(m-1);

% check left column and row down
diff(1:m-1, 2:n) = -diff(2:m, 1:n-1);
diff(m,:) = 0;
diff(:,1) = 0;
update_indici(1-m);

% %% check borders

% bottom row 
diff = f(m,:) - f(1,:);
index = find(diff > bestDiff(m:m:mn));
index2 = index*m;
bestDiff(index2) = diff(index);
bestDirection(index2) = index*m-m+1;

% top row
diff = -diff;
index = find(diff > bestDiff(1:n:mn));
index2 = index*m-m+1;
bestDiff(index2) = diff(index);
bestDirection(index2) = index*m;

% right column
diff = (f(:,n) - f(:,1))';
index = find(diff > bestDiff(mn-m+1:mn) );
index2 = index + mn - m;
bestDiff(index2) = diff(index);
bestDirection(index2) = index;

% left column
diff = -diff;
index = find(diff > bestDiff(1:m));
index2 = index;
bestDiff(index2) = diff(index);
bestDirection(index2) = index -m + mn;

% top left corner
diff = (f(1,1) - f(m,n)) * metric;
if diff > bestDiff(1),
  bestDiff(1) = diff;
  bestDirection(1) = mn;
end;

% bottom right corner
diff = - diff;
if diff > bestDiff(mn),
  bestDiff(mn) = diff;
  bestDirection(mn) = 1;
end

% top right corner
diff = (f(1,n) - f(m,1)) * metric;
if diff > bestDiff(mn - m + 1),
  bestDiff(mn - m + 1) = diff;
  bestDirection(mn - m + 1) = m;
end  

% bottom left corner
diff = -diff;
if diff > bestDiff(m),
  bestDiff(m) = diff;
  bestDirection(m) = mn - m + 1;
end


%%% Divide into basins

% Label all minima pixels;
minima = 255 * reshape((bestDirection == (1:mn)), m, n);
index = find(minima);
basinNumbers = 1:length(index);
basins = -ones(m,n);
basins(index) = basinNumbers;


values = 1:mn;
values(index) = 0;
bestDirMatrix = sparse(1:mn, bestDirection, values, mn, mn);



while ~isempty(index),
  [index dummy] = find(bestDirMatrix(:,index));
  basins(index) = basins(bestDirection(index));
end  


% offset is the sum of 
% the row number and
% the correct multiple of column length that is m (nr of rows)
function update_indici(offset);

global bestDiff bestDirection diff % Used by update_indici

  index = find(diff > bestDiff);
  bestDiff(index) = diff(index);
  bestDirection(index) = index+offset;

  return;










