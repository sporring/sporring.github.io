%function J = random_image(J,I,n)

% I and J are assumed quadratic!
if (size(I,1) ~= size(I,2)) & (size(J,1) ~= size(J,2))
  error('Use only quadratic images for I and J!');
end

N = size(J,1);
M = size(I,1);

if n > 0.7*N^2/M^2
  error('Too many copies.  No room for all non-overlapping copies!');
end

X = rand(2,1)*(N-M-1)+round(M/2)*ones(2,1);

for i = 1:n
  still_space = 0;
  while still_space < 100
    x = round(rand(2,1)*(N-M-1)+1*ones(2,1));
    if sqrt(min(sum((X-x*ones(1,size(X,2))).^2))) > sqrt(2)*M+1
      break;
    end
    still_space = still_space+1;
  end
  if still_space < 100
    X = [X,x];
    J([1:M]+x(1),[1:M]+x(2)) = I;
  else
    break;
  end
end
