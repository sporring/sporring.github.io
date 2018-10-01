% A test program for em.m
%   Jon Sporring, DIKU, 17/11-2003

% Choose image to try: 1 - vietnam; 2 - dug; 3 - car; 4 - synthetic
IMAGE = 4;

switch(IMAGE)
  case 1,
   disp('Vietnam multispectral image');
   examplebands = [3,4,5];
   downsample = 4;
   tmp = imresize(imread(['vietnam_1.tif']),1/downsample);
   Iorig = zeros([size(tmp,1),size(tmp,2),6]);
   Iorig(:,:,1) = tmp;
   Iorig(:,:,2) = imresize(imread('vietnam_2.tif'),1/downsample);
   Iorig(:,:,3) = imresize(imread('vietnam_3.tif'),1/downsample);
   Iorig(:,:,4) = imresize(imread('vietnam_4.tif'),1/downsample);
   Iorig(:,:,5) = imresize(imread('vietnam_5.tif'),1/downsample);
   Iorig(:,:,6) = imresize(imread('vietnam_7.tif'),1/downsample);
 case 2,
  disp('Dug color image');
  examplebands = [1,2,3];
  downsample = 4;
  Iorig = imresize(double(imread('dug4.jpg')),1/downsample);
 case 3,
  disp('Car color image');
  examplebands = [1,2,3];
  downsample = 2;
  Iorig = imresize(double(imread('car1.jpg')),1/downsample);
 case 4,
  disp('Synthetic color image');
  examplebands = [1,2,3];
  Iorig = 255*[repmat(permute([0.75,0.25,0.25],[1,3,2]),[50,25,1]),repmat(permute([0.5,0.5,0.25],[1,3,2]),[50,25,1])];
  Iorig = Iorig+[10*randn(50,25,3),5*randn(50,25,3)];
 otherwise,
  error('This image will be added in the future!');
end

% Show the examplebands of the image;
figure(1)
clf
image(Iorig(:,:,examplebands)/255);
axis image;
drawnow;

% Reshape the data into NxD data and prepare code-minimization loop
I = reshape(Iorig,[size(Iorig,1)*size(Iorig,2),size(Iorig,3)]);
f = figure(2);
clf;
set(f,'DoubleBuffer','on')
min_l = inf;
min_mu = 0;
min_C = 0;
lk = zeros(20,1);
for k = 1:length(lk)
  % Initialize model by random examples
  mu = zeros(k,size(I,2));
  for i = 1:size(mu,1)
    mu(i,:) = I(round(rand*(size(I,1)-1)+1),:);
  end
  C = std(I(:))^2*eye(size(mu,2));
  C = repmat(C(:)',[size(mu,1),1]);
  
  E = 0;
  i = 1;
  dE = 1;
  muold = zeros(size(mu));
  Prior = ones(size(mu,1),1)/size(mu,1);
  while(dE > 0.0001)
    muold = mu;

    % Performe an EM step
    [mu,C,Posteriori,Likelihood,Prior,l] = emGaussianMixture(I,mu,C,Prior);
    E(i) = l;
    fprintf(1,'%d klasser, %g bits per pixel per band         \r',k,l);
    if i < 2
      dE = 1;
    else
      dE = E(i-1)-E(i);
    end
    i=i+1;
    
    % Classify the pixels according to the maximum posteriori
    [val,ind] = max(Posteriori,[],2);

    % Visualize the classifications by the colors of the example bands
    figure(f)
    image(reshape(mu(ind,examplebands)/255,[size(Iorig,1),size(Iorig,2),length(examplebands)])); 
    axis image;
    title([num2str(k),' klasser']);
    drawnow;
  end
  fprintf(1,'\n');
  lk(k) = l;

  % Store the resutlts, if this number of classes is best seen.
  if l < min_l
    min_l = l;
    min_mu = mu;
    min_C = C;
    min_Posteriori = Posteriori;
  end
end

% We're done, visualize the optimal results.
disp(min_l)
disp(min_mu)
disp(min_C)
[val,ind] = max(min_Posteriori,[],2);
image(reshape(min_mu(ind,examplebands)/255,[size(Iorig,1),size(Iorig,2),length(examplebands)])); 
axis image;
title([num2str(size(min_mu,1)),' klasser']);

figure(3);
plot(lk);
xlabel('Number of Classes');
ylabel('Coding Cost Per Pixel Per Band');
