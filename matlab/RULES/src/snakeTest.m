% A small program for testing snake.m.  The parameters alpha, beta,
% gamma, sigma, N, and I are all meant to be played with!
%
%                     Copyright, Jon Sporring, DIKU, December 2003.

% Constants of the method
alpha = 1000;
beta = 10;
gamma = 1000;
sigma = 5;
N = 100;
I = zeros(256,256);
I(128+[-64:64],128+[-64:64]) = 256;
I(64+[-32:32],128+[-32:32]) = 0;
%I = imread('054.tif');

% Prepare the Image term.  The function scale is a convolution with a
% Gaussian and may be found on my homepage.
Is = real(ifft2(scale(fft2(I),sigma,0,0)));
[Isc,Isr] = gradient(Is);
G = sigma^2*(Isc.^2+Isr.^2);

% Make an initial curve
t = linspace(0,2*pi,N)';
t = t(1:end-1);
x = [size(I,1)/2+size(I,1)*7/16*cos(t), size(I,1)/2+size(I,1)*7/16*sin(t)];

% Setup the figure
f = figure(1);
clf
subplot(1,2,1)
imagesc(I);
colormap(gray);
title('Original');
subplot(1,2,2)
imagesc(G);
colormap(gray);
title('Image Energy');
set(f,'DoubleBuffer','on')
h1 = -1;
h2 = -1;
for t = 1:1000
  x = snake(x,alpha,beta,G,gamma,10);
  
  if h1 > -1
    delete(h1);
  end
  if h2 > -1
    delete(h2);
  end
  subplot(1,2,1);
  hold on;
  h1 = plot([x(:,2);x(1,2)],[x(:,1);x(1,1)],'linewidth',sigma);
  hold off;
  subplot(1,2,2);
  hold on;
  h2 = plot([x(:,2);x(1,2)],[x(:,1);x(1,1)],'linewidth',sigma);
  hold off;
  drawnow;
end

