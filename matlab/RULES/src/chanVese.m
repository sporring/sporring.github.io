function phi = tv;
% Total Variation with Curvature regularization!
%
%                     Copyright, Jon Sporring, DIKU, December 2003.
  
% % BW example
% I = C2*ones(128,128);
% I(64+[-32:32],64+[-32:32]) = 50;
% I(32+[-16:16],64+[-16:16]) = 0;
% I = fliplr(I');
% % Constants of the method
% sigma = 10;
% mu = 0;
% nu = 0;
% lambda1 = 1;
% lambda2 = 1;
% C1 = 50;
% C2 = 0;
% gamma = 1000;
% N = 10;

  % MR image example
  I = imresize(double(imread('054.tif')),[256,256]);
  % Constants of the method
  sigma = 5;
  mu = 0.01;
  nu = 0;
  lambda1 = 1;
  lambda2 = 1;
  C1 = 150;
  C2 = 50;
  gamma = 100;
  N = 100;

  % phi = sgnDstFromImg(I>=C1);
  phi = sgnDstFromImg(real(ifft2(scale(fft2(rand(size(I))> 0.5),5,0,0)))>0.5);
  lambda1 = lambda1/prod(size(I));
  lambda2 = lambda2/prod(size(I));
  
  % Setup the figure
  f = figure(1);
  clf
  set(f,'DoubleBuffer','on')
  drawit(I,phi,0);
  
  time = 0;
  mu2 = 0;
  nu2 = 0;
  for i = 1:gamma*N
    % Convergens er hurtigere, hvis de første iterationer lader mu=nu=0
    if i == 100;
      disp('Now using possible non-zero values for nu and mu');
      mu2=mu;
      nu2=nu;
    end
    
    disp([time,gamma]);
    
    phix = (phi([2:end,1],:)-phi([end,1:end-1],:))/2;
    phiy = (phi(:,[2:end,1])-phi(:,[end,1:end-1]))/2;

    phixx = (phi([2:end,1],:)-2*phi+phi([end,1:end-1],:));
    phiyy = (phi(:,[2:end,1])-2*phi+phi(:,[end,1:end-1]));
    phixy = (phix(:,[2:end,1])-phix(:,[end,1:end-1]))/2;
    
    phix2 = phix.^2;
    phiy2 = phiy.^2;
    g2 = phix2+phiy2;
    
    % Introduce randomness at zero gradient.
    ind = find(g2 < 1000*eps);
    phix(ind) = 2*rand(size(ind))-1;
    phiy(ind) = 2*rand(size(ind))-1;
    phix2(ind) = phix(ind).^2;
    phiy2(ind) = phiy(ind).^2;
    g2(ind) = phix2(ind)+phiy2(ind);
    
    k = (phixx.*phiy2+phiyy.*phix2-2*phixy.*phix.*phiy)./g2.^(3/2);
    k(ind) = 0;
    
    phit = exp(-phi.^2/(2*sigma^2)).*(mu2*k-nu2-lambda1*(I-C1).^2+lambda2*(I-C2).^2);
    % phit = sqrt(g2).*(mu2*k-nu2-lambda1*(I-C1).^2+lambda2*(I-C2).^2);
    gamma = 1.1*max(abs(phit(:)));
    phi = phi+(1/gamma)*phit;
    time = time+1/gamma;
    
    if rem(i,10) == 0 
      disp('Reinitializing');
      phi = sgnDstFromImg(phi<0);
    
      drawit(I,phi,phit);
    end
  end
  
function drawit(I,phi,phit)

  disp('drawing');

  subplot(2,2,3);
  imagesc(phit);
  axis equal
  axis tight
  colorbar;
    
  subplot(2,2,2);
  imagesc(phi);
  axis equal
  axis tight
  colorbar;
    
  subplot(2,2,1);
  imagesc(I);
  axis equal
  axis tight
  colormap(gray);
  colorbar;
  hold on;
  [C,h] = contour(phi,[0,0],'r');
  hold off;
  drawnow;
