I = double(imread('noise.tif'));
colormap(gray(256));

help eed

disp('eed(I)');
J = eed(I);
subplot(3,1,1); imagesc(I); subplot(3,1,2); imagesc(J); subplot(3,1,3); imagesc(abs(J-I))

input('Press <enter>');
disp('eed(I,10)');
J = eed(I,10);
subplot(3,1,1); imagesc(I); subplot(3,1,2); imagesc(J); subplot(3,1,3); imagesc(abs(J-I))

input('Press <enter>');
disp('eed(I,10,1)');
J = eed(I,10,0.1);
subplot(3,1,1); imagesc(I); subplot(3,1,2); imagesc(J); subplot(3,1,3); imagesc(abs(J-I))

input('Press <enter>');
disp('eed(I,10,1,1)');
J = eed(I,10,0.1,1);
subplot(3,1,1); imagesc(I); subplot(3,1,2); imagesc(J); subplot(3,1,3); imagesc(abs(J-I))

input('Press <enter>');
disp('eed(I,10,1,7,5)');
J = eed(I,10,0.1,3,5);
subplot(3,1,1); imagesc(I); subplot(3,1,2); imagesc(J); subplot(3,1,3); imagesc(abs(J-I))

