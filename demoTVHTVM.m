 % Demo of TVHTVM
%  close all;
 clear all;
 
 Orig = imread('babyface1.jpg');
 Orig = rgb2gray(Orig);
 I = double(Orig);
 [m,n] = size(Orig);
  
% add noise
std_n=2;    % Gaussian noise standard deviation
rng(0);
In = randn(size(I))*std_n; % White Gaussian noise
f  = I + I.^(1/2).*In;  % noisy input image

%  snrf = SNR(I,f)

 para.mu = 70;
 para.lambda = 0.6;
 para.maxiter = 500;
 para.tol = 5e-4;
 para.niter = 2;
 para.gamma = 1;
 para.ep = 0.0001;
 
 tic
 out = TVHTVM(f,I,para);
 toc
 psnr1=psnr(I,f);
iteration = out.iter;
snr2 = SNR(I,out.u);
psnru = psnr(I,out.u)

 figure;subplot(1,3,1);imshow(uint8(I));title('original');
 subplot(1,3,2);imshow(uint8(f));title('observed');
 subplot(1,3,3);imshow(uint8(out.u));
 title(sprintf('psnr %4.2f,rel err %4.2f%%',out.psnru(end),out.err(end)*100)); 
 
 
 

%%
%  %% test on real ultrasound image
% f=imread('Fig6f.png');
% f = f(:,:,1);
% f = double(f);
% [m,n] = size(f);
%  
%  para.mu = 8e+1;
%  para.lambda = 0.3;
%  para.maxiter = 200;
%  para.tol = 5e-4;
%  para.niter = 2;
%  para.gamma = 1;
%  para.ep = 0.0001;
%  
%  tic
%  out = TVHTVM1(f,para);
%  toc
%  
% iteration = out.iter;
% figure;subplot(1,2,1);imshow(uint8(f));title('observed');
% subplot(1,2,2);imshow(uint8(out.u));title('proposed result');
% 
% figure;plot(1:out.iter,out.relerr)
% title('convergence')
% xlabel('iteration number')
% ylabel('relative error')

% % % % 
% [s1,s2] = ssim(out.u,I);
% figure, imshow(s2),title([sprintf('SSIM = %.4f',s1)   ' SSIM map HTV']);
% 

 
 
 
 
 
