
%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear all
clc;
K=imread('3.bmp');
%K=rgb2gray(K);
K=double(K);
I=K;

%%%%%%%%%%%% Add noise
std_n=3; % Gaussian noise standard deviation
In=normrnd(0,std_n,(size(I)));
%In = randn(size(I))*std_n; % White Gaussian noise
I0 = I + I.^(1/2).*In;  % noisy input image
% show original and noisy images
close all
figure(1); imshow(uint8(I)); title('Original')
figure(2); imshow(uint8(I0)); title('Noisy image')
imwrite(uint8(I0),'Noisy image.bmp');
SNR1=SNR(K,I0);
ReErr1=ReErr(K,I0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  TV denoising with (scalar) data fidelity term  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Here we assume the noise variance is known.
%% Therefore we can compute lambda (instead of fixing
%% iter). TV therfore can run automatically quite well.
%% The process is slower.
J=I0; 
% params
ep_J =0.05; % minimum mean change in image J
lam=0.2; %J_old=0;
std=0;
iter=1; dt=0.2; ep=0.0001; 
var_n = std_n^2; % noise variance
k=0;
tic;
while (abs(std-std_n) > ep_J),  % iterate until convergence
%for i=1:350
  J=fof(J,iter,dt,lam,ep,I0);  
  J =max(ep+0.5,min(J,255)); % scalar 
   lam = clam(J,I0,var_n,ep);
  N=((I0-J).^2)./J;
  std=sqrt(mean(mean(N)));% update lambda (fidelity term) 
   k=k+1;
end % for i

  % iterate until convergence
% for k=1:308
%   J=fof(J,iter,dt,lam,ep,I0);  
%   J =max(ep+0.5,min(J,255)); % scalar 
%    lam = clam(J,I0,var_n,ep);
%   N=((I0-J).^2)./J;
%   std=sqrt(mean(mean(N)));% update lambda (fidelity term) 
%    k=k+1;
% end % for i


toc;
iteration=k
figure(3); imshow(uint8(J)); title('Denoised image with lambda')
imwrite(uint8(J),'TV filter with lambda.bmp');
%compute the SNR
SNR2=SNR(K,J)
ReErr2=ReErr(K,J)
%%%%%%%%%%%%%%%%%%

