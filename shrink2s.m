function y = shrink2s(x,thresh)
%% Isotropic shrinkage operator in R^n spaces
% shrink2(x,thresh) returns 
%      $ y = max(||x||-thresh,0).*x/||x|| $ if $x\neq 0$
%      $ y = 0 $ otherwise.
% Inputs: 
% x - M * N * k matrix
% thresh - thresholding parameter >0.
%
% Outputs:
% y - M * N * k matrix
%
% Si 07/10/2014

k = size(x,3);
r = sqrt(sum(x.^2,3));
y = zeros(size(x));
for i = 1:k
    y(:,:,i) = x(:,:,i).*(1-thresh./r);
end
y(isnan(y)) = 0;

