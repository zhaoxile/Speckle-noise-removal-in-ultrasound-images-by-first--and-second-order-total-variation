function SNR=SNR(K,J)
% J:restored image K:original image
N=J-K;
L=(J-mean(mean(J))).^2;
X=(N-mean(mean(N))).^2;
a=sum(sum(L));
b=sum(sum(X));
SNR=10*log10(a/b);