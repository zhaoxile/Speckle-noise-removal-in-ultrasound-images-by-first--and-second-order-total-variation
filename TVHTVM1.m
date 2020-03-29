function out = TVHTVM1(f,para)
%  write by si 01/19/2015
%  solve TV+HTV multiplicative noise remove based on “Fast reduction of speckle noise in real ultrasound images 2013”
%  solve the model min_{z} \alpha||z||_{TV}+(1-\alpha)||z||_{HTV}+.......
%    ......................\lambda(f.*e^{-z/2}log(f/e^{z})-f.*e^{-z/2}+e^{z/2})
%  min_{x} \alpha||Dx||+(1-\alpha)||D^{2}x||+<\lambda,f.*e^{-x/2}log(f/e^{x})-f.*e^{-x/2}+e^{x/2}>
%   s.t. x=z, p=Dz, y=D^{2}z;

%  u=e^{z}
% update \alpha through \alpha  =   1 ,          if ||z||_{TV}>=c   c:constant
%                               =  1/2cos(2*pi*|Du|/c)+1/2 .   elsewhere 
% reference “Fast reduction of speckle noise in real ultrasound images 2013 Jiehuang Signal Processing”
%            “Iterative image restoration combining total variation minimization and second-order functional 2006 M.Lysaker”


[m,n] = size(f);
mu = para.mu;
lambda = para.lambda;
niter = para.niter;
tol = para.tol;
maxiter = para.maxiter;
gamma = para.gamma;
ep = para.ep;
% z = zeros(m,n);
% x = zeros(m,n);
f = regular(abs(f),ep);   % very important,convex model (f>0)
u = f;
z = log(f);  %  保证log有意义，不出现inf
x = zeros(m,n);
% x = z;
p = zeros(m,n,2);
q = zeros(m,n,4);
b1 = zeros(m,n);
b2 = zeros(m,n,2);
b3 = zeros(m,n,4);
cont = 1;
k = 0;

% out.psnru = [];
out.relerr = [];


eigsDtD = abs(psf2otf([1,-1],[m,n])).^2 + abs(psf2otf([1;-1],[m,n])).^2;
eigsD2tD2 = abs(psf2otf([1,-2,1],[m,n])).^2 +abs(psf2otf([1;-2;1],[m,n])).^2 + 2*abs(psf2otf([1,-1; -1,1],[m,n])).^2;
% DE = mu(2)*eigsDtD + mu(3)*eigsD2tD2; 
DE = eigsDtD +eigsD2tD2; 

H=fspecial('gaussian',5,2);
ff=imfilter(f,H,'circular','conv');
[fx,fy]=grad(ff);
ftemp=fx.^2+fy.^2;
eta=0.0001;
kk=0.01;
alpha=(eta+kk*ftemp)./(1+eta+kk*ftemp);
beta=1-alpha;

h = waitbar(0,'Please wait ...');

while cont
    k = k+1;
    % solve z-subproblem
    uold = u;
    tt1 = p-b2;
    tt2 = q-b3;
%     Rigz = mu(1)*(x-b1)+mu(2)*(div(tt1(:,:,1),tt1(:,:,2)))+mu(3)*(div2(tt2(:,:,1),tt2(:,:,2),tt2(:,:,3),tt2(:,:,4)));
%     z = ifft2(fft2(Rigz)./(mu(1)+DE));
    Rigz = (x-b1)+(div(tt1(:,:,1),tt1(:,:,2)))+(div2(tt2(:,:,1),tt2(:,:,2),tt2(:,:,3),tt2(:,:,4)));
    z = ifft2(fft2(Rigz)./(1+DE));
    z = real(z);
    u = real(exp(z));
    u(u<=0)=ep;u(u>255)=255;

    
    % solve x-subproblem by newton method
    xold = x;
    x = Multinewton(xold,mu(1),lambda,f,b1,niter);
    
    
    % sove p-subproblem
    [zx,zy] = grad(z);
    Dz(:,:,1) = zx;
    Dz(:,:,2) = zy;
    p = shrink2s(Dz+b2,alpha/mu);
    
    
    % solve q-subproblem
    DDz =grad2(z);
    q = shrink2s(DDz+b3,beta/mu);
    
    % update b
    b1 = b1-gamma*(x-z);
    b2 = b2-gamma*(p-Dz);
    b3 = b3-gamma*(q-DDz);
    
    
    relerr = norm(u(:)-uold(:))/norm(uold(:));
    % stopping criterion
    if relerr<=tol || k>maxiter
        cont = 0;
    end 
    out.relerr = [out.relerr,relerr];
    
    
%     % update alpha
%     Mz = sqrt(zx.^2+zy.^2);
%     Mz = Mz./max(Mz(:));   %  scale
%     alpha = 0.5*cos(2*pi*Mz./c)+0.5;
%     alpha(Mz>=c) = 1;
    
 
    % Display the progress
    waitbar(k/maxiter,h,sprintf('iteration %d',k));
    
end
close(h);


out.u = u;
out.iter = k;


        
        
  
  % solve x-subproblem by newton method
% x = arg min_{x} <lambda,f*exp(-x/2)log(f/exp(x))-f*exp(-x/2)+exp(x/2)>+mu/2||x-z-b1||_{2}^{2}
% ..........+mu(2)/2||D^{2}z+b2-y||_{2}^{2}


    function  x = Multinewton(xp,mux,lambda,f,b1,niter)
        kn=0;
        cont1 = 1;
        while cont1
            kn = kn+1;
            temp1 = f.*exp(-xp/2);
            temp2 = exp(xp/2);
            dF = lambda*(-0.5*temp1.*(log(f)-xp+1)+0.5*temp2)+mux*(xp-z-b1);
            ddF = lambda/4*temp1.*(exp(xp)./f+log(f)+3-xp)+mux;
            x = xp-dF./ddF;
          % re_inim = norm(x(:)-xp(:))/norm(xp(:));
            xp = x;
            cont1 = (kn<niter); 
        end
    end
    
  
  function [Ux,Uy] = grad(U)
            %                   ux = Dxplus(u);
            %                   uy = Dyplus(u);
            Ux = [diff(U,1,2)  U(:,1) - U(:,end)];
            Uy = [diff(U,1,1); U(1,:) - U(end,:)];
  end

 
 function DtXY = div(X,Y)
         % Transpose of the forward finite difference operator
         % div = - grad^*
            DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
            DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];
 end
 
 
     function DDu = grad2(U)

       Ux = [diff(U,1,2)  U(:,1) - U(:,end)];
       Uy = [diff(U,1,1); U(1,:) - U(end,:)];
       DDu(:,:,1) = [Ux(:,1) - Ux(:,end) diff(Ux,1,2)];
       DDu(:,:,2) = [diff(Ux,1,1); Ux(1,:) - Ux(end,:)]; 
       DDu(:,:,3) = DDu(:,:,2);
       DDu(:,:,4) = [Uy(1,:) - Uy(end,:); diff(Uy,1,1)];

     end

    function D2tX4= div2(X1,X2,X3,X4)
%         D2tX4 = Dxplus(Dxminus(X1)) + Dyminus(Dxminus(X3)) +Dxminus(Dyminus(X2))+Dyplus(Dyminus(X4));
        DtX4= [X1(:,2) + X1(:,end) - 2*X1(:,1), diff(X1,2,2), X1(:,1) + X1(:,end-1) - 2*X1(:,end)]; % Dxx
        DtX4 = DtX4 + ...
                    [X4(2,:) + X4(end,:) - 2*X4(1,:); diff(X4,2,1); X4(1,:) + X4(end-1,:) - 2*X4(end,:)];  % Dyy
        Temp = [X2(end,:) - X2(1,:); -diff(X3,1,1)];
        D2tX4 = DtX4 + 2*[Temp(:,end) - Temp(:,1), -diff(Temp,1,2)];                                               % Dxy
    end
    
end


