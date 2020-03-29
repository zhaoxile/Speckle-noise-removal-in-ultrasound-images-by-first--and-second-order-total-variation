function reg=regular(M,ep)
[ny,nx]=size(M);
for i=1:ny
        for j=1:nx
           if (abs(M(i,j))<=0)
                M(i,j)=ep;
            else
                M(i,j)=M(i,j);
            end
        end
end
 reg=M;   