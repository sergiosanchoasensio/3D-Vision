function Xtrain = triangulate (x1, x2, P1, P2, imsize)

H = [ 2/imsize(1)     0       -1
         0       2/imsize(2)  -1
         0            0        1 ];

x11=euclid(H*homog(x1));
x22=euclid(H*homog(x2));
P11=H*P1;
P22=H*P2;

A = [x11(1)*P11(3,:)-P11(1,:);x11(2)*P11(3,:)-P11(2,:); x22(1)*P22(3,:)-P22(1,:);x22(2)*P22(3,:)-P22(2,:)];

Aprim = [A(:,1) A(:,2) A(:,3)];

[U,D,V] = svd(Aprim);
UT=transpose(U);
A4prim=UT*A(:,4);
A4prim = homog(A4prim);

for i=1:size(D,2)
    Y(i,1)= (-A4prim(i))/D(i,i);
end

Xtrain = V*Y;
Xtrain=homog(Xtrain);

end