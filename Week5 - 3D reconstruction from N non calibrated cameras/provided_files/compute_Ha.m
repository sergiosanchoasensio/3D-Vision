function [ Ha ] = compute_Ha(Pproj, Hp, u, v, z)
% ToDo: compute the matrix Ha that 
%       upgrades the projective reconstruction to an affine reconstruction
% Use the following vanishing points given by three pair of orthogonal lines
% and assume that the skew factor is zero and that pixels are square

    % From lecture 9, page 35/40
    Amat = ...
    [u(1)*v(1), u(1)*v(2)+u(2)*v(1), u(1)*v(3)+u(3)*v(1), u(2)*v(2), u(2)*v(3)+u(3)*v(2), u(3)*v(3);
     u(1)*z(1), u(1)*z(2)+u(2)*z(1), u(1)*z(3)+u(3)*z(1), u(2)*z(2), u(2)*z(3)+u(3)*z(2), u(3)*z(3);
     v(1)*z(1), v(1)*z(2)+v(2)*z(1), v(1)*z(3)+v(3)*z(1), v(2)*z(2), v(2)*z(3)+v(3)*z(2), v(3)*z(3);
     0        , 1                  , 0                  , 0        , 0                  , 0        ;
     1        , 0                  , 0                  , -1       , 0                  , 0        ];

    [~,~,V] = svd(Amat);
    w = V(:,end);

    W = [w(1), w(2), w(3);
         w(2), w(4), w(5);
         w(3), w(5), w(6)];

    % P[M|m]
    P = Pproj(1:3,:)*inv(Hp);
    M = P(:,1:3);

    B = inv(M'*W*M)
    A = chol(B);

    Ha = [inv(A'), [0;0;0];
          0, 0, 0, 1];
end

