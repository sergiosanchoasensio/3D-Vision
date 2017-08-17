function [ E ] = gs_errfunction( parameters, observedPoints )
%GS_ERRFUNCTION computes the reprojection error

a = length(observedPoints(:))/2;
Xobs = reshape(observedPoints, 2, []);
x = Xobs(:,1:a/2);
x(end+1,:) = 1;
xp = Xobs(:,(a/2+1):end);
xp(end+1,:) = 1;

H = reshape(parameters(1:9), 3, []);
x_tilde = reshape(parameters(10:end), 2, []);
x_tilde(end+1,:) = 1;

for i=1:length(x)
    H_x = H*x_tilde(:,i);
    H_x = H_x/H_x(end,:);
    H_xp = inv(H)*xp(:,i);
    H_xp = H_xp/H_xp(end,:);
    
    aux = sum((xp(:,i) - H_x).^2) + sum((x(:,i) - H_xp).^2);
    d(i) = sqrt(aux);
end

E = d;
end
