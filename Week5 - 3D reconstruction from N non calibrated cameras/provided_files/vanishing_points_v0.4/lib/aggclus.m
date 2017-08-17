function clus = aggclus(X, THRESHOLD);
% agglomerative clustering using single link

N = size(X,1);

%THRESHOLD = .2;

D = pdist2(X,X);

n = sqrt(X(:,1).^2 + X(:,2).^2);
n = repmat(n, [1 N]);

n = max(n,n'); % max norm

D = D./n;
D(logical(eye(N)))=inf;

clus = cell(1,N);
for i=1:N
    clus{i}=i;
end

[V, I] = min(D(:));

while V < THRESHOLD
    [i, j] = ind2sub(size(D),I);
    
    if i>j
        t = i;
        i=j;
        j = t;
    end
        
    clus{i} = [clus{i} clus{j}];
 
    clus = clus([1:j-1 j+1:end]);
    
    d = min([D(i,:); D(j,:)]);
    
    D(i,:) = d;
    D(:,i) = d';
    D(j,:) = [];
    D(:,j) = [];
    
    D(i,i)=inf;
    
    [V, I] = min(D(:));
end