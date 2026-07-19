function U = pclda(data,labels,percent)

ngroups = size(data,1);

X = [];
groupnames = {};

for i=1:ngroups
    groupsize = size(data{i},1);
    X(end+1:end+groupsize,:) = data{i};
    groupnames(end+1:end+groupsize) = labels(i);
end

% PCA
Xcentered = X-mean(X,1);
[~,Xcentered_S,Xcentered_V] = svd(Xcentered,"econ");
    % For reproducibility, make the sign of the element with the smallest absolute value in each
    % column positive.
    % (Signs of singular vectors that 'svd' function of Matlab outputs can be different 
    % across different machines and releases of Matlab.)
    [~,rowind_eachcol] = min(abs(Xcentered_V),[],1);
    Xcentered_V = Xcentered_V .* sign(Xcentered_V(sub2ind(size(Xcentered_V),rowind_eachcol,1:size(Xcentered_V,2))));
ramanscore = Xcentered*Xcentered_V;
ramanexplained = diag(Xcentered_S).^2/sum(diag(Xcentered_S).^2)*100;

hoge = find(cumsum(ramanexplained)>percent);
nPCs = hoge(1);

Mdl = fitcdiscr(ramanscore(:,1:nPCs),groupnames);
[Vmat,eee] = eigs(Mdl.BetweenSigma,Mdl.Sigma,min(nPCs,ngroups-1)); 
    % For reproducibility, make the sign of the element with the largest absolute value in each
    % column positive.
    % (Signs of eigenvectors that 'eigs' function of Matlab outputs can be different 
    % across different machines and releases of Matlab.)
    [~,rowind_eachcol] = max(abs(Vmat),[],1);
    Vmat = Vmat .* sign( Vmat( sub2ind(size(Vmat),rowind_eachcol,1:size(Vmat,2)) ) );

[~,indexeee] = sort(diag(eee));
Vmat = Vmat(:,indexeee(end:-1:1));

V = normalize(Vmat,1,'norm');

Utmp = ramanscore(:,1:nPCs)*V;

U = cell(ngroups,2);
for i=1:ngroups
    index = strcmp(labels{i},groupnames);
    U{i,2} = Utmp(index,:);    
end

U(:,1) = labels;

end

