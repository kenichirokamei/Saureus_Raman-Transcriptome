function [npress,rpresses,significance,rp] = new_permuted_plsr(X,oldY,Y,plsnum,npermutations)

n = size(X,1);

if factorial(n)<npermutations
    npermutations = factorial(n)-1;
    rp = perms(1:n);
    rp(end,:) = [];
else
    rp = zeros(npermutations,n);
    for k=1:npermutations
        rp(k,:) = randperm(n);
    end
end

rpresses = NaN(1,npermutations);
for i=1:npermutations
    [~,rpresses(i)] = new_plsr(X(rp(i,:),:),oldY,Y,plsnum);
end
[~,npress] = new_plsr(X,oldY,Y,plsnum);

% p-value correction according to 
% Phipson, B., and Smyth, G. K. (2010). Permutation p-values should never be zero: calculating exact
% p-values when permutations are randomly drawn. Stat. Appl. Genet. Molec. Biol. Volume 9, Issue
% 1, Article 39.
significance = (sum(rpresses<=npress)+1)/(npermutations+1);

