function [newestY,press] = new_plsr(X,oldY,Y,plsnum)

n = size(X,1);
newestY = NaN(size(Y));
for testdata=1:n
    supervisors = setdiff(1:n,testdata);
    
    [~,~,~,~,beta] = plsregress(X(supervisors,:),oldY{testdata},plsnum);
    estY = [ones(size(X(testdata,:),1),1) X(testdata,:)]*beta;
    
    % change bases
    newestY(testdata,:) = estY*pinv(oldY{testdata})*Y(supervisors,:);
end

press = sum((newestY-Y).^2,'all');

end

