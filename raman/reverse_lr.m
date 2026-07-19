function [estX,press] = reverse_lr(X,Y,metric)

n = size(X,1);
estX = NaN(size(X,1),size(X,2));
for testdata=1:n
    supervisors = setdiff(1:n,testdata);
    beta = lsqminnorm([ones(size(supervisors,2),1) Y(supervisors,:)], X(supervisors,:));        
    estX(testdata,:) = [ones(size(testdata,2),1) Y(testdata,:)]*beta;
end

if nargin==3
    pre = NaN(n,1);
    for testdata=1:n
        switch metric
            case 1
                pre(testdata) = sum((estX(testdata,:)-X(testdata,:)).^2);
            case 2
                pre(testdata) =  1-dot(estX(testdata,:),X(testdata,:))/(norm(estX(testdata,:))*norm(X(testdata,:)));
            case 3
                corrr = 1-corrcoef(estX(testdata,:)',X(testdata,:)');
                pre(testdata) = corrr(1,2);
            case 4
                pre(testdata) = norm(estX(testdata,:)-X(testdata,:),1);
            case 5
                pre(testdata) = median(abs(estX(testdata,:)-X(testdata,:))./(X(testdata,:)+1),2);
        end
    end
    press = sum(pre);
end

end

