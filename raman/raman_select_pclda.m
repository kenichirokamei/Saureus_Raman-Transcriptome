% For reproducibility
rng(1)

% Range of Raman spectra (fingerprint region)
low = 2;
high = 480;

ngroups = size(rawraman,1);

prefilter = @(x) savitzkyGolayFilt(x,4,0,5,[],2);
% Downloaded from 
% https://jp.mathworks.com/matlabcentral/fileexchange/30299-savitzky-golay-smooth-differentiation-filters-and-filter-application
% Note that these argument values cause no smoothing. 
postfilter = @(x) (x-mean(x,2))./std(x,0,2);

trainraman = cell(ngroups,2);
trainraman(:,1) = rawraman(:,1);

% Randomly select 5/6th from data set
for i=1:ngroups
    ramangroup = rawraman{i,2};
    groupsize = size(ramangroup,1);
    
    testn=floor(groupsize/6);
    trainn = groupsize-testn;
    
    remn = sort(randsample(groupsize,testn,false));
    seln = setdiff(1:groupsize,remn);
    
    trainraman{i,2} = postfilter(prefilter(rawraman{i,2}(seln,low:high)));
end

U = pclda(trainraman(:,2),trainraman(:,1),98);

meanU = NaN(ngroups,ngroups-1);
for i=1:ngroups
    meanU(i,:) = mean(U{i,2},1);
end

meanUold = cell(ngroups,1);
for testdata=1:ngroups
    supervisors = setdiff(1:ngroups,testdata);
    nsupervisors = length(supervisors);
    
    Uold = pclda(trainraman(supervisors,2),trainraman(supervisors,1),98);
    
    for i=1:nsupervisors
        meanUold{testdata}(i,:) = mean(Uold{i,2}, 1);
    end
end

