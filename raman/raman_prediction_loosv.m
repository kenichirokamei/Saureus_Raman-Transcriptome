%% estimation

plsnum = ngroups-3;
estmeanU = new_plsr(genesmean,meanUold,meanU,plsnum);

% dimensions to plot
dims = [1,2];
if size(dims,2)~=2
    error('Specify two dimensions using a row vector.')
end

color = ["#3e90ba","#dab24e","#c35866",...
    "#888888","#888888","#888888","#888888","#888888","#888888",...
    "#3e90ba","#3e90ba",...
    "#CCCCCC","#CCCCCC","#CCCCCC"];
colorf = ["none","none","none",...
    "#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC",...
    "none","none",...
    "#EEEEEE","#EEEEEE","#EEEEEE"];
marker = ["*","x","+","d","^","v",">","<","s","o","p","o","p","h"];

figure;
hold on
for i=1:ngroups
    plot(U{i,2}(:,dims(1)),U{i,2}(:,dims(2)),...
        'Marker',marker(mod(i-1,length(marker))+1),...
        'Color',"#CCCCCC",...
        'LineStyle','none',...
        'LineWidth',0.6)
end
set(gca,'DefaultLineLineWidth',0.6)
set(gca,'LineStyleOrder',{'-'},'ColorOrder',hexadecimalcolorcode2rgbtriplet(color),'ColorOrderIndex',1);
for i=1:ngroups
    cova = cov(U{i,2}(:,dims));
    me = mean(U{i,2}(:,dims));
    e = error_ellipse(cova,me,'conf',0.95);
    % Downloaded from https://jp.mathworks.com/matlabcentral/fileexchange/4705-error_ellipse
end
p = [];
for i=1:ngroups
    p(end+1) = plot(estmeanU(i,dims(1)),estmeanU(i,dims(2)), ...
        'Marker',marker(mod(i-1,length(marker))+1),...
        'MarkerEdgeColor',color(mod(i-1,length(color))+1),...
        'MarkerFaceColor',colorf(mod(i-1,length(colorf))+1),...
        'LineStyle','none',...
        'LineWidth',1.5,...
        'MarkerSize',8,...
        'DisplayName',replace(extractBetween(string(U(i,1)),"","_repall"),"d","Δ"));
end
title("Measured (light) + Estimated from transcriptomes (bold)")
xlabel("LDA"+dims(1))
ylabel("LDA"+dims(2))
legend(p,'Interpreter','none')
set(gca,'FontSize',15)
set(gcf,'renderer','painters')


%% permutation test

% For reproducibility
rng(1)

npermutations = 1000;
plsnum = ngroups-3;
[npress,rpress,significance] = new_permuted_plsr(genesmean,meanUold,meanU,plsnum,npermutations);

figure;
hold on
histogram(rpress,'FaceColor',"#bfd147",'FaceAlpha',1)
xline(npress,'-',"Original PRESS_{r'} = "+sprintf("%.2e",npress),...
    'LabelOrientation','horizontal','FontSize',15,...
    'Color',"#111111",'LineWidth',1)
xlabel("PRESS_{r'}")
ylabel("Frequency")
ylim([0 150])
set(gca,'FontSize',15)
set(gcf,'renderer','painters')

significance

