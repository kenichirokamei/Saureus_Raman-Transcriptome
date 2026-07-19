color = ["#3e90ba","#dab24e","#c35866",...
    "#888888","#888888","#888888","#888888","#888888","#888888",...
    "#3e90ba","#3e90ba",...
    "#CCCCCC","#CCCCCC","#CCCCCC"];
colorf = ["none","none","none",...
    "#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC",...
    "none","none",...
    "#EEEEEE","#EEEEEE","#EEEEEE"];
marker = ["*","x","+","d","^","v",">","<","s","o","p","o","p","h"];

pltind = [1:3]; % N315, Mu3, Mu50

% For reproducibility
rng(1)
nsamplings = 100;
rat = 9/10;
nthpc = 2;

pcaalltopexpl = NaN(nsamplings,nthpc);

pcaallstrainlatent = cell(length(pltind),2);
    pcaallstrainlatent(:,1) = rawraman(pltind,1);
    for i=1:length(pltind)
        pcaallstrainlatent{i,2} = NaN(nsamplings,nthpc);
    end

evalrange = low:high;

normalized_rawraman = cell(ngroups,2);
normalized_rawraman(:,1) = rawraman(:,1);
for i=1:ngroups
    normalized_rawraman{i,2} = postfilter(rawraman{i,2}(:,evalrange));
end

for j=1:nsamplings
    normalized_rawraman_selected = cell(ngroups,2);
    normalized_rawraman_selected(:,1) = rawraman(:,1);
    for i=1:ngroups
        ncells = size(rawraman{i,2},1);
        selected_cells = randsample(ncells,round(ncells*rat));
        normalized_rawraman_selected{i,2} = normalized_rawraman{i,2}(selected_cells,:);
    end
    
    normalized_rawraman_selected_cat = [];
    for i=pltind
        normalized_rawraman_selected_cat = [normalized_rawraman_selected_cat; normalized_rawraman_selected{i,2}];
    end
    [~,score,~,~,explained] = pca(normalized_rawraman_selected_cat);
    pcaalltopexpl(j,:) = explained(1:nthpc);
    cellind = 0;
    straincounter = 0;
    for i=pltind
        cellind = cellind(end)+1:cellind(end)+size(normalized_rawraman_selected{i,2},1);
        straincounter = straincounter + 1;
        for k=1:nthpc
            pcaallstrainlatent{straincounter,2}(j,k) = var(score(cellind,k));
        end
    end
end

for k=1:nthpc

    figure;
    hold on
    counter = zeros(1,2);
    for i=1:nsamplings
        tmp = NaN(length(pltind),1);
        for j=1:length(pltind)
            plot(pcaallstrainlatent{j,2}(i,k),pltind(j), ...
                    'MarkerEdgeColor',color(mod(pltind(j)-1,length(color))+1), ...
                    'MarkerFaceColor',colorf(mod(pltind(j)-1,length(colorf))+1), ...
                    'Marker',marker(mod(pltind(j)-1,length(marker))+1), ...
                    'LineStyle','none','LineWidth',1, ...
                    'MarkerSize',15)
            tmp(j) = pcaallstrainlatent{j,2}(i,k);
        end
        if isequal(pltind,1:3)
            if tmp(1)<tmp(2)
                counter(:,1) = counter(:,1) + 1;
            end
            if tmp(3)<tmp(2)
                counter(:,2) = counter(:,2) + 1;
            end
        end
        plot(tmp,pltind,'-','Color',"#AAAAAA")
    end
    if isequal(pltind,1:3)
        xaxrange = xlim;
        text(diff(xaxrange)*0.04+xaxrange(1),1.5,[counter(1)/nsamplings*100+"% of resamples", "show Mu3 > N315"],'FontSize',12)
        text(diff(xaxrange)*0.04+xaxrange(1),2.5,[counter(2)/nsamplings*100+"% of resamples", "show Mu3 > Mu50"],'FontSize',12)
    end
    yticks(pltind)
    yticklabels(replace(string(extractBetween(rawraman(pltind,1),"","_repall")),"d","Δ"))
    ylim([pltind(1)-0.5 pltind(end)+0.5])
    xlabel(["Within-strain variance along PC"+k, ...
        "(Mean variance explained by PC"+k+" axes: "+mean(pcaalltopexpl(:,k))+"%)"])
    set(gca,'FontSize',15)
    set(gca,'YDir','reverse')
    set(gcf,'renderer','painters')

end

figure;
margin = 30;
for i=1:length(pltind)
    subplot(length(pltind),4,(i-1)*4+1)
    hold on
    yline(0,'LineStyle',':')
    s = shadedErrorBar(evalrange,normalized_rawraman{pltind(i),2},{@mean,@std});
    set(s.patch,'FaceColor',color(mod(pltind(i)-1,length(color))+1))
    set(s.edge,'LineStyle','none')
    set(s.mainLine,'LineWidth',1,'Color',color(mod(pltind(i)-1,length(color))+1))
    text(0.08,0.9,replace(string(extractBetween(normalized_rawraman(pltind(i),1),"","_repall")),"d","Δ"), ...
        'Units','Normalized','FontSize',10)
    if pltind(i)==pltind(end)
        ylabel("Normalized intensity")
    end
    xticks(rawramanindex2Ramanshift{1,2})
    if pltind(i)==pltind(end)
        xticklabels(rawramanindex2Ramanshift{2,2})
    else
        xticklabels({})
    end
    xlim([evalrange(1)-margin evalrange(end)+margin])
    set(gca,'FontSize',15)
end
subplot(length(pltind),4,setdiff(1:length(pltind)*4,([1:length(pltind)]-1)*4+1))
hold on
for i=1:length(pltind)
    s = shadedErrorBar(evalrange,normalized_rawraman{pltind(i),2},{@mean,@std});
    set(s.patch,'FaceColor',color(mod(pltind(i)-1,length(color))+1))
    set(s.edge,'LineStyle','none')
    set(s.mainLine,'LineWidth',1,'Color',color(mod(pltind(i)-1,length(color))+1))
end
yline(0,'LineStyle','--')
legend([replace(string(extractBetween(normalized_rawraman(pltind,1),"","_repall")),"d","Δ"); "y=0"])
% ylabel("Normalized intensity")
xlabel("Raman shift (cm^{-1})",'Interpreter','tex')
xticks(rawramanindex2Ramanshift{1,2})
xticklabels(rawramanindex2Ramanshift{2,2})
xlim([evalrange(1)-margin evalrange(end)+margin])
set(gca,'FontSize',15)
set(gcf,'renderer','painters')

