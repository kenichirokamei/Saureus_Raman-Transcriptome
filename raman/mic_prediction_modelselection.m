color = ["#3e90ba","#dab24e","#c35866",...
    "#888888","#888888","#888888","#888888","#888888","#888888",...
    "#3e90ba","#3e90ba",...
    "#CCCCCC","#CCCCCC","#CCCCCC"];
colorf = ["none","none","none",...
    "#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC",...
    "none","none",...
    "#EEEEEE","#EEEEEE","#EEEEEE"];
marker = ["*","x","+","d","^","v",">","<","s","o","p","o","p","h"];


%% AICc

AICc = [];
j = 1;
while ngroups-j-3>0 % The denominator of the penalty term is positive.
    X = lsqminnorm([ones(ngroups,1) meanU(:,1:j)],MIC_VCM_14strains48h);
    RSS = sum((MIC_VCM_14strains48h - [ones(ngroups,1) meanU(:,1:j)]*X).^2);
    AICc(end+1) = string(ngroups*log(2*pi*RSS/ngroups) + ngroups + 2*ngroups*(j+2)/(ngroups-j-3));
    j = j + 1;
end
[~,AICc_min_pos] = min(AICc);

figure;
hold on
plot(1:length(AICc),AICc,'-o',"MarkerSize",12,"LineWidth",1,"Color","#bfd147")
xticks(1:length(AICc))
xlim([1-0.5 length(AICc)+0.5])
xlabel("LDA Axes")
xticklabels(repmat("LDA",1,length(AICc))+"1"+["",repmat(" to LDA",1,length(AICc)-1)]+["",string(2:length(AICc))])
ylabel("AICc")
set(gca,'FontSize',15)

X = lsqminnorm([ones(ngroups,1) meanU(:,1:AICc_min_pos)],MIC_VCM_14strains48h);
axrange = [floor(min([MIC_VCM_14strains48h; [ones(ngroups,1) meanU(:,1:AICc_min_pos)]*X])), ...
         ceil(max([MIC_VCM_14strains48h; [ones(ngroups,1) meanU(:,1:AICc_min_pos)]*X]))];
figure;
    hold on
    for j=1:ngroups
        plot([1 meanU(j,1:AICc_min_pos)]*X,MIC_VCM_14strains48h(j),...
            'Marker',marker(mod(j-1,length(marker))+1),...
            'MarkerEdgeColor',color(mod(j-1,length(color))+1),...
            'MarkerFaceColor',colorf(mod(j-1,length(colorf))+1),...
            'LineStyle','none',...
            'LineWidth',1,...
            'MarkerSize',12,...
            'DisplayName',replace(MIC_VCM_14strains48h_strainnames(j),"d","Δ"))
    end
    plot(axrange,axrange,'-','Color',"#AAAAAA",'DisplayName',"y=x")
    pcc = corrcoef(MIC_VCM_14strains48h,[ones(ngroups,1) meanU(:,1:AICc_min_pos)]*X);
    pcc = pcc(1,2);
    Rsquared = sum(([ones(ngroups,1) meanU(:,1:AICc_min_pos)]*X - mean(MIC_VCM_14strains48h)).^2) / ...
                sum((MIC_VCM_14strains48h - mean(MIC_VCM_14strains48h)).^2);
    text(0.05,0.95,["Multiple correlation coefficient: "+pcc,"Coefficient of determination: "+Rsquared],'FontSize',12,'Units','normalized')
    legend('show','Location','southeast')
    xlim(axrange)
    ylim(axrange)
    ylabel("Vancomycin MIC (μg/mL)")
    xlabel(join(compose("%+.2f",X)+["";repmat("[LDA",AICc_min_pos,1)]+["";split(string(1:AICc_min_pos))'+"]"]))
    title("Minimum MICc")
    set(gca,'FontSize',15)
    set(gcf,'renderer','painters')

    
%% LOOCV

LOOCVMSE = [];
j = 1;
while j<ngroups-1 % The system is not underdetermined. 
    LOOCV_tmp = 0;
    for k=1:ngroups
        teststrain = k;
        trainstrain = setdiff(1:ngroups,teststrain);
        X = lsqminnorm([ones(ngroups-1,1) meanU(trainstrain,1:j)],MIC_VCM_14strains48h(trainstrain));
        LOOCV_tmp = LOOCV_tmp + (MIC_VCM_14strains48h(teststrain) - [1 meanU(teststrain,1:j)]*X).^2;
    end
    LOOCVMSE(end+1) = string(LOOCV_tmp/ngroups);
    j = j + 1;
end
[~,LOOCVMSE_min_pos] = min(LOOCVMSE);

figure;
plot(1:length(LOOCVMSE),LOOCVMSE,'-o',"MarkerSize",12,"LineWidth",1,"Color","#bfd147")
xticks(1:length(LOOCVMSE))
xlim([1-0.5 length(LOOCVMSE)+0.5])
ylim([10^0 10^2*2])
xlabel("LDA axes")
xticklabels(repmat("LDA",1,length(LOOCVMSE))+"1"+["",repmat(" to LDA",1,length(LOOCVMSE)-1)]+["",string(2:length(LOOCVMSE))])
ylabel("MSE of LOOCV")
set(gca,'YScale','log')
set(gca,'FontSize',15)

X = lsqminnorm([ones(ngroups,1) meanU(:,1:LOOCVMSE_min_pos)],MIC_VCM_14strains48h);
axrange = [floor(min([MIC_VCM_14strains48h; [ones(ngroups,1) meanU(:,1:LOOCVMSE_min_pos)]*X])), ...
         ceil(max([MIC_VCM_14strains48h; [ones(ngroups,1) meanU(:,1:LOOCVMSE_min_pos)]*X]))];
figure;
    hold on
    for j=1:ngroups
        plot([1 meanU(j,1:LOOCVMSE_min_pos)]*X,MIC_VCM_14strains48h(j),...
            'Marker',marker(mod(j-1,length(marker))+1),...
            'MarkerEdgeColor',color(mod(j-1,length(color))+1),...
            'MarkerFaceColor',colorf(mod(j-1,length(colorf))+1),...
            'LineStyle','none',...
            'LineWidth',1,...
            'MarkerSize',12,...
            'DisplayName',replace(MIC_VCM_14strains48h_strainnames(j),"d","Δ"))
    end
    plot(axrange,axrange,'-','Color',"#AAAAAA",'DisplayName',"y=x")
    pcc = corrcoef(MIC_VCM_14strains48h,[ones(ngroups,1) meanU(:,1:LOOCVMSE_min_pos)]*X);
    pcc = pcc(1,2);
    Rsquared = sum(([ones(ngroups,1) meanU(:,1:LOOCVMSE_min_pos)]*X - mean(MIC_VCM_14strains48h)).^2) / ...
                sum((MIC_VCM_14strains48h - mean(MIC_VCM_14strains48h)).^2);
    text(0.05,0.95,["Multiple correlation coefficient: "+pcc,"Coefficient of determination: "+Rsquared],'FontSize',12,'Units','normalized')
    legend('show','Location','southeast')
    xlim(axrange)
    ylim(axrange)
    ylabel("Vancomycin MIC (μg/mL)")
    xlabel(join(compose("%+.2f",X)+["";repmat("[LDA",LOOCVMSE_min_pos,1)]+["";split(string(1:LOOCVMSE_min_pos))'+"]"]))
    title("Minimum MSE of LOOCV")
    set(gca,'FontSize',15)
    set(gcf,'renderer','painters')

