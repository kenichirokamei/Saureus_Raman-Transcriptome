%% estimation

estgenesmean = reverse_lr(genesmean,meanU);

figure;
sz = 5;
for testdata=1:ngroups
    subplot(ceil(ngroups/2),2,testdata);
    hold on
    [sortedgenes,sortedgenesindex] = sort(genesmean(testdata,:));
    errorbar(1:length(sortedgenes),sortedgenes+1,...
        genesse(testdata,sortedgenesindex),genesse(testdata,sortedgenesindex),...
        'o','Color',"#777777",'DisplayName',"Measured")
    s = scatter(1:length(sortedgenes),estgenesmean(testdata,sortedgenesindex)+1,sz, ...
        'MarkerFaceColor',"#bfd147",'MarkerEdgeColor',"#bfd147",'DisplayName',"Estimated"); 
    alpha(s,0.5)
    xlabel('Sorted RNA index')
    ylabel('RNA abundance + 1')
    ylim([0.1 100])
    yticks(logspace(-1,+2,4))
    legend('Location','southeast')
    title(replace(genes_strainnames(testdata),"d","Δ"))
    legend('show')
    set(gca, 'YScale', 'log')
    set(gca,'FontSize',12)
    set(gcf,'renderer','painters')
end

figure;
sz = 5;
for testdata=1:ngroups
    subplot(ceil(ngroups/2),2,testdata);
    hold on
    s = scatter(genesmean(testdata,:),estgenesmean(testdata,:),sz, ...
        'MarkerFaceColor',"#bfd147",'MarkerEdgeColor',"#bfd147");
    alpha(s,0.5)
    xlim([1e-4 1e+2])
    ylim([1e-4 1e+2])
    plot(xlim,ylim,'-','Color',"#111111")
    title(replace(genes_strainnames(testdata),"d","Δ"))
    xlabel("Measured RNA abundance")
    ylabel("Estimated RNA abundance")
    xticks(logspace(-4,+2,4))
    yticks(logspace(-4,+2,4))
    set(gca,'XScale','log')
    set(gca,'YScale','log')
    set(gca,'FontSize',12)
    set(gcf,'renderer','painters')
end


%% permutation test

% For reproducibility
rng(1)

npermutations = 1000;

metric = 1;
% 1: PRESS
% 2: cosine distance
% 3: 1 - Peason's correlation coefficient
% 4: L1 norm
% 5: median of relative errors

[npress,rpresses,significance] = reverse_permuted_lr(genesmean,meanU,npermutations,metric);

figure;
hold on
histogram(rpresses,'FaceColor',"#bfd147",'FaceAlpha',1)
xline(npress,'-',"Original PRESS_t = "+sprintf("%.2e",npress),...
    'LabelOrientation','horizontal','FontSize',15,...
    'Color',"#111111",'LineWidth',1)
xlim(round(xlim,-4))
ylim([0 250])
xlabel("PRESS_t")
ylabel("Frequency")
set(gca,'FontSize',15)
set(gcf,'renderer','painters')

significance

