B = genesmean'*(pinv([ones(ngroups,1) meanU]))';

xax = 1;
yax = 2;

% Specify the x-axis range where the genes you want to circle in red are located.
x_thr1 = 1.7;
x_thr2 = -2.6;
highlight_ind = find(B(:,xax+1)./B(:,1)>x_thr1 | B(:,xax+1)./B(:,1)<x_thr2);

figure;
hold on
scatter(B(:,xax+1)./B(:,1),B(:,yax+1)./B(:,1),[],[0.5 0.5 0.5],'.','SizeData',30)
for i=1:length(highlight_ind)
    plot(B(highlight_ind(i),xax+1)./B(highlight_ind(i),1),B(highlight_ind(i),yax+1)./B(highlight_ind(i),1),'o','Color','r')
    text(B(highlight_ind(i),xax+1)./B(highlight_ind(i),1),B(highlight_ind(i),yax+1)./B(highlight_ind(i),1),"   "+genes_description(1,highlight_ind(i)),'FontSize',12)
end
xline(0,'--','Color','#AAAAAA');
yline(0,'--','Color','#AAAAAA');
xlabel("$\beta^{\mathrm{LDA}1}$",'Interpreter','latex')
ylabel("$\beta^{\mathrm{LDA}2}$",'Interpreter','latex')
set(gca,'FontSize',15)
set(gcf,'renderer','painters')

% To plot the expression patterns of the genes circled in red, uncomment this section.
%{
for i=1:length(highlight_ind)
    figure('Position', [10 10 250 250]);
    errorbar(1:ngroups,genesmean(:,highlight_ind(i)),genesse(:,highlight_ind(i)),'-o')
    xline(3.5,'--')
    title([genes_description(1,highlight_ind(i)),"\it "+genes_description(2,highlight_ind(i))])
    xticks(1:ngroups)
    xlim([0.5 ngroups+0.5])
    xticklabels(replace(genes_strainnames,"d","Δ"))
    xtickangle(45)
    ylabel("Abundance (intensity)")
    set(gca,'FontSize',15)
    set(gcf,'renderer','painters')
end
%}

