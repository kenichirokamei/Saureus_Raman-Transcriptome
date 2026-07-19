% For reproducibility
rng(2)

color = ["#3e90ba","#dab24e","#c35866",...
    "#888888","#888888","#888888","#888888","#888888","#888888",...
    "#3e90ba","#3e90ba",...
    "#CCCCCC","#CCCCCC","#CCCCCC"];
colorf = ["none","none","none",...
    "#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC","#CCCCCC",...
    "none","none",...
    "#EEEEEE","#EEEEEE","#EEEEEE"];
marker = ["*","x","+","d","^","v",">","<","s","o","p","o","p","h"];

U_fortsne = [];
U_names_fortsne = [];
for i=1:size(U,1)
    U_fortsne = [U_fortsne; U{i,2}];
    U_names_fortsne = [U_names_fortsne; repmat(string(U{i,1}), size(U{i,2},1), 1)];
end

tsneresult = tsne(U_fortsne);

figure;
g = gscatter(tsneresult(:,1),tsneresult(:,2),...
        replace(extractAfter(extractBefore(U_names_fortsne,"_repall"),""),"d","Δ"),...
        hexadecimalcolorcode2rgbtriplet(color),char(join(marker,'')));
for i=1:length(g)
    set(g(i),'MarkerFaceColor',colorf(i))
end

title("Visualization of LDA space by t-SNE")
xlabel("t-SNE1")
ylabel("t-SNE2")
set(gca,'FontSize',15)
set(gcf,'renderer','painters')

