using CSV, DataFrames, Plots
default(fontfamily="Helvetica",
        titlelocation=:center, 
        framestyle=:box,
        linewidth=1.5, background_color=:white,
        palette = :seaborn_bright, grid=false)
acc_0051 = CSV.read("acc_at0.05nd1.csv", DataFrame)
sep_0051 = CSV.read("spe_at0.05nd1.csv", DataFrame)
sen_0051 = CSV.read("sen_at0.05nd1.csv", DataFrame)
auc_0051 = CSV.read("auc_at0.05nd1.csv", DataFrame)
siz = [1000,500,200,100,50,20]

begin 
    acc_0051 = CSV.read("acc_at0.05nd1.csv", DataFrame)
    sep_0051 = CSV.read("spe_at0.05nd1.csv", DataFrame)
    sen_0051 = CSV.read("sen_at0.05nd1.csv", DataFrame)
    auc_0051 = CSV.read("auc_at0.05nd1.csv", DataFrame)
    p1 = plot(siz,acc_0051.V1, ylims=[0,1], label ="X2",
            xlabel= "Size", ylabel="Score", title="Accuracy");
    plot!(p1,siz,acc_0051.V2, ylims=[0,1], label ="MI");
    plot!(p1,siz,acc_0051.V3, ylims=[0,1], label ="AIC");
    plot!(p1,siz,acc_0051.V4, ylims=[0,1], label ="BIC",color="black");

    p2 = plot(siz,sen_0051.V1, ylims=[0,1], label ="X2",
            xlabel= "Size", ylabel="Score", title="Sensitivity");
    plot!(p2,siz,sen_0051.V2, ylims=[0,1], label ="MI");
    plot!(p2,siz,sen_0051.V3, ylims=[0,1], label ="AIC");
    plot!(p2,siz,sen_0051.V4, ylims=[0,1], label ="BIC", color="black") ;

    p3 = plot(siz,sep_0051.V1, ylims=[0,1], label ="X2",
            xlabel= "Size", ylabel="Score", title="Specificity");
    plot!(p3,siz,sep_0051.V2, ylims=[0,1], label ="MI");
    plot!(p3,siz,sep_0051.V3, ylims=[0,1], label ="AIC");
    plot!(p3,siz,sep_0051.V4, ylims=[0,1], label ="BIC", color="black") ;

    p4 = plot(siz,auc_0051.V1, ylims=[0,1], label ="X2",
            xlabel= "Size", ylabel="Score", title="AUC");
    plot!(p4,siz,auc_0051.V2, ylims=[0,1], label ="MI");
    plot!(p4,siz,auc_0051.V3, ylims=[0,1], label ="AIC");
    plot!(p4,siz,auc_0051.V4, ylims=[0,1], label ="BIC", color="black") ;

    plot(p1,p2,p3,p4, size = (800,800),
        plot_title = "pval < 0.05; threshold < 0.5")
end