using CSV, DataFrames, Plots


x2 = CSV.read("x2.csv", DataFrame)


mi = CSV.read("mi.csv", DataFrame)
aic = CSV.read("aic.csv", DataFrame)
bic = CSV.read("bic.csv", DataFrame)

select!(x2, Not(:Column1))
select!(mi, Not(:Column1))
select!(aic, Not(:Column1))
select!(bic, Not(:Column1))
siz = [500, 1000, 2000]

x2sd = CSV.read("x2sd.csv", DataFrame)


misd = CSV.read("misd.csv", DataFrame)
aicsd = CSV.read("aicsd.csv", DataFrame)
bicsd = CSV.read("bicsd.csv", DataFrame)

select!(x2sd, Not(:Column1))
select!(misd, Not(:Column1))
select!(aicsd, Not(:Column1))
select!(bicsd, Not(:Column1))


begin
    default(fontfamily="Helvetica",
        titlelocation=:center, 
        framestyle=:box,
        linewidth=2.5, background_color=:white,
        palette = :rainbow, grid=false,
        ylims=[0,1], xlabel ="Size", ylabel="Score")
    
    AccPlot = plot(siz, x2.Acc, ribbon = x2sd.Acc,label="X2", title="Accuracy")
    plot!(AccPlot, siz, mi.Acc, label="MI", ribbon = misd.Acc)
    plot!(AccPlot, siz, aic.Acc, label="AIC", ribbon = aicsd.Acc)
    plot!(AccPlot, siz, bic.Acc, label="BIC", ribbon = bicsd.Acc)

    SenPlot =  plot(siz, x2.Sen, label="X2", title="Sensitivity",ribbon = x2sd.Sen)
    plot!(SenPlot, siz, mi.Sen, label="MI", ribbon = misd.Sen)
    plot!(SenPlot, siz, aic.Sen, label="AIC", ribbon = aicsd.Sen)
    plot!(SenPlot, siz, bic.Sen, label="BIC", ribbon = bicsd.Sen)

    SpePlot =  plot(siz, x2.Spe, label="X2", title="Specificity",ribbon = x2sd.Sen)
    plot!(SpePlot, siz, mi.Spe, label="MI", ribbon = misd.Spe)
    plot!(SpePlot, siz, aic.Spe, label="AIC", ribbon = aicsd.Spe)
    plot!(SpePlot, siz, bic.Spe, label="BIC", ribbon = bicsd.Spe)

    AUCPlot =  plot(siz, x2.AUC, label="X2", title="AUC",ribbon = x2sd.Sen)
    plot!(AUCPlot, siz, mi.AUC, label="MI", ribbon = misd.AUC)
    plot!(AUCPlot, siz, aic.AUC, label="AIC", ribbon = aicsd.AUC)
    plot!(AUCPlot, siz, bic.AUC, label="BIC", ribbon = bicsd.AUC)

    l = @layout [a ; c d]

    plot(AccPlot, SenPlot, SpePlot, size= (700,800),
    plot_title = "alpha ≤ 0.05, thresh ≤ 1", layout=l)
end

#bic.Spe= parse.(Float64, bic.Spe)
#bic.Spe[4] = "0" 
#aic
#x2
clipboard(sprint(show, "text/tab-separated-values",  DataFrame(x2mean=x2.AUC,SD=x2sd.AUC)))

clipboard(sprint(show, "text/csv", DataFrame(x2mean=x2.AUC,SD=x2sd.AUC)))

clipboard(sprint(show, "text/tab-separated-values",  DataFrame(mean=mi.AUC,SD=misd.AUC)))
clipboard(sprint(show, "text/tab-separated-values",  DataFrame(mean=bic.AUC,SD=bicsd.AUC)))