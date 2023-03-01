using CSV, DataFrames, Plots


x2 = CSV.read("x2.csv", DataFrame)


mi = CSV.read("mi.csv", DataFrame)
aic = CSV.read("aic.csv", DataFrame)
bic = CSV.read("bic.csv", DataFrame)

select!(x2, Not(:Column1))
select!(mi, Not(:Column1))
select!(aic, Not(:Column1))
select!(bic, Not(:Column1))
siz = [20,50,100,200,500,1000]

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
    
    AccPlot = plot(siz, x2.Acc, label="X2", title="Accuracy")
    plot!(AccPlot, siz, mi.Acc, label="MI")
    plot!(AccPlot, siz, aic.Acc, label="AIC")
    plot!(AccPlot, siz, bic.Acc, label="BIC")

    SenPlot =  plot(siz, x2.Sen, label="X2", title="Sensitivity")
    plot!(SenPlot, siz, mi.Sen, label="MI")
    plot!(SenPlot, siz, aic.Sen, label="AIC")
    plot!(SenPlot, siz, bic.Sen, label="BIC")

    SpePlot =  plot(siz, x2.Spe, label="X2", title="Specificity")
    plot!(SpePlot, siz, mi.Spe, label="MI")
    plot!(SpePlot, siz, aic.Spe, label="AIC")
    plot!(SpePlot, siz, bic.Spe, label="BIC")

    AUCPlot =  plot(siz, x2.AUC, label="X2", title="AUC")
    plot!(AUCPlot, siz, mi.AUC, label="MI")
    plot!(AUCPlot, siz, aic.AUC, label="AIC")
    plot!(AUCPlot, siz, bic.AUC, label="BIC")

    l = @layout [a b; c d]

    plot(AccPlot, SenPlot, SpePlot, size= (700,800),
    plot_title = "alpha ≤ 0.05, thresh ≤ 1", layout=l)
end

#bic.Spe= parse.(Float64, bic.Spe)
#bic.Spe[4] = "0" 
#aic
#x2