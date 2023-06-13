using CSV, DataFrames, Plots


x2 = CSV.read("x2roc.csv", DataFrame)
mi = CSV.read("miroc.csv", DataFrame)
aic = CSV.read("aicroc.csv", DataFrame)
bic = CSV.read("bicroc.csv", DataFrame)

auc = CSV.read("auc.csv", DataFrame)
aucsd = CSV.read("aucsd.csv", DataFrame)
select!(auc, Not(:Column1))
select!(aucsd, Not(:Column1))
siz = [200,500,1000,2000]
default(fontfamily="Helvetica",
titlelocation=:center, 
framestyle=:box,
linewidth=2.5, background_color=:white,
palette = :rainbow, grid=false,
ylims=[0,1], xlabel ="Size", ylabel="Score",title="AU-ROC")
begin  
    plot(siz, auc.X2, ribbon=aucsd.X2, label ="X2")
    plot!(siz,auc.MI ,ribbon=aucsd.MI, label = "MI")
    plot!(siz,auc.AIC ,ribbon=aucsd.AIC, label = "AIC")
    plot!(siz,auc.BIC ,ribbon=aucsd.BIC, label = "BIC")
    plot!(siz,auc.GAIC ,ribbon=aucsd.GIAC, label = "GAIC")
    plot!(siz,auc.XBIC ,ribbon=aucsd.XBIC, label = "XBIC", color=:goldenrod4)
    plot!(siz,auc.DIMT ,ribbon=aucsd.DIMT, label = "D-IMT", color=:deepskyblue4)
    hline!([0.5], linestyle = :dash, lw=1, label=false)
    hline!([0.7], linestyle = :dash, lw=1, label=false)
    hline!([0.8], linestyle = :dash, lw=1, label=false)
end
plot()
select!(x2, Not(:Column1))
select!(mi, Not(:Column1))
select!(aic, Not(:Column1))
select!(bic, Not(:Column1))

x2sd = CSV.read("x2rocsd.csv", DataFrame)
misd = CSV.read("mirocsd.csv", DataFrame)
aicsd = CSV.read("aicrocsd.csv", DataFrame)
bicsd = CSV.read("bicrocsd.csv", DataFrame)

select!(x2sd, Not(:Column1))
select!(misd, Not(:Column1))
select!(aicsd, Not(:Column1))
select!(bicsd, Not(:Column1))

default()

default(fontfamily="Helvetica",
    titlelocation=:center, 
    framestyle=:box,
    linewidth=2.5, background_color=:white,
    palette = :rainbow, grid=false, yaxis = [0,1], xaxis=[0,1],
    xlabel ="1 -  Specificity", ylabel="Sensitivity")
plot()
begin 
    plot(title="Sample Size: 2000")
    plot!(1 .- x2.Spe,x2.Sen, label = "X2", ribbon=(x2sd.Spe, x2sd.Sen), fillalpha=0.2,marker = ([:c :d], 5, 0.8, Plots.stroke(3, :gray)))
    plot!(1 .- mi.Spe,mi.Sen, label = "MI",ribbon=(misd.Spe, misd.Sen), fillalpha=0.2,marker = ([:c :d], 5, 0.8, Plots.stroke(3, :gray)))
    plot!(1 .- aic.Spe,aic.Sen, label = "AIC",ribbon=(aicsd.Spe, aicsd.Sen), fillalpha=0.2,marker = ([:c :d], 5, 0.8, Plots.stroke(3, :gray)))
    plot!(1 .- bic.Spe,bic.Sen, label = "BIC",ribbon=(bicsd.Spe, bicsd.Sen), fillalpha=0.2,marker = ([:c :d], 5, 0.8, Plots.stroke(3, :gray)))
    Plots.abline!(1,0,line=:dash, color="black", alpha=0.2)
end


using Trapz
x2fpr = 1 .- x2.Spe
x2tpr = x2.Sen
trapz(x2fpr,x2tpr)

mifpr = 1 .- mi.Spe
mitpr = mi.Sen
trapz(mifpr,mitpr)

aicfpr = 1 .- aic.Spe
aictpr = aic.Sen
trapz(aicfpr,aictpr)

bicfpr = 1 .- bic.Spe
bictpr = bic.Sen
trapz(bicfpr,bictpr)


x2fpr = 1 .- x2sd.Spe
x2tpr = x2sd.Sen
trapz(x2fpr,x2tpr)

mifpr = 1 .- misd.Spe
mitpr = misd.Sen
trapz(mifpr,mitpr)

aicfpr = 1 .- aicsd.Spe
aictpr = aicsd.Sen
trapz(aicfpr,aictpr)

bicfpr = 1 .- bicsd.Spe
bictpr = bicsd.Sen
trapz(bicfpr,bictpr)