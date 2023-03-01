using CSV, DataFrames, Plots


x2 = CSV.read("x2roc.csv", DataFrame)
mi = CSV.read("miroc.csv", DataFrame)
aic = CSV.read("aicroc.csv", DataFrame)
bic = CSV.read("bicroc.csv", DataFrame)

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
    palette = :rainbow, grid=false, xlims =[0,1], ylims = [0,1],
    xlabel ="1 - Specificity", ylabel="Sensitivity")

begin 
    plot(title="Sample Size: 1000")
    plot!(1 .- x2.Spe,x2.Sen, label = "X2", ribbon=x2sd.Sen)
    plot!(1 .- mi.Spe,mi.Sen, label = "MI")
    plot!(1 .- aic.Spe,aic.Sen, label = "AIC")
    plot!(1 .- bic.Spe,bic.Sen, label = "BIC")
    Plots.abline!(1,-0.0002,line=:dash, color="black", alpha=0.2)
end



