using CSV, DataFrames
Q = CSV.read("q.csv", DataFrame)
i,j = size(Q)

adjMatrix =  [hcat(zeros(20,20), Matrix(Q));
hcat(transpose(Matrix(Q)), zeros(5,5))]

CSV.write("adj.csv",DataFrame(adjMatrix,:auto))
head = "a" .* string.(1:25)

halfadj = hcat(transpose(Matrix(Q)), zeros(5,5))
CSV.write("halfadj.csv", DataFrame(halfadj, :auto))


newadj = [hcat(zeros(5,5), transpose(Matrix(Q)));
            zeros(20,25)]

CSV.write("newadj.csv", DataFrame(newadj, :auto))

fd = CSV.read("fulldata.csv", DataFrame)
select!(fd,Not(:Column1))