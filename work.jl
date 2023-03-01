using RCall, CSV, DataFrames
@rimport base as r
QMat <- CSV.read("q.csv", DataFrame)

r.as_da