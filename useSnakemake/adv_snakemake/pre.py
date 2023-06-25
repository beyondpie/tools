import snapatac2 as sa2
print(sa2.__version__)

a = snakemake.output
print(type(a))
print(a)
print(snakemake.output[0])

print(type(snakemake.output['a1']))

b = snakemake.params
print(type(b))
print(b)
print(snakemake.params[0])
print(snakemake.params["knn"])

