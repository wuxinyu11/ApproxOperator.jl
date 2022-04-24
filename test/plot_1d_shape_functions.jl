using Revise, YAML, ApproxOperator

config = YAML.load_file("fem.yml")

elements,nodes = importmsh("./msh/bar.msh",config)

sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z],n = 2,γ = 1)

sp(elements["Ω"],elements["Γᵍ"])
