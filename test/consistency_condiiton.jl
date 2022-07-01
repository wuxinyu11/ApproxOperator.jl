
using Revise, ApproxOperator, BenchmarkTools, YAML

config = YAML.load_file("./fem.yml")
elements = importmsh("./msh/patchtest.msh",config)

# @btime set𝝭!(elements["Ω"])
# @btime set∇𝝭!(elements["Ω"])
# set∇𝝭!(elements["Ω"])
# set∇̃𝝭!(elements["Ω"])
# set∇̃²𝝭!(elements["Ω"])
# @btime set∇²𝝭!(elements["Ω"])
# set∇³𝝭!(elements["Ω"])
# set∇̂³𝝭!(elements["Ω"])
# set∇̃𝝭!(elements["Ω"])
# set∇̃²𝝭!(elements["Ω"])
# f = checkConsistency(elements["Ω"])
# f = checkConsistency(elements["Ω"],ApproxOperator.get∇𝝭,ApproxOperator.get∇𝒑)
# f = checkConsistency(elements["Ω"],ApproxOperator.get∇²𝝭,ApproxOperator.get∇²𝒑)
# f = checkConsistency(elements["Ω"],ApproxOperator.get∇³𝝭,ApproxOperator.get∇³𝒑)

# err_chol, err, err_x, err_y = ApproxOperator.test_cal𝗠!(elements["Ω"][1],(0.,0.,0.))
# err,err_chol,err_inv,err_I,err1,err2 = ApproxOperator.test_cal𝗠!(elements["Ω"][1],(0.,0.,0.))