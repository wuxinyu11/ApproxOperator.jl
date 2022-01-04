
function efficiency_meshfree()

to = TimerOutput()
x = GaussPoint(0.0,0.0,0.0,1.0)
@timeit to "Point" begin
    @timeit to "construction MFPoint" xₘ = MFPoint(x,3,Val(:∂1))
    print(xₘ)
end
show(to)

end
