
using BenchmarkTools
using TimerOutputs

struct T{S}
    type::Val{S}
end
struct R end
(::T{1})() = 1.0
(::T{2})() = 1.0
(::T{3})() = 1.0
f(::Val{1}) = 1.0
# f(::Val{2}) = 1.0
# f(::Val{3}) = 1.0
g(::R) = 1.0
function test(n::Int)
    to = TimerOutput()
    t = T(Val(1))
    r = R()

    a = 0.0
    val1 = Val(1)
    @timeit to "type" begin
        for i in 1:n
            t()
        end
    end
    @timeit to "struct" begin
        for i in 1:n
            g(r)
        end
    end
    @timeit to "val" begin
        for i in 1:n
            f(val1)
        end
    end
    show(to)
end
