# This file is a part of JuliaFEM/MiniBall.jl.
# License is GPL: see https://github.com/JuliaFEM/MiniBall.jl/blob/master/LICENSE.md

using Miniball
using Base.Test
using TimerOutputs

include("test_unit_circle.jl")
exit(0)

pkg_dir = Pkg.dir("Miniball")
maybe_test_files = readdir(joinpath(pkg_dir, "test"))
is_test_file(fn) = startswith(fn, "test_") & endswith(fn, ".jl")
test_files = filter(is_test_file, maybe_test_files)
info("$(length(test_files)) test files found.")

const to = TimerOutput()
@testset "Miniball.jl" begin
    for fn in test_files
        info("----- Running tests from file $fn -----")
        t0 = time()
        timeit(to, fn) do
            include(fn)
        end
        dt = round(time() - t0, 2)
        info("----- Testing file $fn completed in $dt seconds -----")
    end
end

println()
println("Test statistics:")
println(to)

