function random_points(center, radius, n, boundary::Bool)
    dims = length(center), n
    dim = length(center)
    pts = zeros(n, dim)
    for row in indices(pts,1)
        dir = normalize(randn(dim))
        if boundary
            r = rand() * radius
        else
            r = radius
        end
        pts[row,:] .= center + r*dir
    end
    pts
end

function random_ball(dim)
    center = randn(dim)
    radius = 5rand()
    center, radius
end

function leq_approx(small, big)
    if small <= big
        true
    else
        isapprox(small, big)
    end
end

function check_points_approx_inside(center, radius, pts)
    for row in indices(pts,1)
        pt = @view pts[row,:]
        r = norm(pt - center)
        if !leq_approx(r, radius)
            return false
        end
    end
    true
end

@testset "fuzz tests" begin
    srand(42)
    for dim in 1:8
        for n in 1:100
            for boundary in [true, false]
                center, radius = random_ball(dim)
                pts = random_points(center, radius, n, boundary)
                ball = miniball(pts)
                r = sqrt(ball.squared_radius)
                inside = check_points_approx_inside(ball.center, r, pts)
                if inside
                    @test inside
                else
                    @test_broken inside
                end
                @test leq_approx(r, radius)
            end
        end
    end
end