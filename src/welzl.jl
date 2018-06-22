export MiniballAlgorithm
export Welzl
export Gaertner
using StaticArrays

abstract type MiniballAlgorithm end
struct Welzl <: MiniballAlgorithm

end
struct Gaertner <: MiniballAlgorithm end

function miniball(pts::AbstractMatrix, ::Gaertner)
    miniball(pts)
end

function miniball(pts::AbstractVector, ::Welzl)
    welzl(pts)
end

function miniball(m::AbstractMatrix, w::Welzl)
    pts = vector_of_static_points(m)
    miniball(pts, w)
end

function vector_of_static_points(m::AbstractMatrix)
    V = SVector{size(m,2), eltype(m)}
    ret = V[]
    for i in indices(m,1)
        v = V(@view m[i,:])
        push!(ret,v)
    end
    ret
end

@noinline function welzl(pts::AbstractVector)
    bdry = eltype(pts)[]
    pts = shuffle!(copy(pts))
    welzl(pts, bdry)
end

function welzl(pts, bdry)
    dim = length(eltype(pts))
    if length(pts) == 1 && isempty(bdry)
        return welzl(bdry,pts)
    end
    if isempty(pts) || length(bdry) > dim
        return sphere_bdry(bdry...)
    end
    pts0 = copy(pts)
    pt = pop!(pts0)
    center0, r0 = welzl(pts0, bdry)
    if norm(pt - center0) <= r0 + sqrt(eps(r0))
        center0, r0
    else
        bdry1 = copy(bdry)
        push!(bdry1, pt)
        welzl(pts0, bdry1)
    end
end

const _2D = Length{2}
const _3D = Length{3}
const _ND = Length

function sphere_bdry(::_ND, pt)
    pt, zero(eltype(pt))
end
function sphere_bdry(::_ND, pt1, pt2)
    center = map(middle,pt1, pt2)
    r = norm(pt1 - pt2) /2
    center, r
end
function sphere_bdry_2d_origin(b,c)
    bx, by = b
    cx, cy = c
    B=bx*bx+by*by
    C=cx*cx+cy*cy
    D=bx*cy-by*cx
    center = (cy*B-by*C)/(2*D)
    r = (bx*C-cx*B)/(2*D)
    center, abs(r)
end
function sphere_bdry(::_2D, pt1, pt2, pt3)
    center, r = sphere_bdry_2d_origin(pt2 - pt1, pt3 - pt1)
    center + pt1, r
end

function sphere_bdry(pt::AbstractVector, pts...)
    sphere_bdry(Length(pt), pt, pts...)
end
