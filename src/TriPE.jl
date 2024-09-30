module TriPE

import Base: +, -, *, /, ==, isapprox
import SpecialFunctions: loggamma

export TriParam, TriFun

# Parameter-containing struct
struct TriParam{T<:Real}
	# Base parameters
	α::T
	β::T
	γ::T
	# Shifts
	sα::Int8
	sβ::Int8
	sγ::Int8
	function TriParam(α::T,β::T,γ::T,sα::Int8,sβ::Int8,sγ::Int8) where {T<:Real}
		@assert α > -1 && β > -1 && γ > -1
		@assert sα ≥ 0 && sβ ≥ 0 && sγ ≥ 0
		new{T}(α,β,γ,sα,sβ,sγ)
	end
end

function TriParam(α::T,β::T,γ::T,sα::Int64,sβ::Int64,sγ::Int64) where {T<:Real}
	TriParam(α,β,γ,Int8(sα),Int8(sβ),Int8(sγ))
end

# Dimension of space of bivariate polynomials up to a certain degree
polyDim(deg::Int64) = (deg+1)*(deg+2)÷2

# Polynomial-representing struct
struct TriFun{T}
	κ::TriParam{T}
	degree::Int64
	coefficients::Vector{T}
	function TriFun(κ::TriParam{T},degree::Int64,coefficients::Vector{T}) where T
		@assert polyDim(degree) == length(coefficients)
		new{T}(κ,degree,coefficients)
	end
end

# Unary operations
-(f::TriFun) = TriFun(f.κ, f.degree, -f.coefficients)

# Binary operations
for op = (:+, :-)
	@eval begin
		function ($op)(f::TriFun{T}, g::TriFun{T}) where T
			@assert f.κ == g.κ
			fl = length(f.coefficients)
			gl = length(g.coefficients)
			retl = max(fl, gl)
			retcoefficients = zeros(T, retl)
			retcoefficients[1:fl] = f.coefficients;
			retcoefficients[1:gl] = ($op)(retcoefficients[1:gl], g.coefficients);
			retd = max(f.degree, g.degree)
			TriFun(f.κ, retd, retcoefficients)
		end
	end
end

# Operations with scalars
for op = (:+, :-)
	@eval begin
		function ($op)(f::TriFun{T}, a::Number) where T
			($op)(f, TriFun(f.κ, 0, [convert(T, a)]))
		end
	end
end
for op = (:*, :/)
	@eval begin
		function ($op)(f::TriFun{T}, a::Number) where T
			TriFun(f.κ, f.degree, ($op)(f.coefficients, a))
		end
	end
end
for op = (:+, :*)
	@eval begin
		($op)(a::Number, f::TriFun{T}) where T = ($op)(f, a)
	end
end
-(a::Number, f::TriFun{T}) where T = a + (-f)

# Position range of coefficients of given degree
positionRange(deg::Integer) = (polyDim(deg-1)+1):polyDim(deg)

# Weighted inner product
function h00(κ::TriParam{T}) where T
	α = κ.α+κ.sα; β = κ.β+κ.sβ; γ = κ.γ+κ.sγ;
	exp(loggamma(α+1) + loggamma(β+1) + loggamma(γ+1) - loggamma(α+β+γ+3))
end

end # module
