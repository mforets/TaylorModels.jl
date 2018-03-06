# rpa_functions.jl

"""
   _rpaar(f::Function, x0::Interval, ii::Interval, _order::Int)

Rigurous polynomial approximation (RPA) with absolute remainder
for the function `f` on the interval `ii`,  using a Taylor expansion
around the *interval* `x0` of order `_order`. The bound is computed
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

"""
function _rpaar(f::Function, x0::Interval{T}, ii::Interval{T}, _order::Int) where {T}
    polf  = f( x0+Taylor1(Interval{T}, _order) )
    polfI = f( ii+Taylor1(Interval{T}, _order+1) )
    Δ = boundarem(f, polf, polfI, x0, ii)
    return TMAbsRem(polf, Δ, x0, ii)
end


"""
   _rparr(f::Function, x0::Interval, ii::Interval, _order::Int)

Rigurous polynomial approximation (RPA) with relative remainder
for the function `f` on the interval `ii`,  using a Taylor expansion
around the *interval* `x0` of order `_order`. The bound is computed
exploiting monotonicity if possible, otherwise, it uses the Lagrange
coefficient.

"""
function _rparr(f::Function, x0::Interval{T}, ii::Interval{T}, _order::Int) where {T}
    polf  = f( x0+Taylor1(Interval{T}, _order) )
    polfI = f( ii+Taylor1(Interval{T}, _order+2) )
    Δ = boundrrem(f, polf, polfI, x0, ii)
    return TMRelRem(polf, Δ, x0, ii)
end


"""
   rpa(g::Function, tmf::TMAbsRem)

Rigurous polynomial approximation (RPA) for the function `g` using the
Taylor Model with absolute remainder `tmf`. The bound is computed
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

"""
function rpa(g::Function, tmf::TMAbsRem)
    _order = get_order(tmf)

    # Short-cut if `tmf` is the independent variable
    tmf == TMAbsRem(_order, tmf.x0, tmf.iI) &&
        return _rpaar(g, tmf.x0, tmf.iI, _order)

    f_pol = tmf.pol
    Δf = remainder(tmf)
    x0 = tmf.x0
    iI = tmf.iI

    # Range of tmf including remainder (Δf)
    range_tmf = bound_taylor1(f_pol, iI-x0) + Δf

    # Compute RPA for `g`, around f_pol[0], over range_tmf
    tmg = _rpaar(g, f_pol[0], range_tmf, _order)

    # Use original independent variable
    tm1 = tmf - f_pol[0]
    tmres = tmg( tm1 )

    # Final remainder
    Δ = remainder(tmres) + remainder(tmg)
    return TMAbsRem(tmres.pol, Δ, x0, iI)
end


"""
   rpa(g::Function, tmf::TMRelRem)

Rigurous polynomial approximation (RPA) for the function `g` using the
Taylor Model with absolute remainder `tmf`. The bound is computed
exploiting monotonicity if possible, otherwise, it uses Lagrange bound.

"""
function rpa(g::Function, tmf::TMRelRem)
    _order = get_order(tmf)

    # Do not overestimate if `tmf` is the independent variable
    tmf == TMRelRem(_order, tmf.x0, tmf.iI) &&
        return _rparr(g, tmf.x0, tmf.iI, _order)

    f_pol = tmf.pol
    Δf = remainder(tmf)
    x0 = tmf.x0
    iI = tmf.iI

    # Range of tmf including remainder (Δf)
    range_tmf = bound_taylor1(f_pol, iI-x0) + Δf * (iI-x0)^(_order+1)

    # Compute RPA for `g`, around f_pol[0], over range_tmf
    tmg = _rparr(g, f_pol[0], range_tmf, _order)
    tm1 = tmf - f_pol[0]
    tmres = tmg( tm1 )

    tmn = TMRelRem(Taylor1(copy(tm1.pol.coeffs)), tm1.rem, tm1.x0, tm1.iI)
    for i = 1:_order
        tmn = tmn * tm1
    end
    Δ = remainder(tmres) + remainder(tmn) * remainder(tmg)
    return TMRelRem(tmres.pol, Δ, x0, iI)
end


# evaluate, and function-like evaluation for TMAbsRem
for TM in tupleTMs
    @eval function evaluate(tmg::$TM, tmf::$TM)
        _order = get_order(tmf)
        @assert _order == get_order(tmg)

        tmres = $TM(zero(tmg.pol[0]), _order, tmf.x0, tmf.iI)
        @inbounds for k = _order:-1:0
            tmres = tmres * tmf
            tmres = tmres + tmg.pol[k]
        end
        return tmres
    end

    @eval (tm::$TM)(x::$TM) = evaluate(tm, x)
end


"""
    rpafp(tm::TMAbsRem{T})
    rpafp(tm::TMRelRem{T})

Convert a `tm` TaylorModel to a T-type RPA. It returns the `Taylor1{T}`
polynomial, the accumulated (absolute or relative) error `Δ::Interval{S}`,
and `ξ0` which is
the mid point about which the expansion is obtained. If ξ0 is not exactly
representable, it returns *preferentiably* a rounded-down value.
This function is primarily used for plotting.

"""
function rpafp(tm::TMAbsRem{T}) where {T}
    fT = tm.pol
    Δ = remainder(tm)
    x0 = tm.x0
    ii = tm.iI
    order = get_order(tm)

    # α=0.484375 is used to get preferentially round-down of the mid point
    # when the mid point is not exactly representable
    α = 0.484375
    ξ0 = mid(x0, α)

    b = Taylor1(Interval{T}, order)
    t = Taylor1(T, order)
    for ind=0:order
        t[ind] = mid(fT[ind], α)
        b[ind] = fT[ind] - Interval(t[ind])
    end
    δ = b(ii-x0)
    Δ = Δ + δ
    return t, Δ, ξ0
end

function rpafp(tm::TMRelRem{T}) where {T}
    fT = tm.pol
    Δ = remainder(tm)
    x0 = tm.x0
    ii = tm.iI
    order = get_order(tm)

    # α=0.484375 is used to get preferentially round-down of the mid point
    # when the mid point is not exactly representable
    α = 0.484375
    ξ0 = mid(x0, α)

    b = Taylor1(Interval{T}, order)
    t = Taylor1(T, order)
    for ind=0:order
        t[ind] = mid(fT[ind], α)
        b[ind] = fT[ind] - Interval(t[ind])
    end
    δ = b(ii-x0)
    # Is the following correct for TMRelRem?
    return t, Δ, ξ0, δ
end


# Elementary functions
fnlist = (:inv, :sqrt, :exp, :log, :sin, :cos, :tan,
    :asin, :acos, :atan, :sinh, :cosh, :tanh)

for TM in tupleTMs
    for fn in fnlist
        @eval $fn(tm::$TM) = rpa($fn, tm)
    end
end