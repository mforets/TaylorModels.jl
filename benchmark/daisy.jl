

using TaylorModels, BenchmarkTools, TaylorSeries
using IntervalArithmetic, IntervalOptimisation

# helper function that returns interval enclosures for the global minimum
# and global maximum of a univariate or multivariate polynomial over a given domain

SUITE["sin"] = BenchmarkGroup()
SUITE["bspline0"] = BenchmarkGroup()
SUITE["bspline1"] = BenchmarkGroup()
SUITE["bspline2"] = BenchmarkGroup()
SUITE["bspline3"] = BenchmarkGroup()
SUITE["kepler0"] = BenchmarkGroup()
SUITE["kepler1"] = BenchmarkGroup()
SUITE["kepler2"] = BenchmarkGroup()
SUITE["Doppler"] = BenchmarkGroup()
SUITE["himmilbeau"] = BenchmarkGroup()
SUITE["Rigidbody1"] = BenchmarkGroup()
SUITE["Rigidbody2"] = BenchmarkGroup()
SUITE["turbine1"] = BenchmarkGroup()
SUITE["turbine2"] = BenchmarkGroup()
SUITE["turbine3"] = BenchmarkGroup()


function _minmax(p, dom) #takes time
    global_min, _ = minimise(p, dom)
    minus_global_max, _ = minimise(-p, dom)
    global_max = -minus_global_max
    return global_min, global_max
end

#sin___________________________________________________#

Dx=Interval(-4.5,-0.3)
dom=Dx
x=Taylor1(7)

sin=x - (x*x*x)/6.0 + (x*x*x*x*x)/120.0 - (x*x*x*x*x*x*x)/5040.0

sin_out = evaluate(sin,dom)
SUITE["sin"]["evaluate"] = @benchmarkable evaluate($sin, $dom)

global_min, global_max = _minmax(sin, dom)
x_low      =sin_out.lo
x_high     =sin_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2
relative_precision_sin=Inerval(-(x_low-x_ref_low)/(x_high-x_ref_low)*100,
                              (x_high-x_ref_high)/(x_high-x_ref_low)*100)

#bspline0___________________________________________________#

Du=Interval(-4.5,-0.3)
dom=Du
u=Taylor1(3)

bspline0=(1 - u) * (1 - u) * (1 - u) / 6.0


bspline0_out = evaluate(bspline0, dom)
SUITE["bspline0"]["evaluate"] = @benchmarkable evaluate($bspline0, $dom)
global_min, global_max = _minmax(bspline0, dom)
x_low      =bspline0_out.lo
x_high     =bspline0_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2
relative_precision_bspline0=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                   (x_high-x_ref_high)/(x_high-x_ref_low)*100)


bspline1=(3 * u*u*u - 6 * u*u + 4) / 6.0


bspline1_out = evaluate(bspline1, dom)
SUITE["bspline1"]["evaluate"] = @benchmarkable evaluate($bspline1, $dom)
global_min, global_max = _minmax(bspline1, dom)
x_low      =bspline1_out.lo
x_high     =bspline1_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2
relative_precision_bspline1=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                   (x_high-x_ref_high)/(x_high-x_ref_low)*100)


bspline2=(-3 * u*u*u  + 3*u*u + 3*u + 1) / 6.0

bspline2 = evaluate(bspline2, dom)
SUITE["bspline2"]["evaluate"] = @benchmarkable evaluate($bspline2, $dom)
global_min, global_max = _minmax(bspline2, dom)
x_low      =bspline2_out.lo
x_high     =bspline2_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2
relative_precision_bspline2=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                   (x_high-x_ref_high)/(x_high-x_ref_low)*100)


bspline3=-u*u*u / 6.0

bspline3 = evaluate(bspline3, dom)
SUITE["bspline3"]["evaluate"] = @benchmarkable evaluate($bspline3, $dom)
global_min, global_max = _minmax(bspline3, dom)
x_low      =bspline3_out.lo
x_high     =bspline3_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2

relative_precision_bspline3=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                            (x_high-x_ref_high)/(x_high-x_ref_low)*100)

#Doppler___________________________________________________________________#

Du=Interval(-4.5,-0.3)
Dv=Interval(0.4,0.9)
Dt=Interval(3.8,7.8)
dom=Du×Dv×Dt
v,u,T=set_variables(Float64,["v","u","T"],order=4)

Doppler=(- ((331.4 + 0.6 * T)) *v) /
                            (((331.4 + 0.6 * T) + u)*((331.4 + 0.6 * T) + u))

v,u,T=set_variables(Float64,["v","u","T"],order=4)

Doppler_out = evaluate(Doppler, dom)
SUITE["Doppler"]["evaluate"] = @benchmarkable evaluate($Doppler, $dom)
global_min, global_max = _minmax(Doppler, dom)
x_low      =Doppler_out.lo
x_high     =Doppler_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2

relative_precision_Doppler=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                  (x_high-x_ref_high)/(x_high-x_ref_low)*100)

#himmilbeau___________________________________________________________________#

Dx1=Interval(-4.5,-0.3)
Dx2=Interval(0.4,0.9)
dom=Dx1×Dx2
x1,x2=set_variables(Float64,["x1","x2"],order=5)

himmilbeau=(x1*x1 + x2 - 11)*(x1 * x1 + x2 - 11)
                                        + (x1 + x2*x2 - 7)*(x1 + x2*x2 - 7)

himmilbeau_out = evaluate(himmilbeau, dom)
SUITE["himmilbeau"]["evaluate"] = @benchmarkable evaluate($himmilbeau, $dom)
global_min, global_max = _minmax(himmilbeau , dom)
x_low      =himmilbeau_out.lo
x_high     =himmilbeau_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2
relative_precision_himmilbeau=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                     (x_high-x_ref_high)/(x_high-x_ref_low)*100)

#kepler________________________________________________________________________#

Dx1=Interval(-4.5,-0.3)
Dx2=Interval(0.4,0.9)
Dx3=Interval(3.8,7.8)
Dx4=Interval(8.0,10.0)
Dx5=Interval(-10.0,8.0)
Dx6=Interval(1.0,2.0)
dom=Dx1×Dx2×Dx3×Dx4×Dx5×Dx6


x1,x2,x3,x4,x5,x6=set_variables(Float64,["x1","x2","x3","x4","x5","x6"],order=3)
kepler0=  x2 * x5 + x3 * x6 - x2 * x3 - x5 * x6
                                           + x1 * (-x1 + x2 + x3 - x4 + x5 + x6)

kepler0 = evaluate(kepler0, dom)
SUITE["kepler0"]["evaluate"] = @benchmarkable evaluate($kepler0, $dom)
global_min, global_max = _minmax(kepler0, dom)
x_low      =kepler0_out.lo
x_high     =kepler0_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2
relative_precision_Kepler0=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                  (x_high-x_ref_high)/(x_high-x_ref_low)*100)

dom=Dx1×Dx2×Dx3×Dx4
x1,x2,x3,x4=set_variables(Float64,["x1","x2","x3","x4"],order=3)
kepler1=  x1 * x4 * (-x1 + x2 + x3 - x4) + x2 * (x1 - x2 + x3 + x4)
          + x3 * (x1 + x2 - x3 + x4) -x2 * x3 * x4 - x1 * x3 - x1 * x2 - x4

kepler1_out = evaluate(kepler1, dom)
SUITE["kepler1"]["evaluate"] = @benchmarkable evaluate($kepler1, $dom)
global_min, global_max = _minmax(kepler1, dom)
x_low      =kepler1_out.lo
x_high     =kepler1_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2

relative_precision_Kepler1=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                  (x_high-x_ref_high)/(x_high-x_ref_low)*100)

dom=Dx1×Dx2×Dx3×Dx4×Dx5×Dx6
x1,x2,x3,x4,x5,x6=set_variables(Float64,["x1","x2","x3","x4","x5","x6"],order=3)
kepler2=  x1 * x4 * (-x1 + x2 + x3 - x4 + x5 + x6)
+ x2 * x5 * (x1 - x2 + x3 + x4 - x5 + x6) +x3* x6 * (x1 + x2 - x3 + x4 + x5 - x6)
- x2 * x3 * x4 -x1* x3* x5 - x1 * x2 * x6 - x4 * x5 * x6


kepler2_out = evaluate(kepler2, dom)
SUITE["kepler2"]["evaluate"] = @benchmarkable evaluate($kepler2, $dom)
global_min, global_max = _minmax(kepler2, dom)
x_low      =kepler2_out.lo
x_high     =kepler2_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2

relative_precision=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                            (x_high-x_ref_high)/(x_high-x_ref_low)*100)

#Rigidbody________________________________________________________#
Dx1=Interval(-4.5,-0.3)
Dx2=Interval(0.4,0.9)
Dx3=Interval(3.8,7.8)
dom=Dx1×Dx2×Dx3
x1,x2,x3=set_variables(Float64,["x1","x2","x3"],order=3)

Rigidbody1=-x1*x2 - 2*x2*x3 - x1 - x3

Rigidbody1_out = evaluate(Rigidbody1, dom)
SUITE["Rigidbody1"]["evaluate"] = @benchmarkable evaluate($Rigidbody1, $dom)
global_min, global_max = _minmax(Rigidbody1, dom)
x_low      =Rigidbody1_out.lo
x_high     =Rigidbody1_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2

relative_precision_Rigidbody1=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                     (x_high-x_ref_high)/(x_high-x_ref_low)*100)

x1,x2,x3=set_variables(Float64,["x1","x2","x3"],order=4)
Rigidbody2=2*(x1*x2*x3) + (3*x3*x3) - x2*(x1*x2*x3) + (3*x3*x3) - x2


Rigidbody2_out = evaluate(Rigidbody2, dom)
SUITE["Rigidbody2"]["evaluate"] = @benchmarkable evaluate($Rigidbody2, $dom)
global_min, global_max = _minmax(Rigidbody2, dom)
x_low      =Rigidbody2_out.lo
x_high     =Rigidbody2_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2

relative_precision_Rigidbody2=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                      (x_high-x_ref_high)/(x_high-x_ref_low)*100)

#turbine_________________________________________________________#

Dv=Interval(-4.5,-0.3)
Dw=Interval(0.4,0.9)
Dr=Interval(3.8,7.8)
dom=Dv×Dw×Dr
v,w,r=set_variables(Float64,["v","w","r"],order=6)
#=
turbine1= 3+ 2/(r*r) - (0.125*(3-2*v)*(w*w*r*r))/(1-v) - 4.5
turbine1_out = evaluate(turbine1, dom)
SUITE["turbine1"]["evaluate"] = @benchmarkable evaluate($turbine1, $dom)
global_min, global_max = _minmax(turbine1, dom)
x_low      =turbine1_out.lo
x_high     =turbine1_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2

relative_precision_turbine1=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                   (x_high-x_ref_high)/(x_high-x_ref_low)*100)
=#

turbine2=6*v - 0.5 * v * (w*w*r*r) / (1-v) - 2.5


turbine2_out = evaluate(turbine2, dom)
SUITE["turbine2"]["evaluate"] = @benchmarkable evaluate($turbine2, $dom)
global_min, global_max = _minmax(turbine2, dom)
x_low      =turbine2_out.lo
x_high     =turbine2_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2
relative_precision_turbine2=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                   (x_high-x_ref_high)/(x_high-x_ref_low)*100)

#=
turbine3= 3 - 2/(r*r) - 0.125 * (1+2*v) * (w*w*r*r) / (1-v) - 0.5
turbine3_out = evaluate(turbine3, dom)
SUITE["turbine3"]["evaluate"] = @benchmarkable evaluate($turbine3, $dom)
global_min, global_max = _minmax(turbine3, dom)
x_low      =turbine3_out.lo
x_high     =turbine3_out.hi
x_ref_low  =(global_min.hi+global_min.lo)/2
x_ref_high =(global_max.hi+global_max.lo)/2

relative_precision_turbine3=Interval(-(x_low-x_r_low)/(x_high-x_ref_low)*100,
                                    (x_high-x_ref_high)/(x_high-x_ref_low)*100)
=#
