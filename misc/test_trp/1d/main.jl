#!/usr/bin/env julia

include("../../../../misc/StochasticArithmetic.jl/src/StochasticArithmetic.jl")
using Main.StochasticArithmetic
using Random

function twoSum(a, b)
  x = a + b
  z = x - a
  y = (a - (x-z)) + (b-z)
  return (x,y)
end

function split(a)
  # New methods should be added to compute the correct
  # constant for each FP representation
  #
  # factor = 1+2^s, where s is the mantissa length
  factor(x::Float64) = 134217729
  factor(x::Float32) = 4097
  
  c = factor(a) * a
  x = c - (c - a)
  y = a - x
  return (x, y)
end

function twoProd(a, b)
  x = a * b
  (a1, a2) = split(a)
  (b1, b2) = split(b)
  y = a2*b2 - (((x - a1*b1) - a2*b1) - a1*b2)
  return (x,y)
end

function twoDiv(a, b)
  x = a / b
  (p1, p2) = twoProd(x, b)
  y = (a-p1-p2)/b
  return (x, y)
end


function round_dir(op, a, b, mode)
  ops = Dict(:+ => twoSum,
             :* => twoProd,
             :/ => twoDiv)
  mode_pred   = Dict(:up   => (err -> (err>0)),
                     :down => (err -> (err<0)))
  mode_change = Dict(:up   => nextfloat,
                     :down => prevfloat)
  (res,err) = ops[op](a,b)
  if (err != 0 && mode_pred[mode](err))
      res = mode_change[mode](res)
  end
  
  return res
end

round_rnd = let rng = MersenneTwister()
  function(op, a, b)
      return round_dir(op, a, b, rand(rng, [:up, :down]))    
  end
end

my_prod(x, y) = x * y
my_prod(x::Float64, y::Float64) = round_rnd(:*, x, y)
my_prod(x::Float32, y::Float32) = round_rnd(:*, x, y)

my_sum(x, y) = x + y
my_sum(x::Float64, y::Float64) = round_rnd(:+, x, y)
my_sum(x::Float32, y::Float32) = round_rnd(:+, x, y)

macro mca(defun)
    transform!(x) = Nothing()
    transform!(e::Expr) = begin
        for i in 1:length(e.args)
            if typeof(e.args[i]) == Expr && e.args[i].head == :(+=)
                ee = e.args[i]
                transform!.(ee.args)
                e.args[i] = :($(ee.args[1]) = $(ee.args[1]) + $(ee.args[2]))
            elseif typeof(e.args[i]) == Expr && e.args[i].head == :(*=)
                ee = e.args[i]
                transform!.(ee.args)
                e.args[i] = :($(ee.args[1]) = $(ee.args[1]) * $(ee.args[2]))
            elseif e.args[i] == :*
                e.args[i] = :my_prod
            elseif e.args[i] == :+
                e.args[i] = :my_sum
            end
            transform!(e.args[i])
        end
    end
        
    transform!(defun)
    println(defun)
    defun
end

struct Params{T}
  dt::T
  dx::T
  v::Vector{T}
end

function _unpack(p::Params{T}) where{T}
  return (p.dt,p.dx)
end
function export_data(filename::String,X::Vector{T},Y::Vector{U}) where{T,U}
  open(filename,"w") do f
    for i in 1:size(X,1)
      x = X[i]
      y = Y[i]
      write(f,"$x $y\n")
    end
  end
end
function export_data_tuple(filename::String,X::Vector{T},Y::Vector{U}) where{T,U}
  open(filename,"w") do f
    for i in 1:size(X,1)
      x = X[i]
      y1,y2 = Y[i]
      write(f,"$x $y1 $y2\n")
    end
  end
end

@mca function euler(un::Vector{T},L,tn::T,param::Params{T}) where{T}
  un1 = Vector{T}(undef,size(un))
  un1 = un - param.dt*L(un,tn,param)
  return un1
end

@mca function rk33(un::Vector{T},L,tn::T,param::Params{T}) where{T}
  dt = param.dt
  u1 = Vector{T}(undef,size(un))
  u2 = Vector{T}(undef,size(un))
  un1 = Vector{T}(undef,size(un))

  Lun = L(un,tn,param)
  for i in 1:size(un,1)
    u1[i] = un[i] - dt*Lun[i]
  end
  Lu1 = L(u1,tn+dt,param)
  for i in 1:size(un,1)
    u2[i] = 0.75*un[i] + 0.25*u1[i] - 0.25*(dt*Lu1[i])
  end
  Lu2 = L(u2,tn+0.5*dt,param)
  for i in 1:size(un,1)
    un1[i] = (1.0/3.0)*un[i] + (2.0/3.0)*u2[i] - (2.0/3.0)*(dt*Lu2[i])
  end
  return un1
end

@mca function cd2(un::Vector{T},t::T,param::Params{T}) where{T}
  dxu = Vector{T}(undef,size(un))

  dxu[1] = un[1+1]-un[size(un,1)]
  dxu[size(un,1)] = un[1]-un[size(un,1)-1]
  for i in 2:size(un,1)-1
    dxu[i] = un[i+1]-un[i-1]
  end

  return dxu/(2.0*param.dx)
end

@mca function flux(u::Vector{T},im2::Int,im1::Int,i::Int,ip1::Int,ip2::Int,ip3::Int) where{T}
  eps = 1e-6

  b0p = (13.0/12.0)*( u[im2] - 2.0*u[im1] + u[i]   )^2 + 0.25*( u[im2] - 4.0*u[im1] + 3.0*u[i] )^2
  b1p = (13.0/12.0)*( u[im1] - 2.0*u[i]   + u[ip1] )^2 + 0.25*( u[im1] - u[ip1] )^2
  b2p = (13.0/12.0)*( u[i]   - 2.0*u[ip1] + u[ip2] )^2 + 0.25*( 3.0*u[i] - 4.0*u[ip1] + u[ip2] )^2

  b0m = (13.0/12.0)*( u[ip1] - 2.0*u[ip2] + u[ip3] )^2 + 0.25*( 3.0*u[ip1] - 4.0*u[ip2] + u[ip3] )^2
  b1m = (13.0/12.0)*( u[i]   - 2.0*u[ip1] + u[ip2] )^2 + 0.25*( u[i] - u[ip2] )^2
  b2m = (13.0/12.0)*( u[im1] - 2.0*u[i]   + u[ip1] )^2 + 0.25*( u[im1] - 4.0*u[i] + 3.0*u[ip1] )^2

  a0p = 0.1/(eps+b0p)^2
  a1p = 0.6/(eps+b1p)^2
  a2p = 0.3/(eps+b2p)^2

  a0m = 0.1/(eps+b0m)^2
  a1m = 0.6/(eps+b1m)^2
  a2m = 0.3/(eps+b2m)^2

  sum_ap = a0p+(a1p+a2p)
  sum_am = a0m+(a1m+a2m)

  w0p = a0p/sum_ap
  w1p = a1p/sum_ap
  w2p = a2p/sum_ap

  w0m = a0m/sum_am
  w1m = a1m/sum_am
  w2m = a2m/sum_am

  uip12p = (w0p*( (2.0/6.0)*u[im2] + ((-7.0/6.0)*u[im1] + (11.0/6.0)*u[i]) ) + w1p*( (-1.0/6.0)*u[im1] + ((5.0/6.0)*u[i] + (2.0/6.0)*u[ip1]) )) + w2p*( (2.0/6.0)*u[i] + ((5.0/6.0)*u[ip1] + (-1.0/6.0)*u[ip2]) )
  uip12m = (w2m*( (-1.0/6.0)*u[im1] + ((5.0/6.0)*u[i] + (2.0/6.0)*u[ip1]) ) + w1m*( (2.0/6.0)*u[i] + ((5.0/6.0)*u[ip1] + (-1.0/6.0)*u[ip2]) )) + w0m*( (11.0/6.0)*u[ip1] + ((-7.0/6.0)*u[ip2] + (2.0/6.0)*u[ip3]) )

  return (uip12p,uip12m)
end

@mca function weno(un::Vector{T},t::T,param::Params{T}) where{T}
  dxu = Vector{T}(undef,size(un))
  vm = Vector{T}(undef,size(param.v))
  vp = Vector{T}(undef,size(param.v))
  for i in 1:size(param.v,1)
    vp[i] = max(param.v[i],0.0)
    vm[i] = min(param.v[i],0.0)
  end

  N = size(un,1)
  i=-1
  uim12p,uim12m = flux(un, ((i-2)+N)%N+1 , ((i-1)+N)%N+1 , (i+N)%N+1 , ((i+1)+N)%N+1 , (i+2)+1 , (i+3)+1 )
  for i in 0:2
    uip12p,uip12m = flux(un, ((i-2)+N)%N+1 , ((i-1)+N)%N+1 , (i+N)%N+1 , ((i+1)+N)%N+1 , (i+2)+1 , (i+3)+1 )
    dxu[i+1] = vp[i+1]*(uip12p-uim12p) + vm[i+1]*(uip12m-uim12m)
    uim12p = uip12p
    uim12m = uip12m
  end
  for i in 3:size(un,1)-4
    uip12p,uip12m = flux(un, (i-2)+1 , (i-1)+1 , (i)+1 , (i+1)+1 , (i+2)+1 , (i+3)+1 )
    dxu[i+1] = vp[i+1]*(uip12p-uim12p) + vm[i+1]*(uip12m-uim12m)
    uim12p = uip12p
    uim12m = uip12m
  end
  for i in size(un,1)-3:size(un,1)-1
    uip12p,uip12m = flux(un, (i-2)+1 , (i-1)+1 , (i)+1 , (i+1)%N+1 , (i+2)%N+1 , (i+3)%N+1 )
    dxu[i+1] = vp[i+1]*(uip12p-uim12p) + vm[i+1]*(uip12m-uim12m)
    uim12p = uip12p
    uim12m = uip12m
  end

  return dxu/param.dx
end

@mca function simu(time_method,space_method,Tf::T,u0::Vector{T},p::Params{T}) where{T}
  current_time = T(0.0)
  n_iter = 0
  un = Vector{T}(u0)
  while ( n_iter*p.dt < Tf )
    un = time_method(un,space_method,current_time,p)
    current_time += p.dt
    n_iter += 1
  end
  return un
end

function launcher(label::String,func,X::Vector{T},Tf::T,p::Params{T}) where{T}
  print("\033[94m>\033[0m $label\n")
  u0 = func(X)
  export_data("init_$label.dat",X,u0)

  print("  \033[95m+\033[0m cd2\n")
  un = @reliable_digits simu(rk33,cd2,Tf,Float64.(u0),p)
  export_data_tuple("end_cd2_$label.dat",X,un)

  print("  \033[95m+\033[0m weno\n")
  un = @reliable_digits simu(rk33,weno,Tf,Float64.(u0),p)
  export_data_tuple("end_weno_$label.dat",X,un)
end

function gate(X::Vector{T}) where{T}
  u = Vector{T}(undef,size(X))
  i = 1
  for x in X
    if x < 0.25 || x > 0.75
      u[i] = 0.0
    else
      u[i] = 1.0
    end
    i += 1
  end
  return u
end

function hat(X::Vector{T}) where{T}
  u = Vector{T}(undef,size(X))
  i = 1
  for x in X
    if x < 0.25 || x > 0.75
      u[i] = 0.0
    elseif x<0.5
      u[i] = 4*x-1
    else
      u[i] = -4*x+3
    end
    i += 1
  end
  return u
end

function sinus(X::Vector{T}) where{T}
  return sin.((2.0*π)*X)
end

function cosinus(X::Vector{T}) where{T}
  return cos.((2.0*π)*X)
end

function sinusHF(X::Vector{T}) where{T}
  return sin.(6*(2.0*π)*X)
end

function main()
  N = 100
  X = Vector{Float64}(Vector([0:N-1])[1]/N)
  dx = Float64(1.0/N)
  v = Vector{Float64}(fill(1.0,N))
  p = Params{Float64}(0.1*dx,dx,v)
  Tf = .1

  #--- cos -------------------------------
  launcher("cos",cosinus,X,Tf,p)
  #--- sin -------------------------------
  launcher("sin",sinus,X,Tf,p)
  #--- sinHF -------------------------------
  launcher("sinHF",sinusHF,X,Tf,p)
  #--- hat -------------------------------
  launcher("hat",hat,X,Tf,p)
  #--- gat -------------------------------
  launcher("gat",gate,X,Tf,p)
end

main()
