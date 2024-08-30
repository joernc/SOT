# reference time
global const tref = DateTime(2000, 1, 1, 0, 0, 0)

"""
    BFGS(f,x0,∇tol,maxiter)

Minimize f with initial guss x0, gradient tolerence ∇tol, and maximum iteration 
number maxiter using the BFGS scheme. 
"""
function BFGS(f,x0,∇tol,maxiter; xrtol=0)

    # dimension of problem 
    d = length(x0) 
    
    # initial gradient 
    old_fval,∇fk = f(x0) 
    
    # initial inverse hessian approx
    Hk = I(d)
    
    # initial search direction
    @printf("-H∇f = %s\n", -Hk*∇fk)
    
    # norm of gradient
    ∇norm = norm(∇fk, Inf)

    # sets the initial step guess to dx ~ 1
    old_old_fval = old_fval + norm(∇fk) / 2

    # initialize iteration
    xk = x0[:]

    # initialize step
    k = 0 

    # iteration
    while ∇norm > ∇tol && k < maxiter
    
        # print at each iteration
        @printf("k = %d, ek = %.2e\n", k, ∇norm)
        @printf("xk = %s\n", xk)
        
        # break for maximum iterations
        if k > maxiter
            @printf("Maximum iterations reached!\n")
            break
        end
        
        # search direction (Newton Method)
        pk = -Hk*∇fk 
        
        # initialize line search
        αk, ∇fkp1 = NaN,∇fk[:]
        try
        
            # line search with strong wolfe condition
            αk, old_fval, old_old_fval, ∇fkp1 = line_search_swolfe(f, xk, pk; ∇fk,old_fval, old_old_fval, amax=1e100)
        catch e
        
            # ine search failed to find a better solution.
            @printf("Line search failed to find a better solution...\n")
            rethrow(e)
           break
        end

        #αk = line_search(f,∇f,xk,pk,∇fk) # line search 
        
        # print line search parameter
        @printf("αk = %s\n", αk)
        if isnan(αk)
            break
        end
        
        # update line search output
        sk = vec(αk * pk)
        xkp1 = xk .+ sk

        # update iteration variable, function, and gradient
        xk = xkp1[:]
        yk = vec(∇fkp1 .- ∇fk)
        ∇fk = ∇fkp1[:] 
        
        # update step and gradient norm
        k += 1
        ∇norm = norm(∇fk, Inf)
        
        # break for convergence
        if ∇norm ≤ ∇tol
            break
        end

        # break for relative covergence
        if αk*norm(pk) ≤ xrtol*(xrtol + norm(xk))
            @printf("Relative tolerence reached to terminate\n")
            @printf("k = %d, ek = %.2e\n", k, ∇norm)
            break
        end

        # hessian update factor 1
        ρk_inv = yk⋅sk
        if ρk_inv == 0.
            ρk = 1000.0
            @printf("Divide-by-zero encountered: ρk assumed large\n")
        else
            ρk = 1 / ρk_inv
        end

        # hessian update factor 2
        A1 = I(d)-ρk*(sk*yk')
        A2 = I(d)-ρk*(yk*sk')
        if k < 1
             Hk = ρk_inv/(yk⋅yk)*I(d)
        end
        
        # BFGS hessian update
        Hk = A1*Hk*A2 + ρk*(sk*sk') 
    end
    
    return xk,Hk
end

function getE(tpairs,ppairs,tstation,pstations;hydro=false,bathy=false)
    l = 1

    # number of T- and P-wave pairs
    nt = size(tpairs, 1)
    np = size(ppairs, 1)
    # find unique events
    t = sort(unique([tpairs.event1; tpairs.event2]))

    # number of unique events
    m = length(t)

    # real time (days)
    tr = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24

    # T-wave pair matrix
    tidx1 = indexin(tpairs.event1, t)
    tidx2 = indexin(tpairs.event2, t)
    Xt = sparse([1:nt; 1:nt], [tidx1; tidx2], [-ones(nt); ones(nt)])

    # P-wave pair matrix
    pidx1 = indexin(ppairs.event1, t)
    pidx2 = indexin(ppairs.event2, t)
    Xp = sparse([1:np; 1:np], [pidx1; pidx2], [-ones(np); ones(np)])

    selpairs = select(ppairs, [:station, :event1, :event2])
    selpairs.avglat = 0.5(ppairs.latitude1+ppairs.latitude2)
    selpairs.avglon = 0.5(ppairs.longitude1+ppairs.longitude2)
    leftjoin!(selpairs, pstations, on=:station)
    ppazms = azimuth.(selpairs.slat,selpairs.slon,selpairs.avglat,selpairs.avglon)
    Xpcs = sparse([1:np; 1:np], [pidx1; pidx2], [-cosd.(ppazms); cosd.(ppazms)])
    Xpsn = sparse([1:np; 1:np], [pidx1; pidx2], [-sind.(ppazms); sind.(ppazms)])
    selpairs = select(selpairs, [:event1, :event2, :avglat, :avglon])
    unique!(selpairs)
    tpairs = select(tpairs, [:event1, :event2])
    leftjoin!(tpairs, selpairs, on=[:event1, :event2])
    tpazms = azimuth.(tstation[1],tstation[2],tpairs.avglat,tpairs.avglon)
    Xtcs = sparse([1:nt; 1:nt], [tidx1; tidx2], [-cosd.(tpazms); cosd.(tpazms)])
    Xtsn = sparse([1:nt; 1:nt], [tidx1; tidx2], [-sind.(tpazms); sind.(tpazms)])
    Ecs,Esn = [Xtcs;Xpcs],[Xtsn;Xpsn]

    # full design matrix
    E = blockdiag([Xt for i = 1:l]..., Xp)
    tm = tr[1]+(tr[m]-tr[1])/2
    E = [E [kron(I(l), Xt*(tr.-tm)); zeros(np, l)]]
    ω = 2π/SOT.meanyear
    E = [E [kron(I(l), Xt*cos.(ω*tr)) kron(I(l), Xt*sin.(ω*tr)); zeros(np, 2l)]]
    ω = 4π/SOT.meanyear
    E = [E [kron(I(l), Xt*cos.(ω*tr)) kron(I(l), Xt*sin.(ω*tr)); zeros(np, 2l)]]
    #D = [I(nt) -(E[1:nt,1:m]*pinv(Array(E[l*nt+1:l*nt+np,l*m+1:(l+1)*m])))]
    if hydro
        #Ih = Diagonal([ones(nt);zeros(np)])
        Ecsh,Esnh = [Xtcs;zeros(size(Xp))],[Xtsn;zeros(size(Xp))]
        return t,E,Ecs,Esn,Ecsh,Esnh
    elseif bathy 
        aref = azimuth.(tstation[1],tstation[2],38.10, 142.85)+1
        nb = cosd.(tpazms.-aref) 
        #nb = exp.(-abs.(tpazms.-aref))
        Xtb = sparse([1:nt; 1:nt], [tidx1; tidx2], [-nb; nb])
        Eb = [Xtb;zeros(size(Xp))]
        return t,E,Ecs,Esn,Eb
    else
        return t,E,Ecs,Esn
    end
end


function loglikelihood(x, y, σtrend, σannual, σsemiannual, t, θ, E, Ecs, Esn, nt, np;nλt=1,cgrid=0,θgrid=0,stn6=0,Ecsh=I,Esnh=I,Eb=I,hydro=false,bathy=false,grad=true)

    σtrend /= meanyear

    l = 1
    m = length(t)
    # real time (days)
    trd = Dates.value.(t - DateTime(2000, 1, 1, 12, 0, 0))/1000/3600/24
    if nλt == 1
        λt,λθ,στ,σp,σx,σn = exp.(x[1:6])
        σnp = exp.(x[7])
        if length(stn6) > 1
            σnp2 = exp.(x[8])
        end
        # solution covariance in time
        if length(cgrid) == 1
            C = exp.(-abs.(trd.-trd')/λt-(θ.-θ').^2/λθ.^2/2)
        else
            interp_linear = linear_interpolation(θgrid, cgrid, extrapolation_bc=0)
            C = exp.(-abs.(trd.-trd')/λt).*interp_linear.(abs.((θ.-θ')))
        end
        C *= στ^2/C[1,1]
    else
        λta,λtb,λθ = exp.(x[1:3])
        στ,wa,σp,σx,σn = exp.(x[4:8]) 
        Ca,Cb = exp.(-((trd.-trd')/λta).^2/2),exp.(-((trd.-trd')/λtb).^2/2)
        Cθ = exp.(-((θ.-θ')/λθ).^2/2)
        C = στ^2*(wa*Ca.+(1-wa)*Cb).*Cθ
    end    

    # covariance matrix assuming no correlation between singular vector expansion coefficients
    R = [C zeros(l*m, m); zeros(m, (l+1)*m)] + σp^2*kron(sparse(ones(l+1, l+1)), I(m))
    R = [R zeros(size(R, 1), l); zeros(l, size(R, 2)) (σtrend.^2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σannual.^2/2)]
    R = [R zeros(size(R, 1), 2l); zeros(2l, size(R, 2)) I(2)*(σsemiannual.^2/2)]

    Rcsn = Ecs*Ecs'+Esn*Esn'

    # noise covariance
    if hydro
        σh = exp(x[end])
        Ch = I
        Rcsnh = Ecsh*Ch*Ecsh'+Esnh*Ch*Esnh'
        N = σn^2*spdiagm([ones(nt);zeros(np)])+σx^2*Rcsn+σh^2*Rcsnh
        if length(stn6) == 1
            N += σnp^2*spdiagm([zeros(nt);ones(np)])
        else
            N += σnp^2*spdiagm([zeros(nt);abs.(stn6.=="IU.MAJ")])
            N += σnp2^2*spdiagm([zeros(nt);abs.(stn6.!="IU.MAJ")])
        end
    elseif bathy
        σb = exp.(x[7])
        #Cδ = σδ^2*exp.(-(θ.- 1)*(θ.-1)'/λδ^2)
        Rb = Eb*Eb'
        N = σn^2*I+σx^2*Rcsn+σb^2*Rb
    else
        N = σn^2*spdiagm([ones(nt);zeros(np)])+σx^2*Rcsn
        N += σnp^2*spdiagm([zeros(nt);ones(np)])
    end
  
  #try
    Ryy = E*R*E' + N
    #Ryy = D*Ryy*D'
    iRyy = inv(Ryy)
    ll = -0.5*(logdet(Ryy)+y'*iRyy*y+length(y)*log(2π))
    if grad
        d = length(x)
        Δll = zeros(d)
        α = iRyy*y
        δll(α,iRyy,δRyy) = 0.5*tr((α*α'-iRyy)*δRyy)
        #δll(α,iRyy,δRyy) = 0.5*tr((α*α'-iRyy)*D*δRyy*D')
        if nλt == 1
            δR1 = E*[(abs.(trd.-trd')/λt).*C zeros(l*m, m+5); zeros(m+5, (l+1)*m+5)]*E'
            Δll[1] = δll(α,iRyy,δR1)
            δR2 = E*[((θ.-θ').^2/λθ.^2).*C zeros(l*m, m+5); zeros(m+5, (l+1)*m+5)]*E'
            Δll[2] = δll(α,iRyy,δR2)
            δR3 = E*[2*C zeros(l*m, m+5); zeros(m+5, (l+1)*m+5)]*E'
            Δll[3] = δll(α,iRyy,δR3)
            δR4 = E*[2*σp^2*kron(sparse(ones(l+1, l+1)), I(m)) zeros((l+1)*m, 5l); zeros(5l, (l+1)*m+5)]*E'
            Δll[4] = δll(α,iRyy,δR4)
            δR5 = 2*σx^2*Rcsn
            Δll[5] = δll(α,iRyy,δR5)
            δR6 = 2*σn^2*spdiagm([ones(nt);zeros(np)])
            Δll[6] = δll(α,iRyy,δR6)
            if length(stn6) == 1
                δR7 = 2*σnp^2*spdiagm([zeros(nt);ones(np)])
            else
                δR7 = 2*σnp^2*spdiagm([zeros(nt);abs.(stn6.=="IU.MAJ")])
                δR8 = 2*σnp2^2*spdiagm([zeros(nt);abs.(stn6.!="IU.MAJ")])
                Δll[8] = δll(α,iRyy,δR8)
            end
            Δll[7] = δll(α,iRyy,δR7)
        else
            δR1 = E*[στ^2*wa*(((trd.-trd')/λta).^2).*Ca.*Cθ zeros(l*m, m+1); zeros(m+1, (l+1)*m+1)]*E'
            δR2 = E*[στ^2*(1-wa)*(((trd.-trd')/λtb).^2).*Cb.*Cθ zeros(l*m, m+1); zeros(m+1, (l+1)*m+1)]*E'
            δRθ = E*[((θ.-θ').^2/λθ.^2).*C zeros(l*m, m+1); zeros(m+1, (l+1)*m+1)]*E'
            Δll[1] = δll(α,iRyy,δR1)
            Δll[2] = δll(α,iRyy,δR2)
            Δll[3] = δll(α,iRyy,δRθ)
            δRτ = E*[2*C zeros(l*m, m+1); zeros(m+1, (l+1)*m+1)]*E'
            Δll[4] = δll(α,iRyy,δRτ)
        end
        
        #δRw = E*[στ^2*wa*(Ca.-Cb).*Cθ zeros(l*m, m+1); zeros(m+1, (l+1)*m+1)]*E'
        #Δll[3] = δll(α,iRyy,δRw)
        #δR5 = 2*σy^2*RcsEa
        #Δll[5] = δll(α,iRyy,δR5)
        #δR6 = -σxy*RcxEa
        #Δll[5] = δll(α,iRyy,δR6)
        
        if hydro & nλt==1
            δRh = 2*exp(x[end])^2*Rcsnh
            Δll[end] = δll(α,iRyy,δRh)
            #δRht = Ecsh*((abs.(trd.-trd')/λh).*Ch)*Ecsh'+Esnh*((abs.(trd.-trd')/λh).*Ch)*Esnh'
            #Δll[8] = δll(α,iRyy,σh^2*δRht)
        elseif bathy
            δR7 = 2*exp(x[end])^2*Rb
            Δll[7] = δll(α,iRyy,δR7)
            #δR8 = Eb*2*((θ.- 1)*(θ.-1)'/λδ.^2 .*Cδ)*Eb'
            #Δll[8] = δll(α,iRyy,δR8)
        end
        
        return ll,Δll
    else
        return ll
    end
  #catch y
  #  @printf("singular covariance\n")
  #  rethrow(y)
  #  return NaN
  #end
end

"""
lon. and lat. of n equally spaced points along the great circle path that 
originates at p1d = (λ1, θ1) and passes through p2d = (λ2, θ2)
"""
function findpath(p1d, p2d, n)
  λ1 = deg2rad(p1d[1])
  θ1 = deg2rad(p1d[2])
  λ2 = deg2rad(p2d[1])
  θ2 = deg2rad(p2d[2])
  λ2p = atan(cos(θ2)*sin(λ2-λ1), sin(θ1)*cos(θ2)*cos(λ2-λ1) - cos(θ1)*sin(θ2))
  θ2p = asin(cos(θ1)*cos(θ2)*cos(λ2-λ1) + sin(θ1)*sin(θ2))
  R2 = [cos(λ2p) -sin(λ2p) 0; sin(λ2p) cos(λ2p) 0; 0 0 1]
  Ry = [sin(θ1) 0 cos(θ1); 0 1 0; -cos(θ1) 0 sin(θ1)]
  Rz = [cos(λ1) -sin(λ1) 0; sin(λ1) cos(λ1) 0; 0 0 1]
  θp = π/2 .- (0:n-1)/(n-1)*(π/2 - θ2p)
  p = [cos.(θp)'; zeros(n)'; sin.(θp)']
  q = Rz*Ry*R2*p
  λq = atan.(q[2,:], q[1,:])
  θq = asin.(q[3,:])
  return rad2deg.(λq), rad2deg.(θq)
end


function exp_and_normalise(lw)
    w = exp.(lw .- max(lw...))
    return w ./ sum(w)
end

"""
credit: GeoMapping.jl

    az = azimuth(lat1,lon1,lat2,lon2)
Compute azimuth, i.e. the angle at (`lat1`,`lon1`) between the point (`lat2`,`lon2`) and the North, counted clockwise starting from the North.
The units of all input and output parameters are degrees.
```
          North
            ↑
            | .
            |   . az
(lat1,lon1) +   .
             ╲ ↙
              ╲
               ╲
                * (lat2,lon2)
```
"""

function azimuth(lat1,lon1,lat2,lon2)
    # azimuth(stalat,stalon,evtlat,evtlon)
    # https://en.wikipedia.org/w/index.php?title=Azimuth&oldid=750059816#Calculating_azimuth

    Δλ = π/180 * (lon2 - lon1)
    ϕ1 = π/180 * lat1
    ϕ2 = π/180 * lat2

    α = atan(sin(Δλ), cos(ϕ1)*tan(ϕ2) - sin(ϕ1)*cos(Δλ))
    return 180/π * α
end

"""
Compute azimuth for pair catalog
"""
function getazimuth(tpairs,ppairs,stalon,stalat,evtlon,evtlat)
    θ0 = azimuth(stalat,stalon,evtlat,evtlon)
  
    uevent = unique(vcat(rename(ppairs[:,2:4],[:event,:latitude,:longitude]),rename(ppairs[:,7:9],[:event,:latitude,:longitude])))
    utime = sort(unique([tpairs.event1; tpairs.event2]))
    
    # number of unique events
    m = length(utime)
    θevent = Array{Float64}(undef, m)
    xevent = Array{Float64}(undef, m)
    yevent = Array{Float64}(undef, m)
    @printf("T-wave event number: %d, pair number: %d\n", m,size(tpairs, 1))
    
    # iterate over tpairs
    for i = 1:m
      eidx = findfirst(x -> x==utime[i],uevent.event)
      xevent[i],yevent[i] = uevent.longitude[eidx],uevent.latitude[eidx]
      θevent[i] = azimuth(stalat,stalon,uevent.latitude[eidx],uevent.longitude[eidx]) - θ0
    end

    return xevent,yevent,θevent
end

# CENTRAL FINITE DIFFERENCE CALCULATION
function grad(f,x)
    h = cbrt(eps())
    d = length(x)
    nabla = zeros(d)
    for i = 1:d 
        x_for = copy(x) 
        x_back = copy(x)
        x_for[i] += h 
        x_back[i] -= h 
        nabla[i] = (f(x_for) - f(x_back))/(2*h) 
        if isnan(nabla[i])
            @printf("NaN gradient\n")
            return nabla*NaN
        end
    end
    return nabla 
end

# CENTRAL FINITE DIFFERENCE CALCULATION
function hess(f,x)
    h = cbrt(eps())
    d = length(x)
    hess = zeros(d,d)
    for i = 1:d 
        x_for = copy(x) 
        x_back = copy(x)
        x_for[i] += h 
        x_back[i] -= h 
        _,∇f_for = f(x_for)
        _,∇f_back = f(x_back)
        hess[:,i] = (∇f_for - ∇f_back)/(2*h)
        if isnan(hess[i,i])
            @printf("NaN hessian\n")
            return hess*NaN
        end
    end
    return 0.5*(hess + hess') 
end

# BACKTRACK LINE SEARCH WITH WOLFE CONDITIONS
function line_search(f,∇f,x,p,∇)
    αi = 1
    ϕ(α) = f(x+α*p)
    dϕ(α) = ∇f(x+α*p)⋅p
    ϕ0 = ϕ(0)
    dϕ0 = ∇⋅p
    ϕi,dϕi = ϕ(αi),dϕ(αi)
    c1 = 1e-4 
    c2 = 0.9 

    while ϕi > ϕ0 + (c1*αi*dϕ0) || dϕi < c2*dϕ0
        αi *= 0.5
        ϕi,dϕi = ϕ(αi),dϕ(αi)
    end
    return αi
end
    
"""
The following functions are Julia implementation of those from Scipy:
https://github.com/scipy/scipy/blob/v1.13.1/scipy/optimize/_linesearch.py#L187-L327
"""
function line_search_swolfe(f, xk, pk; ∇fk=NaN, old_fval=NaN, old_old_fval=NaN, c1=1e-4, c2=0.9, amax=NaN,maxiter=10)
    #fc = [0]
    #gc = [0]
    gval = NaN

    function ϕ(α)
        #fc[1] += 1
        fval,_ = f(xk + α * pk)
        return fval
    end

    function derϕ(α)
        #gc[1] += 1
        _,gval = f(xk + α * pk)
        return dot(gval, pk)
    end

    derϕ0 = dot(∇fk, pk)

    α_star, ϕ_star, old_fval, derϕ_star = scalar_search_swolfe(ϕ, derϕ; ϕ0=old_fval, old_ϕ0=old_old_fval, derϕ0, c1, c2, amax, maxiter=maxiter)

    if isnan(derϕ_star)
        @printf("The line search algorithm did not converge\n")
    else
        # derϕ_star is a number (derϕ) -- so use the most recently
        # calculated gradient used in computing it derϕ = ∇fk*pk
        # this is the gradient at the next step no need to compute it
        # again in the outer loop.
        derϕ_star = gval
    end

    return α_star, ϕ_star, old_fval, derϕ_star
end

"""
reference: Scipy
"""
function scalar_search_swolfe(ϕ, derϕ; ϕ0=NaN, old_ϕ0=NaN, derϕ0=NaN, c1=1e-4, c2=0.9, amax=NaN, maxiter=10)
    if isnan(ϕ0)
        ϕ0 = ϕ(0.)
    end

    if isnan(derϕ0)
        derϕ0 = derϕ(0.)
    end

    α0 = 0
    if !isnan(old_ϕ0) && (derϕ0 != 0)
        α1 = min(1.0, 1.01*2*(ϕ0 - old_ϕ0)/derϕ0)
    else
        α1 = 1.0
    end

    if α1 < 0
        α1 = 1.0
    end

    if !isnan(amax)
        α1 = min(α1, amax)
    end

    ϕ_α1 = ϕ(α1)

    ϕ_α0 = ϕ0
    derϕ_α0 = derϕ0

    extra_condition(α, ϕ) = true

    i = 1
    while i ≤ maxiter
        if α1 == 0 || (!isnan(amax) && α0 == amax)
            # α1 == 0: This shouldn't happen. Perhaps the increment has
            # slipped below machine precision?
            α_star = NaN
            ϕ_star = ϕ0
            ϕ0 = old_ϕ0
            derϕ_star = NaN

            if α1 == 0
                @printf("Rounding errors prevent the line search from converging\n")
            else
                @printf("The line search algorithm could not find a solution ≤ amax: %s\n", amax)
            end

            break
        end
        not_first_iteration = i > 1
        if (ϕ_α1 > ϕ0 + c1 * α1 * derϕ0) || ((ϕ_α1 >= ϕ_α0) && not_first_iteration)
            α_star, ϕ_star, derϕ_star = zoom(α0, α1, ϕ_α0, ϕ_α1, derϕ_α0, ϕ, derϕ, ϕ0, derϕ0, c1, c2, extra_condition)
            break
        end
        derϕ_α1 = derϕ(α1)
        if (abs(derϕ_α1) ≤ -c2*derϕ0)
            if extra_condition(α1, ϕ_α1)
                α_star = α1
                ϕ_star = ϕ_α1
                derϕ_star = derϕ_α1
                break
            end
        end
        if (derϕ_α1 ≥ 0)
            α_star, ϕ_star, derϕ_star = zoom(α1, α0, ϕ_α1,ϕ_α0, derϕ_α1, ϕ, derϕ, ϕ0, derϕ0, c1, c2, extra_condition)
            break
        end
        α2 = 2 * α1  # increase by factor of two on each iteration
        if !isnan(amax)
            α2 = min(α2, amax)
        end
        α0 = copy(α1)
        α1 = copy(α2)
        ϕ_α0 = copy(ϕ_α1)
        ϕ_α1 = ϕ(α1)
        derϕ_α0 = copy(derϕ_α1)

        i += 1
    end

    if i > maxiter
        # stopping test maxiter reached
        α_star = copy(α1)
        ϕ_star = copy(ϕ_α1)
        derϕ_star = NaN
        @printf("Maximum Wolfe iteration reached\n")
        @printf("The line search algorithm did not converge\n")
    end
    
    return α_star, ϕ_star, ϕ0, derϕ_star
end

"""
reference: Scipy
"""
function zoom(α_lo, α_hi, ϕ_lo, ϕ_hi, derϕ_lo, ϕ, derϕ, ϕ0, derϕ0, c1, c2, extra_condition)
    maxiter = 10
    i = 0
    delta1 = 0.2  # cubic interpolant check
    delta2 = 0.1  # quadratic interpolant check
    ϕ_rec = ϕ0
    a_rec = 0
    a_star, val_star, valprime_star = NaN,NaN,NaN
    while true
        # interpolate to find a trial step length between α_lo and
        # α_hi Need to choose interpolation here. Use cubic
        # interpolation and then if the result is within delta *
        # dalpha or outside of the interval bounded by α_lo or α_hi
        # then use quadratic interpolation, if the result is still too
        # close, then use bisection

        dalpha = α_hi - α_lo
        if dalpha < 0 
            a, b = α_hi, α_lo
        else 
            a, b = α_lo, α_hi
        end
        # minimizer of cubic interpolant
        # (uses ϕ_lo, derϕ_lo, ϕ_hi, and the most recent value of ϕ)
        #
        # if the result is too close to the end points (or out of the
        # interval), then use quadratic interpolation with ϕ_lo,
        # derϕ_lo and ϕ_hi if the result is still too close to the
        # end points (or out of the interval) then use bisection

        if (i > 0) 
            cchk = delta1 * dalpha
            αj = cubicmin(α_lo, ϕ_lo, derϕ_lo, α_hi, ϕ_hi, a_rec, ϕ_rec)
        end
        if (i == 0) || isnan(αj) || (αj > b - cchk) || (αj < a + cchk) 
            qchk = delta2 * dalpha
            αj = quadmin(α_lo, ϕ_lo, derϕ_lo, α_hi, ϕ_hi)
            if isnan(αj) || (αj > b-qchk) || (αj < a+qchk) 
                αj = α_lo + 0.5*dalpha
            end
        end
αj
        # Check new value of αj

        ϕ_αj = ϕ(αj)
        if (ϕ_αj > ϕ0 + c1*αj*derϕ0) || (ϕ_αj >= ϕ_lo) 
            ϕ_rec = copy(ϕ_hi)
            a_rec = copy(α_hi)
            α_hi = copy(αj)
            ϕ_hi = copy(ϕ_αj)
        else 
            derϕ_αj = derϕ(αj)
            if abs(derϕ_αj) ≤ -c2*derϕ0 && extra_condition(αj, ϕ_αj) 
                a_star = copy(αj)
                val_star = copy(ϕ_αj)
                valprime_star = copy(derϕ_αj)
                break
            end
            if derϕ_αj*(α_hi - α_lo) ≥ 0 
                ϕ_rec = copy(ϕ_hi)
                a_rec = copy(α_hi)
                α_hi = copy(α_lo)
                ϕ_hi = copy(ϕ_lo)
            else 
                ϕ_rec = copy(ϕ_lo)
                a_rec = copy(α_lo)
            end
            α_lo = copy(αj)
            ϕ_lo = copy(ϕ_αj)
            derϕ_lo = copy(derϕ_αj)
        end
        i += 1
        if (i > maxiter) 
            # Failed to find a conforming step size
            @printf("Maximum iteration reached\n")
            @printf("Failed to find a conforming step size\n")
            a_star = NaN
            val_star = NaN
            valprime_star = NaN
            break
        end
    end
    return a_star, val_star, valprime_star
end

"""
reference: Scipy
"""
function cubicmin(a, fa, fpa, b, fb, c, fc)
    xmin = NaN
    try
        C = fpa
        db = b - a
        dc = c - a
        denom = (db * dc)^2 * (db - dc)
        d1 = [dc^2 -db^2; -dc^3 db^3]
        A, B = d1*vec([fb - fa - C * db, fc - fa - C * dc])
        A /= denom
        B /= denom
        radical = B * B - 3 * A * C
        xmin = a + (-B + sqrt(radical)) / (3 * A)
    catch y
        return NaN
    end

    if !isfinite(xmin)
        return NaN
    end

    return xmin
end

"""
reference: Scipy
"""
function quadmin(a, fa, fpa, b, fb)
    xmin = NaN
    try
        D = fa
        C = fpa
        db = b - a * 1.0
        B = (fb - D - C * db) / (db * db)
        xmin = a - C / (2.0 * B)
    catch y
        return NaN
    end
    if !isfinite(xmin)
        return NaN
    end
    return xmin
end