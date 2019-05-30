
using Plots

t_repr(s) = s^(3/4)
function t_repr(S,X,Patches, vpatche = 2)
    T = Float64[]
    for i in 1:length(S)
        A = (vpatche-Patches[X[i][1]+1,X[i][2]+1, X[i][3]+1])
        if A>(S[i])^3
            push!(T, t_repr(S[i]))
        else
            push!(T, Inf)
        end
    end
    return T
end
function vecinos(x, L)
    x1 = mod.(x+[0,0,1], L)
    x2 = mod.(x+[0,1,0], L)
    x3 = mod.(x-[0,0,1], L)
    x4 = mod.(x-[0,1,0], L)
    x5 = mod.(x+[1,0,0], L)
    x6 = mod.(x-[1,0,0], L)
    return [x1,x2,x3,x4,x5,x6]
end

function reproduccion_comida(Depredadores,Comida,r, ComidaMax = 2, DepredadorMax = 2 )
    n1,n2,n3 = size(Comida)
    Presa2 = zeros(n1,n2,n3)
    for i in 1:n1
        for j in 1:n2
            for k in 1:n3
                x = Comida[i,j,k]/ComidaMax
                if x < 1e-5
                    x = 1e-5
                end
                y = Depredadores[i,j,k]/DepredadorMax
                if y>1 
                    y = 1
                end
                x = r*x*(1-2y/(y+1))
                Presa2[i,j,k] = x*ComidaMax
                if Presa2[i,j,k]> ComidaMax
                    Presa2[i,j,k] = ComidaMax
                end
            end
        end
    end
    return Presa2
end

function condiciones_iniciales(N, L = 10, Maxsize = 0.5)
    Patches = zeros(L,L,L)
    X = Array{Int,1}[]
    S = Float64[]
    t_vida = zeros(N)
    t_muerte = Float64[]
    t_rep = zeros(N)
    for i in 1:N
        push!(X, rand(0:(L-1), 3))
        push!(S, rand()*Maxsize)
        push!(t_muerte, 10*S[end]^1/4)
        Patches[X[i][1]+1,X[i][2]+1, X[i][3]+1] += (S[i])^3
    end
    return X,S,t_vida,t_muerte,t_rep, Patches
end

function movimiento(X,S,t_vida,t_muerte, t_rep, Patches, Comida, time, L=20, r = 2, vpatche = 2)
    if length(X) == 0
        Comida = reproduccion_comida(Patches,Comida,r, vpatche, vpatche)
        return X,S,t_vida,t_muerte,t_rep, Patches, Comida, time+.1
    end
    t_reproduccion = t_repr(S,X, Patches, vpatche).-t_rep
    t_muertes = t_muerte-t_vida
    t_2rep = copy(t_reproduccion)
    t_2muertes = copy(t_muertes)    
    test = true
    t = 0
    t3 = (ceil(10*time)-10*time)/10.
    if t3 == 0
        t3 = 0.1
    end
    while test
        i1 = argmin(t_2rep)
        i2 = argmin(t_2muertes)
    #     @show Comida[X[i1][1]+1,X[i1][2]+1, X[i1][3]+1]> Patches[X[i1][1]+1,X[i1][2]+1, X[i1][3]+1]
    #     @show t_2muertes[i2], t_reproduccion[i1]
        if t3<t_2muertes[i2] && t3<t_2rep[i1]
            t = t3
            Comida = reproduccion_comida(Patches,Comida,r,vpatche,vpatche)
            break
        end 
        if t_2muertes[i2]<t_2rep[i1] 
            t = t_muertes[i2]
            if t<0
                t = 0
            end
            Patches[X[i2][1]+1,X[i2][2]+1, X[i2][3]+1] -= (S[i2])^3
            if Patches[X[i2][1]+1,X[i2][2]+1, X[i2][3]+1] < 0
                Patches[X[i2][1]+1,X[i2][2]+1, X[i2][3]+1]  = 0.
            end
            deleteat!(X, i2)
            deleteat!(S, i2)
            deleteat!(t_vida, i2)
            deleteat!(t_muerte, i2)
            deleteat!(t_rep, i2)
            test = false
        else
            t = t_2rep[i1]
            if t<0
                t = 0
            end
               
            if Comida[X[i1][1]+1,X[i1][2]+1, X[i1][3]+1]> Patches[X[i1][1]+1,X[i1][2]+1, X[i1][3]+1]
                push!(X,vecinos(X[i1], L)[rand([1,2,3,4])])
                push!(S, max(1e-10,S[i1]+(rand()-0.5)*1e-8))
                push!(t_vida, -t)
                push!(t_muerte, 10*S[end]^1/4)
                push!(t_rep, -t)
                t_rep[i1] = -t
                Patches[X[end][1]+1,X[end][2]+1, X[end][3]+1] += (S[end])^3
                test = false
            else
                t_2rep[i1] = Inf
            end
        end
    end
    t_vida = t_vida .+ t
    time += t
    t_rep = t_rep .+ t
    return X,S,t_vida,t_muerte,t_rep, Patches, Comida, time 
end

function dispersion(N, t_max, r = 1.2, L = 20, vpatche = 2., Maxsize = 2. )
    time = 0
    @show "hola"
    X,S,t_vida,t_muerte,t_rep, Patches = condiciones_iniciales(N, L, Maxsize)
    n1,n2,n3 = size(Patches)
    Comida = 1*ones(n1,n2,n3)
    contador = 0
    T = []
    NN = []
    contador2 = 0
    while time<t_max
        contador += 1
        if mod(contador, 1000) == 1
            histogram(S, nbins = 20, key = false, show = :ijulia, title = "t = $time, N = $(length(S))")
        end
        X,S,t_vida,t_muerte,t_rep, Patches, Comida, time = movimiento(X,S,t_vida,t_muerte, t_rep, Patches, Comida, time, L, r, vpatche)    
        push!(T, time)
        push!(NN,length(S))
        if mod(contador, 1000) == 500
    #        plot(T,NN, key = false, show = :ijulia, xlabel = "time", ylabel = "N")
        end
    end
    return X,S,t_vida,t_muerte,t_rep, Patches, Comida, NN, T
end

function dispersion(X,S,t_vida,t_muerte,t_rep, Patches, Comida, t_max, NN, T, r = 1.2, L = 20, vpatche = 2)
    time = 0
    contador = 0
    com = []
    contador2 = 0
    while time<t_max
        contador += 1
        if mod(contador, 1000) == 1
            histogram(S, nbins = 30, key = false, show = :ijulia, title = "t = $time, N = $(length(S))")
        end
        X,S,t_vida,t_muerte,t_rep, Patches, Comida, time = movimiento(X,S,t_vida,t_muerte, t_rep, Patches, Comida, time, L, r, vpatche)
        push!(T, time)
        push!(NN,length(S))
        push!(com, sum(Comida))
        if length(S)<10
            break
        end
       # @show length(S), time
        if mod(contador, 1000) == 500
 #           plot(T,NN, key = false, show = :ijulia, xlabel = "time", ylabel = "N")
  #          plot(T,com, key = false, show = :ijulia, xlabel = "time", ylabel = "comida")
        end
    end
    return X,S,t_vida,t_muerte,t_rep, Patches, Comida, NN, T
end

function promedio(T, N, dT)
    Np = []
    Tp = []
    Tmax = maximum(T)
    t = 0
    NNp = 0
    contador = 1 
    while t< Tmax
        t += dT
        push!(Tp, t)
        contador2 = contador
        while T[contador]<t
            NNp += N[contador]
            contador += 1
        end
        dN = contador-contador2
        push!(Np, NNp/dN)
        NNp = 0
    end
    return Tp, Np
end     
