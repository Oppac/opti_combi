using JuMP, Cbc, Base

# generate random distance matrix
function generate_random_tsp(N, max_dist)
    srand()
    points = rand(0:max_dist, 1, N)
    c = Int64[norm(points[:,x] - points[:,y])
        for x in 1:N, y in 1:N]
    return c
end

# create the distance matrix
function generate_c(data)
    if isa(data, String)
        # see data_test.txt for correct format
        try
            c = readdlm("./$data", Int64)
            N = size(c, 1)
            return N, c
        catch e
            println("File must exist and be of the correct format")
            println("(see data_test.txt as an example)")
        end
    else
        # usage: max_dist = maximum possible distance
        N = data; max_dist = 25
        c = generate_random_tsp(N, max_dist)
        return N, c
    end
end

# compute the cycle taken
function solved(N, m, x, cycle)
    x_val = getvalue(x)
    push!(cycle, 1)
    while true
        v, idx = findmax(x_val[f=cycle[end],t=1:N])
        if idx == cycle[1]
            break
        else
            push!(cycle, idx)
        end
    end

    if length(cycle) < N
        @constraint(m, sum(x[f=cycle, t=cycle]) <= length(cycle)-1)
        return false
    end

    println("cycle: ", cycle)
    println("lenght: ", length(cycle))
    println("objective value: ", getobjectivevalue(m))
    #Base.showarray(STDOUT, val, false) #print xij matrix
    return true
end

function solve_tsp(N, c, p, mode)
    m = Model(solver = CbcSolver())

    # N^2 binary variables / xij = 1 if connect two cities, 0 otherwise
    # i city coming from / j city going to
    @variable(m, x[i=1:N, j=1:N], Bin)

    # minimize distance of selected connections
    @objective(m, Min, sum(c[i,j] * x[i,j]
               for i=1:N for j=1:N))

    # a city cannot connect with itself
    @constraint(m, [i=1:N], x[i,i] == 0)
    # enter and leave each city only once (cyclic path)
    @constraint(m, [i=1:N], sum(x[i, j=1:N]) == 1)
    @constraint(m, [j=1:N], sum(x[i=1:N, j]) == 1)

    if mode == "MTZ"
        # What is the "potential" of city i
    	@variable(m, u[1:N-1] >= 0)

        for i = 1:N-1, j = 1:N-1
            @constraint(m, u[i] - u[j] + p*x[i, j] <= p-1)
        end
    end

    if mode == "FCG"
                        
        @variable(m, y[i=1:N, j=1:N] >= 0)
        @variable(m, z[i=1:N, j=1:N] >= 0)

        @constraint(m, sum(y[1, j=1:N]) - sum(y[j=1:N, 1]) == N-1)
        @constraint(m, [i=2:N], sum(y[i, j=1:N]) - sum(y[j=1:N, i]) == -1)

        @constraint(m, sum(z[1, j=1:N]) - sum(z[j=1:N, 1]) == -(N-1))
        @constraint(m, [i=2:N], sum(z[i, j=1:N]) - sum(z[j=1:N, i]) == 1)

        @constraint(m, [i=1:N], sum(y[i, j=1:N]) + sum(z[i, j=1:N]) == N-1)

         for j = 1:N, i = 1:N
            @constraint(m, y[i, j] + z[i, j] == (N-1) * x[i, j])
        end
    end


    cycle = Array{Int}(0)
    status = solve(m)

    while !solved(m, x, cycle, N)
      status = solve(m)
    end
end


#-----------------------------------------------------#

data = readline(STDIN)
if parse(data) isa Number
    data = parse(Int64, data)
end

if isa(data, String) || isa(data, Int)
    # N = Number of cities
    # c = Distances between cities
    N, c = generate_c(data)
    Base.showarray(STDOUT, c, false) #print c in stdout
    println(" ")

    # p = Maximal number of cities that can be visited in one tour
    p = N-1

    mode = "MTZ"
    #mode = "FCG"

    solve_tsp(N, c, p, mode)
    #print(@time solve_tsp(N, c, p, mode))

else
    println("usage: text file matrix or number of cities")
end


#JuMP Citation
# @article{DunningHuchetteLubin2017,
# author = {Iain Dunning and Joey Huchette and Miles Lubin},
# title = {JuMP: A Modeling Language for Mathematical Optimization},
# journal = {SIAM Review},
# volume = {59},
# number = {2},
# pages = {295-320},
# year = {2017},
# doi = {10.1137/15M1020575},
# }
