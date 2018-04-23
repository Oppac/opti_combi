using JuMP, Cbc, Base

# generate random distance matrix
function generate_random_tsp(n, max_dist)
    srand()
    points = rand(0:max_dist, 1, n)
    c = Int64[norm(points[:,x] - points[:,y])
        for x in 1:n, y in 1:n]
    return c
end

# create the distance matrix
function generate_c(data)
    if isa(data, String)
        # see data_test.txt for correct format
        try
            c = readdlm("./$data", Int64)
            n = size(c, 1)
            return n, c
        catch e
            println("File must exist and be of the correct format")
            println("(see data_test.txt as an example)")
        end
    else
        # usage: max_dist = maximum possible distance
        n = data; max_dist = 25
        c = generate_random_tsp(n, max_dist)
        return n, c
    end
end

# compute the cycle taken
function solved(m, x, cycle)
    x_val = getvalue(x)
    push!(cycle, 1)
    while true
        v, idx = findmax(x_val[f=cycle[end],t=1:n])
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

    # n^2 binary variables / xij = 1 if connect two cities, 0 otherwise
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
        # what is the "potential" of city i
    	@variable(m, u[1:N-1] >= 0)

        for i = 1:N-1
            for j = 1:N-1
                @constraint(m, u[i] - u[j] + p*x[i, j] <= p - 1)
            end
        end
    end

    if mode == "FCG"

    end

    cycle = Array{Int}(0)
    status = solve(m)

    while !solved(m, x, cycle)
      status = solve(m)
    end
end


#-----------------------------------------------------#

data = readline(STDIN)
if parse(data) isa Number
    data = parse(Int64, data)
end

if isa(data, String) || isa(data, Int)
    # n = number of cities
    # c = distances between cities
    n, c = generate_c(data)
    #Base.showarray(STDOUT, c, false) #print c in stdout

    # p = maximal number of cities that can be visited in one tour
    p = n-1

    mode = "MTZ"
    #mode = "FCG"

    solve_tsp(n, c, p, mode)
    #print(@time solve_tsp(n, c, p, mode))

else
    println("usage: text file matrix or number of cities")
end
