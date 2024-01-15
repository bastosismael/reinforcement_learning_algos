using Statistics
using StatsBase
using Distributions

# Monte Carlo Exploring Starts( Estimating π optimal)
struct State
    x::Int
    y::Int
end

# Definition of the environment
n = 4
𝒮 = [State(x,y) for x=1:n, y=1:n]
terminal_states = [State(2,2), State(3,3)]
inbounds(s::State) = 1 ≤ s.x ≤ n && 1 ≤ s.y ≤ n 
Base.:+(s1::State, s2::State) = State(s1.x + s2.x, s1.y +s2.y)
Base.:*(s1::State, b::Bool) = State(s1.x*b, s1.y*b)
isterminal(s) = s ∈ terminal_states

@enum Actions UP DOWN LEFT RIGHT STAY
const 𝒜 = Dict(UP => State(0, 1), DOWN => State(0, -1), RIGHT => State(1, 0), LEFT => State(-1, 0), STAY => State(0,0))

function R(s)
    if s ∈ terminal_states
        return 10
    else
        return -1
    end
end

function π₀(s)
    return rand(𝒜)[1]
end


# This code is hand-made just for the monte carlos ES (The worst algorithm ever)
function generate_episode_es(𝒮, 𝒜, π::Dict, πₒ::Function)
    s₀ = rand(𝒮)
    a₀ = rand(𝒜)[1]
    if ~inbounds(s₀ + 𝒜[a₀])
        a₀ = STAY
        episode = [(s₀, a₀)]
    end
    episode = [(s₀, a₀)]
    if isterminal(s₀)
        return episode
    end
    s = s₀ +  𝒜[a₀]
    while ~isterminal(s)
        a = get(π, s, π₀(s))
        if ~inbounds(s + 𝒜[a])
            a = STAY
        end
        push!(episode, (s, a))
        s = s + 𝒜[a]   
    end
    push!(episode, (s, STAY))
    return episode
end



# Problem: The montecarlo ES has the problem that
# if the π is not irreducible, so the function ends 
# in a eternal loop. 

function monte_carlo_es(𝒮, 𝒜, γ=1)
    # Initializaing π and Q
    πᵧ = Dict{State, Actions}()
    Q = Dict{State, Dict{Actions, Float64}}()
    returns = Dict((s,a) => [] for s in 𝒮, a ∈ keys(𝒜))
    for i in 1:5
        println(i)
        episode = generate_episode_es(𝒮, 𝒜, πᵧ, π₀)
        G = 0
        for t = reverse(1:length(episode)-1)
            println(t)
            sₜ, aₜ = episode[t][1], episode[t][2] 
            s₊ = episode[t+1][1]
            G = γ*G + R(s₊)
            if (sₜ, aₜ) ∉ episode[1:t-1]
                push!(returns[(sₜ, aₜ)], G)
                Q[sₜ] = Dict(aₜ => mean(returns[(sₜ, aₜ)]))
                πᵧ[sₜ] = findmax(Q[sₜ])[2]
            end
         end
    end
    return πᵧ
end

function generate_episode(𝒮, 𝒜, π::Dict)
    s₀ = rand(𝒮)
    a₀ = rand(𝒜)[1]
    if ~inbounds(s₀ + 𝒜[a₀])
        a₀ = STAY
        episode = [(s₀, a₀)]
    end
    episode = [(s₀, a₀)]
    if isterminal(s₀)
        return episode
    end
    s = s₀ +  𝒜[a₀]
    while ~isterminal(s)
        a = sample(collect(keys(π[s])), Weights(collect(values(π[s]))), 1)[1]
        if ~inbounds(s + 𝒜[a])
            a = STAY
        end
        push!(episode, (s, a))
        s = s + 𝒜[a]   
    end
    push!(episode, (s, STAY))
    return episode
end

function monte_carlo_control(ϵ, γ=1, n=100)
    # The comments are the case when we store the reutrns, less efficient
    π = Dict(s => Dict(a => ifelse(s ∉ terminal_states, 1/length(𝒜), ifelse(a == STAY, 1.0, 0.0)) for a ∈ keys(𝒜)) for s ∈ 𝒮)
    # Q = Dict{State, Dict{Actions, Float64}}()
    # returns = Dict((s,a) => [] for s ∈ 𝒮, a ∈ keys(𝒜))
    C = Dict((s,a) => 0 for s ∈ 𝒮, a ∈ keys(𝒜))
    Q = Dict{State, Dict{Actions, Float64}}(s => Dict(a => 0.0 for a ∈ keys(𝒜)) for s ∈ 𝒮)
    for i=1:n
        println(i)
        episode = generate_episode(𝒮, 𝒜, π)
        G = 0
        for t=reverse(1:length(episode)-1)
            sₜ, aₜ = episode[t][1], episode[t][2] 
            s₊ = episode[t+1][1]
            G = γ*G + R(s₊)
            if (sₜ, aₜ) ∉ episode[1:t-1]
                # push!(returns[(sₜ, aₜ)], G)
                # Q[sₜ]  = Dict(aₜ => mean(returns[(sₜ, aₜ)]))
                C[sₜ,aₜ] += 1
                Q[sₜ][aₜ] = Q[sₜ][aₜ] + 1/C[sₜ,aₜ] * (G - Q[sₜ][aₜ])
                A = findmax(Q[sₜ])[2]
                for a ∈ keys(𝒜)
                    π[sₜ][a] = ifelse(a == A, 1 - ϵ + ϵ/length(𝒜), ϵ/length(𝒜))
                end
            end
        end
    end
    return(π)
end

function off_policy_mc(γ=1, n=500)
    Q = Dict(s => Dict(a => rand(Uniform(0,10)) for a ∈ keys(𝒜)) for s ∈ 𝒮)
    C = Dict(s => Dict(a => 0.0 for a ∈ keys(𝒜)) for s ∈ 𝒮)
    π = Dict(s => findmax(Q[s])[2] for s ∈ 𝒮)

    for i=1:n
        b = Dict(s => Dict(a => ifelse(s ∉ terminal_states, 1/length(𝒜), ifelse(a==STAY, 1.0, 0.0)) for a ∈ keys(𝒜)) for s ∈ 𝒮)
        episode = generate_episode(𝒮, 𝒜, b)
        G = 0
        W = 1
        for t=reverse(1:length(episode) -1)
            sₜ, aₜ = episode[t][1], episode[t][2]
            s₊ = episode[t+1][1]
            G = γ*G + R(s₊)
            C[sₜ][aₜ] = C[sₜ][aₜ] + W
            Q[sₜ][aₜ] = Q[sₜ][aₜ] + W/C[sₜ][aₜ] * (G - Q[sₜ][aₜ])
            π[sₜ] = findmax(Q[sₜ])[2]
            if aₜ ≠ π[sₜ] 
                break
            end
            W = W * 1/b[sₜ][aₜ]
        end
    end
    return π
end


function sarsa(α, ϵ, γ=1)
    Q = Dict(s => Dict(a => -200.0 for a ∈ keys(𝒜)) for s ∈ 𝒮)
    π = Dict(s => Dict(a => ifelse(s ∉ terminal_states, 1/length(𝒜), ifelse(a == STAY, 1.0, 0.0)) for a ∈ keys(𝒜)) for s ∈ 𝒮)
    for i=1:10000
        episode = generate_episode(𝒮, 𝒜, π)
        for t=1:length(episode)-1
            s, a = episode[t][1], episode[t][2]
            s′, a′ = episode[t+1][1], episode[t+1][2]
            r = R(s′)
            Q[s][a] = Q[s][a] + α*(r + γ*Q[s′][a′] - Q[s][a])
            s = s′
            a = a′
        end
        π = Dict(s => Dict(a => ifelse(a == findmax(Q[s])[2], 1 - ϵ + ϵ/length(𝒜), ϵ/length(𝒜)) for a ∈ keys(𝒜)) for s ∈ 𝒮)
    end
    return π
end

function q_learning(α, ϵ, γ=1)
    Q = Dict(s => Dict(a => -100.0 for a ∈ keys(𝒜)) for s ∈ 𝒮)
    π = Dict(s => Dict(a => ifelse(s ∉ terminal_states, 1/length(𝒜), ifelse(a == STAY, 1.0, 0.0)) for a ∈ keys(𝒜)) for s ∈ 𝒮)
    for i=1:8000
        episode = generate_episode(𝒮, 𝒜, π)
        for t=1:length(episode)-1
            s, a = episode[t][1], episode[t][2]
            s′ = episode[t+1][1]
            r = R(s′)
            Q[s][a] = Q[s][a] + α*(r + γ*findmax(Q[s′])[1] - Q[s][a])
        end
        π = Dict(s => Dict(a => ifelse(a == findmax(Q[s])[2], 1 - ϵ + ϵ/length(𝒜), ϵ/length(𝒜)) for a ∈ keys(𝒜)) for s ∈ 𝒮)
    end
    return π
end


#π = Dict(s => rand(𝒜) for s in 𝒮)

# A = generate_episode(𝒮, 𝒜, Dict{State, Actions}(), π₀)
#mt = monte_carlo_es(𝒮, 𝒜)
#m = monte_carlo_control(1e-2, 1, 1000)
#ms = monte_carlo_es(𝒮, 𝒜)
# off_mc = off_policy_mc()
#print(m[State(1,2)])
s =  sarsa(0.5, 0.01)
d = q_learning(0.5, 0.01)

