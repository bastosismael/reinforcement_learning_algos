using Statistics
# Monte Carlo Exploring Starts( Estimating π optimal)
struct State
    x::Int
    y::Int
end

# Definition of the environment
n = 4
𝒮 = [State(x,y) for x=1:n, y=1:n]
terminal_states = [State(1,4), State(1,1)]
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


function generate_episode(𝒮, 𝒜, π::Dict, πₒ::Function)
    s₀ = rand(𝒮)
    a₀ = rand(𝒜)[1]
    episode = [(s₀, a₀)]
    s = s₀ +  (𝒜[a₀] *  inbounds(s₀ + 𝒜[a₀]))
    while ~isterminal(s)
        a = get(π, s, π₀(s))
        if ~inbounds(s + 𝒜[a])
            a = STAY
        end
        s = s + 𝒜[a]
        push!(episode, (s, a))
    end
    return episode
end

function monte_carlo_es(𝒮, 𝒜, γ=1)
    # Initializaing π and Q
    πᵧ = Dict{State, Actions}()
    Q = Dict(s => Dict(a => 0.0 for a ∈ keys(𝒜)) for s in 𝒮)
    returns = Dict((s,a) => 0.0 for s in 𝒮, a ∈ keys(𝒜))
    for i in 1:3
        episode = generate_episode(𝒮, 𝒜, πᵧ, π₀)
        G = 0
        for t = reverse(1:length(episode))
            print(t)
            sₜ, aₜ = episode[t][1], episode[t][2] 
            G = γ*G + R(sₜ)
            if (sₜ, aₜ) ∉ episode[1:t-1]
                returns[(sₜ, aₜ)] = G
                   Q[sₜ][aₜ] = mean(returns[(sₜ, aₜ)])
                 πᵧ[sₜ] = findmax(Q[sₜ])[2]
            end
         end
    end
    return πᵧ
end


#π = Dict(s => rand(𝒜) for s in 𝒮)

# A = generate_episode(𝒮, 𝒜, Dict{State, Actions}(), π₀)
mt = monte_carlo_es(𝒮, 𝒜)