using Parameters
using Statistics

struct State
	x::Int
	y::Int
end

@with_kw struct params
	N::Int = 4
	null_state::Array{State} = [State(1,1), State(4,4)]
	p::Real = 0.5
end

par = params()
𝒮 = [State(i,j) for i=1:par.N for j=1:par.N]
@enum Action LEFT RIGHT UP DOWN
𝒜 = [LEFT, RIGHT, UP, DOWN]

inbounds(s::State) = 1 ≤ s.x ≤ par.N && 1 ≤ s.y ≤ par.N

const MOVEMENTS = Dict(LEFT => State(-1,0),
						RIGHT => State(1,0),
						UP => State(0,1),
						DOWN => State(0,-1));
Base.:+(s1::State, s2::State) = State(s1.x + s2.x, s1.y + s2.y)

function T(s::State, a::Action, s′::State)
	return 1
end

function R(s::State)
	if s ∈ par.null_state
		return -1
	else
		return -1
	end
end

struct MDP
	γ::Real
	𝒮::Array{State}
	𝒜::Array{Action}
	T
	R
end

function π(a, s)
	if s ∈ par.null_state
		return 0
	else
		return 1/length(𝒜)
	end
end

γ = 1
𝒫 = MDP(γ, 𝒮, 𝒜, T, R)

function lookahead(𝒫::MDP, 𝒱::Dict, s)
	𝒮, 𝒜, T, R, γ = 𝒫.𝒮, 𝒫.𝒜, 𝒫.T, 𝒫.R, 𝒫.γ
	value = 0
	for a ∈ 𝒜
		if ~inbounds(s + MOVEMENTS[a])
			s′ = s
		else
			s′ = s + MOVEMENTS[a]
		end
		value += π(a, s) * T(s,a,s′) * (γ*𝒱[s′] + R(s′))
	end
	return value
end

function iterative_policy_evaluation(𝒫::MDP, π, k_max=true)
	𝒮, T, R, γ = 𝒫.𝒮, 𝒫.T, 𝒫.R, 𝒫.γ
	𝒱 = Dict(s => 0 for s in 𝒮)
	for k in 1:k_max
		𝒱 = Dict(ifelse(s ∉  par.null_state, s => lookahead(𝒫, 𝒱, s), s => 0) for s in 𝒮)
	end
	return 𝒱
end


function policy_iteration(𝒫::MDP, π, θ)
	𝒮, T, R, γ = 𝒫.𝒮, 𝒫.T, 𝒫.R, 𝒫.γ
	𝒱 = Dict(s => 0 for s in 𝒮)
	Δ = 1
	while Δ > θ 
		𝒱ₐ = 𝒱
		𝒱 = Dict(ifelse(s ∉  par.null_state, s => lookahead(𝒫, 𝒱, s), s => 0) for s in 𝒮)
		Δ = maximum([abs(x) for x ∈ collect(values(𝒱)) - collect(values(𝒱ₐ))])
	end
	return 𝒱
end

function return_valid_state(s, a)
	if ~inbounds(s + MOVEMENTS[a])
		s′ = s
	else
		s′ = s + MOVEMENTS[a]
	end
	return s′
end




function policy_improvement(𝒫::MDP, π)
	𝒮, T, R, γ = 𝒫.𝒮, 𝒫.T, 𝒫.R, 𝒫.γ
	π₀ = Dict(s => UP for s ∈ 𝒮)
	# Policy Evaluation
	𝒱 = policy_iteration(𝒫, π, 0.01)
	#Policy Improvement
	policy_stable = false
	π₁ = π₀
	while ~policy_stable
		for s ∈ 𝒮 
			action_values = Dict(a => T(s, a, return_valid_state(s,a))*(R(return_valid_state(s,a)) + γ*𝒱[return_valid_state(s,a)]) for a ∈ 𝒜)
			v, a₀ = findmax(action_values)
			π₁[s] = a₀
		end
		if π₀ == π₁
			policy_stable = true
		else
			π₀ = π₁
		end
	end
	return π₀
end


function gen_episode(𝒫::MDP)
	𝒜, 𝒮, R = 𝒫.𝒜, 𝒫.𝒮, 𝒫.R
	s₀ = rand(𝒮)
	trajectory = [s₀]
	s = s₀
	while s ∉ par.null_state
		s = rand(Set([ifelse(inbounds(s+MOVEMENTS[a]), s+MOVEMENTS[a], s) for a ∈ 𝒜]) )
		push!(trajectory, s)
	end
	return trajectory
end
function mc_policy_evaluation(π, k_max=100)
	𝒱 = Dict(s => 0.0 for s ∈ 𝒮)
	returns = Dict(s => [] for s ∈ 𝒮)
	for i in 1:k_max
		episode = gen_episode(𝒫)
		G = 0
		for i in reverse(1:length(episode))
			G = γ*G + R(episode[i])
			if episode[i] ∉ episode[1:i-1] && episode[i] ∉ par.null_state
				push!(returns[episode[i]], G)
				𝒱[episode[i]] = mean(returns[episode[i]])
			end
		end
	end
	return 𝒱
end
		




# 𝒱₀ = Dict(s => 0 for s in 𝒮)
# lookahead(𝒫, 𝒱₀, State(2,2))

gen_episode(𝒫)
policy_iteration(𝒫, π, 0.0001)
policy_improvement(𝒫, π)
mc_policy_evaluation(𝒫, 20000)