#!/bin/env julia
using Random
using Combinatorics

import Base: isless

# SACO contants based on
# the results of the paper
m = 40
ρ = .04
β = 2
ω = .1
λ = 1
τmin = .001
stp = .001
τmax = .999
Nmax = 500
nants = 20

mutable struct VisibleArc
    satellite       # sa
    antenna         # ea
    start_time      # ts
    end_time        # te
    heuristic       # ηi
end

struct Satellite    # sa
    start_time
    end_time
end

struct Antenna      # ea
    start_time
    end_time
end

mutable struct Ant
    pheromone_value
    solution
    fitness
end

isless(s1::Satellite, s2::Satellite) = isless(s1.start_time, s2.start_time)
isless(a1::Antenna, a2::Antenna) = isless(a1.start_time, a2.start_time)


function sort!(va::VisibleArc)
    i = sortperm(va.satellite)
    va.satellite = va.satellite[i]
    i = sortperm(va.antenna)
    va.antenna = va.antenna[i]
end

function canConnect(a::Antenna, s::Satellite)
    return a.start_time <= s.end_time && s.start_time <= a.end_time
end

function conflict(s1, s2)
    return s2.start_time <= s1.end_time && s1.start_time <= s2.end_time
end

function feasible(c, s, e)
    return c.start_time >= s && c.start_time <= e
end

function generateRandomDataset()
    # Create our visible arc
    va = VisibleArc(Vector{Satellite}(), Vector{Antenna}(), τmin, τmax, [])
    # Random number of satellites in [10,20]
    nsatellites = rand(10:1:20)
    # Random number of antennas in [5,12]
    nantennas = rand(5:1:12)
    while nantennas >= nsatellites
        nantennas = rand(5:1:12)
    end
    # Add random satellites throughout the visible arc
    for i in 1:nsatellites
        s = rand(τmin:stp:τmax-stp)
        e = rand(τmin:stp:τmax)
        while e <= s
            e = rand(τmin:stp:τmax)
        end
        push!(va.satellite, Satellite(s, e))
    end
    # Add random antennas throughout the visible arc
    for i in 1:nantennas
        s = rand(τmin:stp:τmax-stp)
        e = rand(τmin:stp:τmax)
        while e <= s
            e = rand(τmin:stp:τmax)
        end
        push!(va.antenna, Antenna(s, e))
    end

    sort!(va)

    return va
end

# Generate all antenna-satellite pairs first
# (No triplets or more for now)
# function generateFeasiblePairs!(result, antennas, satellites, s, e)
#     for antenna in antennas
#         if feasible(antenna, s, e)
#             for sat in satellites
#                 if canConnect(antenna, sat)
#                     push!(result, VisibleArc([sat], [antenna], max(antenna.start_time, sat.start_time), min(antenna.end_time, sat.end_time), 0))
#                 end
#             end
#         end
#     end
# end

# TODO
# function checkForMoreCombinations!(result, antennas, satellites, s, e)
#     found = false
#     for sat in satellites
#         for va in result
#             for antenna in va.antenna
#                 if canConnect(antenna, sat)
#                     conf = false
#                     for s in va.satellite
#                         if conflict(sat,s)
#                             conf = true
#                             break
#                         end
#                         if !conf
#                             found = true
#                             push!(va.satellite, sat)
#                         end
#                     end
#                 end
#             end
#         end
#     end
#     return found
# end

# A not-optimized-at-all combinations search
# function combinations(antennas, satellites, working_lengths)
#     result = Vector{VisibleArc}()
#     end_times = working_lengths .+ (τmax/length(working_lengths) - stp)
#     for t=1:length(working_lengths)
#         s = working_lengths[t]
#         e = end_times[t]
#         generateFeasiblePairs!(result, antennas, satellites, s, e)
#         work = true
#         while work
#             work = checkForMoreCombinations!(result, antennas, satellites, s, e)
#         end
#     end
#     return result
# end

function validCombination(c, s, e)
    for a in c.antenna
        # Check if this antenna is
        # outside of this working period
        if !feasible(a, s, e)
            return false
        end
        for sat in c.satellite
            # Check if this satellite is
            # outside of this working period
            if !feasible(sat, s, e)
                return false
            end
            # Check if this antenna can connect
            # to this satellite
            # if !canConnect(a, sat)
            #     return false
            # end
            for sat2 in c.satellite
                # Check if this satellite is in
                # conflict with another satellite
                if sat != sat2 && conflict(sat, sat2)
                    return false
                end
            end
        end
        # TODO not sure how to handle this kind of conflicts
        # Same for sats
        for a2 in c.antenna
            # Check if this antenna is in
            # conflict with another antenna
            if a != a2 && conflict(a, a2)
                return false
            end
        end
    end
    return true
end

function fixCombination!(v, working_lengths)
    mins = τmax
    maxe = τmin
    for a in v.antenna
        if a.start_time < mins
            mins = a.start_time
        end
        if a.end_time < maxe
            maxe = a.end_time
        end
    end
    for s in v.satellite
        if s.start_time < mins
            mins = s.start_time
        end
        if s.end_time < maxe
            maxe = s.end_time
        end
    end
    v.start_time = mins
    v.end_time = maxe
    η!(v, working_lengths)
end

function legalCombinations(antennas, satellites, working_lengths)
    lc = Vector{VisibleArc}()
    end_times = working_lengths .+ (τmax/length(working_lengths) - stp)
    antenna_combs = collect(combinations(antennas))
    satellite_combs = collect(combinations(satellites))
    total = length(working_lengths) * length(antenna_combs) * length(satellite_combs)
    progress = 0
    inc = length(antenna_combs) * length(satellite_combs)
    for t=1:length(working_lengths)
        println("Progress = " * string(progress/total*100) * "%")
        for a in antenna_combs
            for s in satellite_combs
                tmpv = VisibleArc(s, a, 0, 0, [])
                if validCombination(tmpv, working_lengths[t], end_times[t])
                    fixCombination!(tmpv, working_lengths)
                    push!(lc, tmpv)
                end
            end
        end
        progress += inc
    end

    return lc
end

# function candicatesPerWorkingLength(visible_arc, working_lengths)
#     candidates = Vector{VisibleArc}()
#     end_times = working_lengths .+ (τmax/length(working_lengths) - stp)
#     for t=1:length(working_lengths)
#         current_wl = Vector{VisibleArc}()

#         # Associate antennas with all feasible sat combinations
#         for a in visible_arc.antenna
#             va = VisibleArc([], [], 0, 0, [])
#             if feasible(a, working_lengths[t], end_times[t])
#                 push!(va.antenna, a)
#                 for i=1:length(visible_arc.satellite)
#                     if canConnect(a, visible_arc.satellite[i])
#                         push!(va.antenna, visible_arc.satellite[i])
#                     else

#                     end
#                 end
#             else # Antennas are sorted by start_time
#                 break
#             end
#         end
#         # for i = 1:length(visible_arc.satellite)
#         #     if visible_arc.satellite[i].start_time >= working_lengths[t]
#         #         push!(va.satellite, visible_arc.satellite[i])
#         #     else # Satellites are sorted
#         #         break
#         #     end
#         # end
#         push!(candidates, va)
#     end
#     return candidates
# end

function initAnts(nants, visible_arcs)
    ants = Vector{Ant}()
    for n=1:nants
        tmpant = Ant(ω*τmax, Vector{VisibleArc}(), 0)
        for va in visible_arcs
            tmpva = VisibleArc(Vector{Satellite}(), Vector{Antenna}(), 0, 0, [])
            # How many (and which) satellites are we
            # going to use for this working period?
            sat_idx = sort(randperm(rand(0:1:length(va.satellite))))
            sats = va.satellite[sat_idx]
            append!(tmpva.satellite, sats)
            # How many (and which) antennas are we
            # going to use for this working period?
            antns_idx = sort(randperm(rand(0:1:length(va.antenna))))
            antns = va.antenna[antns_idx]
            append!(tmpva.antenna, antns)

            push!(tmpant.solution, tmpva)
        end
        push!(ants, tmpant)
    end
    return ants
end

function η!(visible_arc, working_lengths)
    for l=1:length(working_lengths)
        tsi = visible_arc.start_time
        tsl = working_lengths[l]
        tei = visible_arc.end_time
        tel = l == length(working_lengths) ? τmax : working_lengths[l+1]-stp
        Δti1 = tsl > tsi ? tsl-tsi : 0
        Δti2 = tel < tei ? tei-tel : 0
        Δti = Δti1 + Δti2
        push!(visible_arc.heuristic, ℯ^(-λ*Δti))
    end
end

function update!(ants)

end

function score!(ants, working_lengths)
    for ant in ants

    end
end

function saco(candicates, working_lengths)
    # Paper step 1
    ants = initAnts(nants, candidates)
    update!(ants)
    score!(ants, working_lengths)
    iteration = 0
    # Paper step 2
    while iteration <= Nmax
        # Paper step 3
        for i = 1:length(ants)
            ants[i].solution = constructSolution(ants[i])
        end
        iteration += 1
    end
end

function main()
    l = 20 # working periods
    working_lengths = Array(τmin:τmax/l:τmax)
    visible_arc = generateRandomDataset()
    candidates = legalCombinations(visible_arc.antenna, visible_arc.satellite, working_lengths)
    println(candidates)
    println(length(candidates))
    # TODO
    # Visualize the visible arc here
    # solution = saco(candidates, working_lengths)
    # TODO
    # Visualize the solution here
end

main()
