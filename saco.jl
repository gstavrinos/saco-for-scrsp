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

function generateRandomDataset()
    # Create our visible arc
    va = VisibleArc(Vector{Satellite}(), Vector{Antenna}(), τmin, τmax, [])
    # Random number of satellites in [10,17]
    nsatellites = rand(10:1:17)
    # Random number of antennas in [5,13]
    nantennas = rand(5:1:13)
    while nantennas >= nsatellites
        nantennas = rand(5:1:13)
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

function validCombination(c)
    for a in c.antenna
        # Check if this antenna is
        # outside of this working period
        for sat in c.satellite
            if !canConnect(a, sat)
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

function feasibleCombinations(antennas, satellites, working_lengths)
    lc = Vector{VisibleArc}()
    antenna_combs = collect(combinations(antennas))
    satellite_combs = collect(combinations(satellites))
    total = length(antenna_combs) * length(satellite_combs)
    progress = 0
    prints_left = ceil(0.05*length(antenna_combs))
    for a in antenna_combs
        for s in satellite_combs
            tmpv = VisibleArc(s, a, 0, 0, [])
            if validCombination(tmpv)
                fixCombination!(tmpv, working_lengths)
                push!(lc, tmpv)
            end
        end
        progress += length(satellite_combs)
        prints_left -= 1
        if prints_left % 10 == 0
            println("Progress = " * string(progress/total*100) * "%")
        end
    end

    return lc
end

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

function saco(candidates, working_lengths)
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
    println("Generating random dataset...")
    visible_arc = generateRandomDataset()
    println("Number of antennas = " * string(length(visible_arc.antenna)))
    println("Number of satellites = " * string(length(visible_arc.satellite)))
    println("Generating candidate combinations...")
    candidates = feasibleCombinations(visible_arc.antenna, visible_arc.satellite, working_lengths)
    println("Generated " * string(length(candidates)) * " candidate combinations.")
    println("Now starting the SACO algorithm...")
    # TODO
    # Visualize the visible arc here
    # solution = saco(candidates, working_lengths)
    # TODO
    # Visualize the solution here
end

@time main()
