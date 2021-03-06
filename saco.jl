#!/bin/env julia
# Algorithm imports
using Random
using StatsBase
using ProgressMeter
using Combinatorics
import Base: isless

# Visualization imports
using Makie
using FileIO
using MeshIO
using Colors
using GeometryTypes


# Algorithm code
# ------------------

# SACO contants based on
# the results of the paper
nants = 20
ρ = .04
β = 2
ω = .1
λ = 1
τmin = .001
stp = .001
τmax = .999
Nmax = 500

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
    index_pool
end

# Store all candidates here
# because my PC can't handle
# adding them to each ant
candidates = Vector{VisibleArc}()

isless(s1::Satellite, s2::Satellite) = isless(s1.start_time, s2.start_time)
isless(a1::Antenna, a2::Antenna) = isless(a1.start_time, a2.start_time)


function sort!(va::VisibleArc)
    i = sortperm(va.satellite)
    va.satellite = va.satellite[i]
    i = sortperm(va.antenna)
    va.antenna = va.antenna[i]
end

function canConnect(a::Antenna, s::Satellite)
    return (s.start_time >= a.start_time && s.start_time <= a.end_time) ||
            (s.end_time >= a.start_time && s.end_time <= a.end_time) ||
            (a.start_time <= s.start_time && a.end_time >= s.end_time) ||
            (a.start_time >= s.start_time && a.end_time <= s.end_time)
end

function feasible(c, s, e)
    return c.start_time < e && c.end_time > s
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
    if !isempty(ARGS) && ARGS[1] == "--quick-test"
        nsatellites = 5
        nantennas = 4
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
        if a.end_time > maxe
            maxe = a.end_time
        end
    end
    for s in v.satellite
        if s.start_time < mins
            mins = s.start_time
        end
        if s.end_time > maxe
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
    @showprogress 0.1 "Generating candidate combinations..." for a in antenna_combs
        for s in satellite_combs
            tmpv = VisibleArc(s, a, 0, 0, [])
            if validCombination(tmpv)
                fixCombination!(tmpv, working_lengths)
                push!(lc, tmpv)
            end
        end
    end

    return lc
end

function initAnts(nants, l)
    ants = [Ant(ω*τmax, Vector{VisibleArc}(), 0, collect(1:length(candidates))) for _=1:nants]
    for ant in ants
        randind = sample(ant.index_pool, l, replace=false)
        append!(ant.solution, candidates[randind])
        f!(ant)
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

function f(solution::Vector{VisibleArc})
    s = 0
    for i=1:length(solution)
        s += solution[i].heuristic[i]
    end
    return s
end

function f!(ant::Ant)
    ant.fitness = f(ant.solution)
    return ant.fitness
end

function update!(ants, curr_best)
    for ant in ants
        ant.pheromone_value *= (1-ρ)
        if ant.solution == curr_best
            ant.pheromone_value += ρ*τmax
        end
        ant.pheromone_value = min(max(ant.pheromone_value, τmin), τmax)
    end
end

function constructSolution!(ant, working_lengths, ants)
    end_times = working_lengths .+ (τmax/length(working_lengths) - stp)
    ant.solution = Vector{VisibleArc}()
    for l=1:length(working_lengths)
        probs = Vector{Float64}()
        for ai in ant.index_pool
            if feasible(candidates[ai], working_lengths[l], end_times[l])
                pd = 0
                for a in ants
                    if a != ant
                        for tmpai in a.index_pool
                            if feasible(candidates[tmpai], working_lengths[l], end_times[l])
                                pd += a.pheromone_value*a.fitness^β
                            else
                                break
                            end
                        end
                    end
                end
                push!(probs, ant.pheromone_value*ant.fitness^β/pd)
            else
                push!(probs, .0)
            end
        end
        if sum(probs) == 0
            push!(ant.solution, candidates[sample(ant.index_pool)])
        else
            push!(ant.solution, candidates[sample(ant.index_pool, Weights(probs))])
        end
        delete_indices = []
        # for a in ant.solution[end].antenna
        #     for i=1:length(ant.index_pool)
        #         for a_p in candidates[ant.index_pool[i]].antenna
        #             if a == a_p
        #                 push!(delete_indices, i)
        #                 break
        #             end
        #         end
        #     end
        # end
        #deleteat!(ant.pool, findall(in(ant.solution[end].antenna), ant.pool[:].antenna))
        #deleteat!(ant.index_pool, delete_indices)
        #println(length(ant.index_pool))
    end
end

function saco(candidates, working_lengths)
    # Paper step 1
    ants = initAnts(nants, length(working_lengths))
    solution = Vector{VisibleArc}()
    update!(ants, solution)
    # Paper step 2
    @showprogress 0.1  "Running the SACO algorithm..." for iteration=1:Nmax
        # Paper step 3
        for i = 1:length(ants)
            # Paper step 4
            constructSolution!(ants[i], working_lengths, ants)
            # Paper step 7
            if isempty(solution) || f(solution) > f!(ants[i])
                # Paper step 8
                solution = ants[i].solution
            end
        end
        # Paper step 10
        update!(ants, solution)
        iteration += 1
    end
    return solution, f(solution)
end

# ------------------


# Visualization code
# ------------------

# Helper function to add the antennas on the globe
function placeAntennas!(scene, antennas, a)
    antenna_pos = getAntennaPositions()
    satellite_pos = getSatellitePositions()
    ΔΤmax = τmax-τmin
    nantennas = length(antennas)
    antp = Vector{Point{3,Float32}}()
    ant_colours = [RGBAf0(rand(), rand(), rand(), 1.0) for i = 1:nantennas]
    arrow_colours = []
    for c in ant_colours
        push!(arrow_colours, c)
        push!(arrow_colours, c)
    end
    directions = Vector{Point3f0}()
    start_pos = Vector{Point3f0}()
    for antenna in antennas
        # Go from [τmin,τmax] to [1,360]
        # which is the range of the earth mesh slices
        p = ceil(Int, 359/ΔΤmax*((antenna.end_time+antenna.start_time)/2-τmin)+1)
        push!(antp, antenna_pos[p])
        starti = ceil(Int, 359/ΔΤmax*(antenna.start_time-τmin)+1)
        endi = ceil(Int, 359/ΔΤmax*(antenna.end_time-τmin)+1)
        push!(directions, satellite_pos[starti]-antenna_pos[p])
        push!(directions, satellite_pos[endi]-antenna_pos[p])
        push!(start_pos, antenna_pos[p])
        push!(start_pos, antenna_pos[p])
    end
    arrows!(scene, start_pos, directions, arrowsize=0, linecolor=arrow_colours)
    ant_sizes = [(0.01,0.01,0.01) for i = 1:nantennas]
    ant_rot = [Vec4f0([0,0,-0.7,0.7]) for i = 1:nantennas]
    scene = meshscatter!(antp, color = ant_colours, markersize = ant_sizes, marker=a, rotation=ant_rot)
end

# Helper function to add the satellites over the globe
function placeSatellites!(scene, satellites, s)
    satellite_pos = getSatellitePositions()
    ΔΤmax = τmax-τmin
    nsats = length(satellites)
    satp = Vector{Point{3,Float32}}()
    sat_colours = [RGBAf0(rand(), rand(), rand(), 1.0) for i = 1:nsats]
    arrow_colours = []
    for c in sat_colours
        push!(arrow_colours, c)
        push!(arrow_colours, c)
    end
    directions = Vector{Point3f0}()
    start_pos = Vector{Point3f0}()
    for satellite in satellites
        p = ceil(Int, 359/ΔΤmax*((satellite.end_time+satellite.start_time)/2-τmin)+1)
        push!(satp, satellite_pos[p])
        starti = ceil(Int, 359/ΔΤmax*(satellite.start_time-τmin)+1)
        endi = ceil(Int, 359/ΔΤmax*(satellite.end_time-τmin)+1)
        push!(directions, satellite_pos[starti]-satellite_pos[p])
        push!(directions, satellite_pos[endi]-satellite_pos[p])
        push!(start_pos, satellite_pos[p])
        push!(start_pos, satellite_pos[p])
    end
    arrows!(scene, start_pos, directions, arrowsize=0, linecolor=arrow_colours)
    sat_sizes = [(0.2,0.2,0.2) for i = 1:nsats]
    scene = meshscatter!(satp, color = sat_colours, markersize = sat_sizes, marker=s)
end

# Helper function to generate possible antenna positions
function getAntennaPositions()
    return decompose(Point3f0, GLNormalUVMesh(Sphere(Point3f0(0), 1.05f0), 360))
end

# Helper function to generate possible antenna positions
function getSatellitePositions()
    return decompose(Point3f0, GLNormalUVMesh(Sphere(Point3f0(0), 1.5f0), 360))
end

# Helper function to load 3D models
function getModels()
    sat = load(string(@__DIR__)*"/res/satellite.obj", GLNormalUVMesh)
    antenna = load(string(@__DIR__)*"/res/antenna.obj", GLNormalUVMesh)
    return sat, antenna
end

# Helper function to load a scene with Earth inside it
function blueMarble♥()
    scene = Scene(resolution = (500, 500), backgroundcolor = :black, center=false)
    earth = load(string(@__DIR__)*"/res/bluemarble-2048.png")
    m = GLNormalUVMesh(Sphere(Point3f0(0), 1f0), 60)
    mesh!(scene, m, color = earth, shading = true, show_axis = false)
    return scene
end

# Helper function to add stars to the scene
function myGodItsFullOfStars!(scene)
    stars = 100_000
    scatter!(scene, (rand(Point3f0, stars) .- 0.5) .* 10,
        glowwidth = 0.005, glow_color = :white, color = RGBA(0.8, 0.9, 0.95, 0.4),
        markersize = rand(range(0.0001, stop=0.005, length=100), stars))
end

# Function to visualize the generated problem
function viz(v::VisibleArc)
    scene = blueMarble♥()
    myGodItsFullOfStars!(scene)
    sat, ant = getModels()
    placeAntennas!(scene, v.antenna, ant)
    placeSatellites!(scene, v.satellite, sat)
    display(scene)
    return scene
end

# Function to visualize the generated solution
function viz(scene::Scene, s::Vector{VisibleArc})
    antenna_pos = getAntennaPositions()
    satellite_pos = getSatellitePositions()
    ΔΤmax = τmax-τmin
    directions = Vector{Point3f0}()
    start_pos = Vector{Point3f0}()
    for va in s
        for antenna in va.antenna
            # Go from [τmin,τmax] to [1,360]
            # which is the range of the earth mesh slices
            anti = ceil(Int, 359/ΔΤmax*((antenna.end_time+antenna.start_time)/2-τmin)+1)
            for satellite in va.satellite
                sati = ceil(Int, 359/ΔΤmax*((satellite.end_time+satellite.start_time)/2-τmin)+1)
                push!(start_pos, antenna_pos[anti])
                push!(directions, satellite_pos[sati]-antenna_pos[anti])
            end
        end
    end
    arrows!(scene, start_pos, directions, arrowsize=0, linecolor=:green)
end

# ------------------

function main()
    l = 20 # working periods
    if !isempty(ARGS) && ARGS[1] == "--quick-test"
        l = 2
    end
    working_lengths = Array(τmin:τmax/l:τmax)
    println("Generating random dataset...")
    visible_arc = generateRandomDataset()
    println("Number of antennas = " * string(length(visible_arc.antenna)))
    println("Number of satellites = " * string(length(visible_arc.satellite)))
    scene = viz(visible_arc)
    println("Press enter to continue...")
    readline()
    global candidates = feasibleCombinations(visible_arc.antenna, visible_arc.satellite, working_lengths)
    println("Generated " * string(length(candidates)) * " candidate combinations.")
    solution, fitness = saco(candidates, working_lengths)
    #println(solution)
    viz(scene, solution)
    println("Press enter to quit...")
    readline()
end

@time main()
