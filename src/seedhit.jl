#          |<-- seed hit -->|
# ~~~~~~~~~^^^^^^^^^^^^^^^^^^~~~~~~~~~~~
#          ^                ^
#          seed_start       seed_stop
immutable SeedHit
    # read range
    rrange::UnitRange{Int}
    # suffix array range
    sarange::UnitRange{Int}
    # wheather forward or reverse complement
    forward::Bool
end

seed_start(seedhit::SeedHit) = start(seedhit.rrange)
seed_stop(seedhit::SeedHit) = last(seedhit.rrange)
seed_length(seedhit::SeedHit) = length(seedhit.rrange)
count(seedhit::SeedHit) = length(seedhit.sarange)
getindex(seedhit::SeedHit, i::Integer) = seedhit.sarange[i]
isforward(seedhit::SeedHit) = seedhit.forward


#              |<-left->|<-- seed hit -->|<-right->|
# read:        ~~~~~~~~~^^^^^^^^^^^^^^^^^^~~~~~~~~~~
# genome: --------------------------------------------------
#                       ^
#                       location
immutable SeedHitExt
    seedhit::SeedHit
    # genome location
    location::Int
    # alignment scores
    lscore::Score  # left
    hscore::Score  # seedhit
    rscore::Score  # right
end

total_score(x::SeedHitExt) = x.lscore + x.hscore + x.rscore
isless(x::SeedHitExt, y::SeedHitExt) = isless(total_score(x), total_score(y))

# delegete
seed_start(seedhit::SeedHitExt) = seed_start(seedhit.seedhit)
seed_stop(seedhit::SeedHitExt) = seed_stop(seedhit.seedhit)
seed_length(seedhit::SeedHitExt) = seed_length(seedhit.seedhit)
isforward(seedhit::SeedHitExt) = isforward(seedhit.seedhit)
