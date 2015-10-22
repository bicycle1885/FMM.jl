type AlignmentProfile
    seed_length::Int
    seed_interval::Int
    max_trials_per_seedhit::Int
    max_seedcut_multiplier::Int
    score_model::AffineGapScoreModel{Score}

    function AlignmentProfile(;kwargs...)
        dict = Dict(kwargs)
        seed_length = dict[:seed_length]
        seed_interval = dict[:seed_interval]
        max_trials_per_seedhit = dict[:max_trials_per_seedhit]
        max_seedcut_multiplier = dict[:max_seedcut_multiplier]
        submat = Array{Score}(5, 5)
        fill!(submat, dict[:mismatching_score])
        submat[diagind(submat)] = dict[:matching_score]
        score_model = AffineGapScoreModel(
            submat,
            gap_open_penalty=dict[:gap_open_penalty],
            gap_extend_penalty=dict[:gap_extend_penalty]
        )
        return new(
            seed_length,
            seed_interval,
            max_trials_per_seedhit,
            max_seedcut_multiplier,
            score_model
        )
    end
end
