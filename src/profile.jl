type AlignmentProfile
    seed_length::Int
    seed_interval::Int
    max_trials_per_seedhit::Int
    score_params::ScoreParams
    function AlignmentProfile(;kwargs...)
        dict = Dict(kwargs)
        seed_length = dict[:seed_length]
        seed_interval = dict[:seed_interval]
        max_trials_per_seedhit = dict[:max_trials_per_seedhit]
        score_params = ScoreParams(
            dict[:matching_score],
            dict[:mismatching_score],
            dict[:gap_open_penalty],
            dict[:gap_ext_penalty]
        )
        return new(
            seed_length,
            seed_interval,
            max_trials_per_seedhit,
            score_params
        )
    end
end
