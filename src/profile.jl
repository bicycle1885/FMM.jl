type AlignmentProfile
    seed_interval::Int
    max_seed_extension::Int
    max_seed_try::Int
    score_model::AffineGapScoreModel{Int}
    score_model16::AffineGapScoreModel{Int16}

    function AlignmentProfile(;kwargs...)
        dict = Dict(kwargs)
        seed_interval = dict[:seed_interval]
        max_seed_extension = dict[:max_seed_extension]
        max_seed_try = dict[:max_seed_try]
        # create affine gap model
        submat = Array{Int}(5, 5)
        fill!(submat, dict[:mismatching_score])
        submat[diagind(submat)] = dict[:matching_score]
        score_model = AffineGapScoreModel(
            submat,
            gap_open_penalty=dict[:gap_open_penalty],
            gap_extend_penalty=dict[:gap_extend_penalty]
        )
        submat = Array{Int16}(5, 5)
        fill!(submat, dict[:mismatching_score])
        submat[diagind(submat)] = dict[:matching_score]
        score_model16 = AffineGapScoreModel(
            submat,
            gap_open_penalty=dict[:gap_open_penalty],
            gap_extend_penalty=dict[:gap_extend_penalty]
        )
        return new(
            seed_interval,
            max_seed_extension,
            max_seed_try,
            score_model,
            score_model16
        )
    end
end
