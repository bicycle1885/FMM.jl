type AlignmentProfile
    seed_interval::Int
    seed_extend_length::Int
    max_seed_extension::Int
    score_model::AffineGapScoreModel{Int}
    score_model8::AffineGapScoreModel{Int8}

    function AlignmentProfile(;kwargs...)
        dict = Dict(kwargs)
        seed_interval = dict[:seed_interval]
        seed_extend_length = dict[:seed_extend_length]
        max_seed_extension = dict[:max_seed_extension]
        # create affine gap model
        submat = Array{Int}(5, 5)
        fill!(submat, dict[:mismatching_score])
        submat[diagind(submat)] = dict[:matching_score]
        score_model = AffineGapScoreModel(
            submat,
            gap_open_penalty=dict[:gap_open_penalty],
            gap_extend_penalty=dict[:gap_extend_penalty]
        )
        submat = Array{Int8}(5, 5)
        fill!(submat, dict[:mismatching_score])
        submat[diagind(submat)] = dict[:matching_score]
        score_model8 = AffineGapScoreModel(
            submat,
            gap_open_penalty=dict[:gap_open_penalty],
            gap_extend_penalty=dict[:gap_extend_penalty]
        )
        return new(
            seed_interval,
            seed_extend_length,
            max_seed_extension,
            score_model,
            score_model8
        )
    end
end
