# score type (sync with `score_t` in alnscore.c)
typealias Score Int16

# score parameters (sync with `struct params_s` in alnscore.c)
type ScoreParams
    matching_score::Score
    mismatching_score::Score
    gap_open_penalty::Score
    gap_ext_penalty::Score
end

function ScoreParams(;
        matching_score=0,
        mismatching_score=0,
        gap_open_penalty=0,
        gap_ext_penalty=0,
    )
    return ScoreParams(matching_score, mismatching_score, gap_open_penalty, gap_ext_penalty)
end

# DNA sequence (sync with `struct dnaseq_s` in alnscore.c)
immutable DNASeq
    seq::Ptr{DNANucleotide}
    len::Csize_t
end

function convert(::Type{DNASeq}, seq::Vector{DNANucleotide})
    return DNASeq(pointer(seq), length(seq))
end

function convert(::Type{DNASeq}, seq::DNASequence)
    return convert(DNASeq, [x for x in seq])
end

# sync with `N_PAR` in alnscore.c
const N_PAR = 8

const libalnscore = Pkg.dir("FMM", "deps", "libalnscore.so")

# global buffer
const _buffer = ccall((:make_buffer, libalnscore), Ptr{Void}, ())

n_call = 0

function alignment_scores!(scores, params, query, seqs, n_seqs=length(seqs))
    @assert length(scores) == N_PAR
    @assert 0 ≤ n_seqs ≤ N_PAR
    ret = ccall((:alignment_score, libalnscore),
        Cint,
        (Ptr{Void}, Ptr{ScoreParams}, DNASeq, Ptr{DNASeq}, Cint, Ptr{Score}),
        _buffer, pointer_from_objref(params), query, seqs, n_seqs, scores
    )
    global n_call
    n_call += 1
    @assert ret == 0 "alignment error"
    return scores
end

function alignment_scores(params, query, seqs)
    return alignment_scores!(Vector{Score}(N_PAR), params, query, seqs)
end
