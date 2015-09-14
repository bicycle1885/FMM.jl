using FMM
using Bio.Seq
using FactCheck

facts("alignment score") do
    params = FMM.ScoreParams(
        matching_score=0,
        mismatching_score=-6,
        gap_open_penalty=5,
        gap_ext_penalty=3
    )
    query::FMM.DNASeq = dna"ACGT"

    # same length
    seqs = FMM.DNASeq[
        dna"ACGT",
        dna"ACGA",
        dna"ACAA",
        dna"AAAA",
        dna"ACGT",
        dna"ACGA",
        dna"ACAA",
        dna"AAAA",
    ]
    scores = FMM.alignment_scores(params, query, seqs)
    @fact scores --> [0, -6, -11, -14, 0, -6, -11, -14]

    # longer
    seqs = FMM.DNASeq[
        dna"ACGTA",
        dna"ACGAT",
        dna"ACAGT",
        dna"AACGT",
        dna"ACGTA",
        dna"ACGAT",
        dna"ACAGT",
        dna"AACGT",
    ]
    scores = FMM.alignment_scores(params, query, seqs)
    @fact scores --> [0, -6, -8, -8, 0, -6, -8, -8]
end

