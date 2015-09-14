module FMM

importall Base

using Bio.Seq
#using Bio.Align
using IntArrays
using FMIndices
import Bio
#import Bio.Align: Character

include("index.jl")
include("alignscore.jl")
include("readstate.jl")
include("align.jl")

function run_alignment(index, read_file)
    # same as --sensitive option of Bowtie2 if nt = 100
    profile = Profile(
        # seed search
        seed_length=22,
        seed_interval=ceil(Int, 1 + 1.15 * sqrt(100)),
        # effort limits
        #max_consective_fails=15,
        max_trials_per_seedhit=2N_PAR,
        # alignment scores
        matching_score=0,
        mismatching_score=-6,
        gap_open_penalty=5,
        gap_ext_penalty=3,
    )
    format = endswith(read_file, ".fa") ? FASTA :
             endswith(read_file, ".fq") ? FASTQ :
             error("unknown format")
    reads = open(read_file, format)
    readstate = ReadState()
    t = @elapsed for rec in reads
        setread!(readstate, rec.seq)
        align_read!(readstate, index, profile)
        if isaligned(readstate)
            println(alignment(readstate))
        end
    end
    info(t, "s")
    global n_call
    @show n_call
end

end # module
