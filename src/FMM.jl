module FMM

importall Base

using Bio.Seq
using IntArrays
using FMIndices
import Bio

include("index.jl")
include("alignscore.jl")
include("profile.jl")
include("readstate.jl")
include("align.jl")

function run_alignment(profile::AlignmentProfile, index, read_file)
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
