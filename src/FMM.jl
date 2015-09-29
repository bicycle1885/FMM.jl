module FMM

importall Base

using Bio.Seq
using IntArrays
using FMIndexes
using PairwiseAlignment
import Bio

include("index.jl")
include("small_pque.jl")
include("seediter.jl")
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
            #score, hit = readstate.best[1]
            #chr, loc = locus(index.genome, hit.location)
            #println(rec.name, '\t', chr, '\t', loc, '\t', score)
        end
    end
    info(t, "s")
end

end # module
