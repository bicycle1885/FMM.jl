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
    info("aligning reads")
    t = @elapsed for rec in reads
        setread!(readstate, rec.seq)
        align_read!(readstate, index, profile)
        println(rec.name)
        if isaligned(readstate)
            aln = alignment(readstate)
            chr, loc = locus(index.genome, aln[2].startpos)
            println(chr, ':', loc)
            println(aln)
        else
            println("not aligned")
        end
        println()
    end
    info("finished: ", t, " s")
end

end # module
