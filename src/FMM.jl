module FMM

import Base:
    getindex,
    setindex!,
    length,
    endof


using Bio.Seq
using Bio.Align
using IntArrays
using FMIndices
import Bio
import Bio.Align: Character

include("index.jl")
include("align.jl")

function run_alignment(index, read_file)
    format = endswith(read_file, ".fa") ? FASTA :
             endswith(read_file, ".fq") ? FASTQ :
             error("unknown format")
    reads = read(read_file, format)
    # same as --sensitive option of Bowtie2 if nt = 100
    profile = Profile(
        seed_length=22,
        seed_interval=ceil(Int, 1 + 1.15 * sqrt(100)),
        max_consective_fails=15
    )
    readstate = ReadState()
    t = @elapsed for rec in reads
        setread!(readstate, rec.seq)
        align_read!(readstate, index, profile)
        if isaligned(readstate)
            println(alignment(readstate))
        end
    end
    info(t, "s")
end

end # module
