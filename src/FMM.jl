module FMM

importall Base

using Bio.Seq
using Bio.Align
using SIMDAlignment
using IntArrays
using FMIndexes
using IndexableBitVectors
using IntervalHeaps
import Bio

include("nmask.jl")
include("genomicseq.jl")
include("genome.jl")
include("index.jl")
include("seediter.jl")
include("alignscore.jl")
include("profile.jl")
include("seedhit.jl")
include("readstate.jl")
include("align.jl")

end # module
