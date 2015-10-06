module FMM

importall Base

using Bio.Seq
using IntArrays
using FMIndexes
using PairwiseAlignment
using IndexableBitVectors
import Bio

include("nmask.jl")
include("genomicseq.jl")
include("genome.jl")
include("index.jl")
include("small_pque.jl")
include("seediter.jl")
include("alignscore.jl")
include("profile.jl")
include("seedhit.jl")
include("readstate.jl")
include("align.jl")

end # module
