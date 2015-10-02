module FMM

importall Base

using Bio.Seq
using IntArrays
using FMIndexes
using PairwiseAlignment
import Bio

include("genome.jl")
include("index.jl")
include("small_pque.jl")
include("seediter.jl")
include("alignscore.jl")
include("profile.jl")
include("readstate.jl")
include("align.jl")

end # module
