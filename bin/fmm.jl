using DocOpt
import FMM

const doc = """
Usage:
    fmm.jl index <genome.fa>
    fmm.jl align <genome.fa> <reads.{fa,fq}>
"""

function main()
    srand(1234)
    args = docopt(doc)
    genome_file = args["<genome.fa>"]
    if args["index"]
        run_index(genome_file)
    elseif args["align"]
        read_file = args["<reads.{fa,fq}>"]
        run_align(genome_file, read_file)
    else
        @assert false
    end
end

function run_index(genome_file)
    index = FMM.build_index(genome_file)
    open(string(genome_file, ".index"), "w+") do io
        serialize(io, index)
    end
end

function run_align(genome_file, read_file)
    profile = FMM.AlignmentProfile(
        seed_interval=13,
        max_seed_hit=8,
        max_seed_try=8,
        # alignment scores
        matching_score=0,
        mismatching_score=-6,
        gap_open_penalty=5,
        gap_extend_penalty=3,
    )
    info("loading index")
    t = @elapsed index = open(deserialize, string(genome_file, ".index"))
    info("finished: ", t, " s")
    FMM.run_alignment(profile, index, read_file, STDOUT)
end

main()
