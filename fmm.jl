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
    # same as --sensitive option of Bowtie2 if nt = 100
    profile = FMM.AlignmentProfile(
        # seed search
        seed_length=22,
        seed_interval=ceil(Int, 1 + 1.15 * sqrt(100)),
        # effort limits
        #max_consective_fails=15,
        max_trials_per_seedhit=4 * FMM.N_PAR,
        # alignment scores
        matching_score=0,
        mismatching_score=-6,
        gap_open_penalty=5,
        gap_ext_penalty=3,
    )
    index = open(deserialize, string(genome_file, ".index"))
    FMM.run_alignment(profile, index, read_file)
    FMM.run_alignment(profile, index, read_file)
    @profile FMM.run_alignment(profile, index, read_file)
    Profile.print(STDERR, format=:tree, cols=1000)
end

main()
