using FMM
using Bio.Seq
using Bio.Align
using SIMDAlignment
using FMIndexes

if VERSION >= v"0.5-"
    using Base.Test
else
    using BaseTestNext
    const Test = BaseTestNext
end

@testset "data structures" begin
    @testset "NMask" begin
        bv = BitVector([true, false, true, false, false])
        nmask = FMM.NMask(bv)
        @test length(nmask) === 5
        for i in 1:endof(bv)
            @test nmask[i] === bv[i]
        end

        bv = falses(10000)
        nmask = FMM.NMask(bv)
        @test length(nmask) === 10000
        for i in 1:endof(bv)
            @test nmask[i] === bv[i]
        end
        @test findnext(nmask,    1) === 0
        @test findnext(nmask, 1000) === 0
        @test findprev(nmask,    1) === 0
        @test findprev(nmask, 1000) === 0

        bv[4] = true
        bv[601] = true
        bv[9900:10000] = true
        nmask = FMM.NMask(bv)
        @test findnext(nmask, 1)     === 4
        @test findnext(nmask, 5)     === 601
        @test findnext(nmask, 602)   === 9900
        @test findnext(nmask, 9901)  === 9901
        @test findnext(nmask, 10000) === 10000
        @test findprev(nmask, 10000) === 10000
        @test findprev(nmask, 9900)  === 9900
        @test findprev(nmask, 9899)  === 601
        @test findprev(nmask, 600)   === 4
        @test findprev(nmask, 3)     === 0
    end

    @testset "GenomicSequence" begin
        seq = dna"ACGTATGTNTNNCGAT"
        gseq = FMM.GenomicSequence(seq)
        @test length(gseq) == length(seq)
        for i in 1:length(seq)
            @test seq[i] === gseq[i]
        end
        # TODO: more tests for subseq
        subseq = FMM.subseq(gseq, 3, 0, false)
        @test isa(subseq, SIMDAlignment.seq_t)
    end

    @testset "Genome" begin
        mini = FMM.load_genome(Pkg.dir("FMM", "test", "mini.fasta"))
        @test length(mini) === 52
        @test collect(FMM.eachname(mini)) == ["chr1", "chr2", "chr3"]

        # chr1
        @test mini[1] === DNA_N
        @test mini[2] === DNA_N
        @test mini[5] === DNA_A
        @test mini[6] === DNA_C
        @test FMM.locus(mini, 1) == ("chr1", 1)
        @test FMM.locus(mini, 2) == ("chr1", 2)
        @test FMM.locus(mini, 26) == ("chr1", 26)
        @test FMM.locus(mini, 27) == ("chr1", 27)
        @test length(mini, "chr1") === 27

        # chr2
        @test mini[28] === DNA_N
        @test mini[31] === DNA_A
        @test mini[32] === DNA_T
        @test FMM.locus(mini, 28) == ("chr2", 1)
        @test FMM.locus(mini, 31) == ("chr2", 4)
        @test length(mini, "chr2") === 16

        # chr3
        @test mini[44] === DNA_A
        @test mini[52] === DNA_T
        @test FMM.locus(mini, 44) == ("chr3", 1)
        @test FMM.locus(mini, 52) == ("chr3", 9)
        @test length(mini, "chr3") === 9
    end

    @testset "GenomeIndex" begin
        k = 4
        index = FMM.build_index(Pkg.dir("FMM", "test", "mini.fasta"), k)

        sarange = index.sa_table[kmer(dna"GATT")]
        @test isa(sarange, UnitRange)
        @test length(sarange) === 1
        @test count(dna"GATT",  index.fmindex) === 1
        @test count(dna"GATTT", index.fmindex) === 1
        @test count(dna"GATTG", index.fmindex) === 0

        sarange = index.sa_table[kmer(dna"TTGA")]
        @test length(sarange) === 2
        @test count(dna"TTGA",  index.fmindex) === 2
        @test count(dna"TTGAC", index.fmindex) === 1
        @test count(dna"TTGAG", index.fmindex) === 1
    end

    @testset "AlignmentProfile" begin
        profile = FMM.AlignmentProfile(
            seed_interval=10,
            max_seed_hit=16,
            max_seed_try=4,
            matching_score=0,
            mismatching_score=-3,
            gap_open_penalty=5,
            gap_extend_penalty=3
        )
        @test profile.seed_interval === 10
        @test profile.max_seed_hit === 16
        @test profile.max_seed_try === 4
        @test profile.score_model.submat[DNA_A,DNA_A] ===  0
        @test profile.score_model.submat[DNA_A,DNA_C] === -3
        @test profile.score_model16.submat[DNA_A,DNA_A] === Int16(0)
        @test profile.score_model16.submat[DNA_A,DNA_C] === Int16(-3)
    end
end
