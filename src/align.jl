using StatsBase: WeightVec, sample, sample!

function run_alignment(profile::AlignmentProfile, index, read_file, output)
    format = endswith(read_file, ".fa") ? FASTA :
             endswith(read_file, ".fq") ? FASTQ :
             error("unknown format")
    reads = open(read_file, format)
    readstate = ReadState()
    n_reads = 0
    info("aligning reads")
    t = @elapsed for rec in reads
        setread!(readstate, rec.seq)
        align_read!(readstate, index, profile)
        n_reads += 1
        println(output, rec.name, " ", rec.metadata)
        if isaligned(readstate)
            aln = alignment(readstate)
            chr, loc = locus(index.genome, Bio.Align.first(aln))
            println(output, chr, ':', loc)
            println(output, aln)
            println(output, cigar(aln.aln))
        else
            println(output, "not aligned")
        end
        println(output)
    end
    info("finished: ", t, " s")
    info(@sprintf("%.1f", n_reads / t), " reads/s")
end

function align_read!{T,k}(rs::ReadState, index::GenomeIndex{T,k}, profile)
    # sarch seed hits
    search_seed!(rs,  true, index, profile.seed_interval, profile.seed_extend_length, profile.max_seed_extension)
    search_seed!(rs, false, index, profile.seed_interval, profile.seed_extend_length, profile.max_seed_extension)

    if !hashit(rs)
        # no clue
        return rs
    end

    # find best alignment from matching seeds
    for seedhit in each_forward_seedhit(rs)
        score_seed!(rs, seedhit, index, profile.seed_extend_length, profile.score_model8)
    end
    for seedhit in each_reverse_seedhit(rs)
        score_seed!(rs, seedhit, index, profile.seed_extend_length, profile.score_model8)
    end

    if !has_aligned_seed(rs)
        return rs
    end

    aln = align_hit(rs, index.genome, profile.score_model)
    set_alignment!(rs, aln)
    return rs
end

function search_seed!{T,k}(rs, forward, index::GenomeIndex{T,k}, interval, extlen, maxseed)
    read = forward ? forward_read(rs) : reverse_read(rs)
    L = extlen
    for s in endof(read)-L:-interval:k+L
        kmer = read[s-k+1:s]
        if hasn(kmer)
            continue
        end
        sa_range::UnitRange{Int} = index.sa_table[convert(DNAKmer{k}, kmer)]
        i = s - k
        nt = read[i]
        while i > L && nt != DNA_N && length(sa_range) > maxseed
            byte::UInt8 = nt
            sa_range = FMIndexes.sa_range(byte, index.fmindex, sa_range)
            i -= 1
            nt = read[i]
        end
        if 1 ≤ length(sa_range) ≤ maxseed
            seed_range = i+1:s
            push!(rs, SeedHit(seed_range, sa_range, forward))
        end
    end
end

# calculate alignment scores of all hit locations
function score_seed!(rs::ReadState, seedhit::SeedHit, index, extlen, model)
    n_hits = count(seedhit)
    read = isforward(seedhit) ? forward_read(rs) : reverse_read(rs)

    locs = Vector{Int}(n_hits)
    for i in 1:n_hits
        locs[i] = FMIndexes.sa_value(seedhit[i], index.fmindex) + 1
    end

    refseqs = Vector{seq_t}(n_hits)

    # align left sequences
    for i in 1:n_hits
        offset = locs[i] - 2
        refseqs[i] = subseq(index.genome, extlen, offset, true)
    end

    left_scores = Vector{Score}(n_hits)
    right_scores = Vector{Score}(n_hits)

    alns = paralign_score(
        model.submat,
        model.gap_open_penalty,
        model.gap_extend_penalty,
        seq_t(read[seed_start(seedhit)-extlen:seed_start(seedhit)-1], true),
        refseqs
    )
    for i in 1:n_hits
        left_scores[i] = alns[i].score
    end

    # align right sequences
    for i in 1:n_hits
        offset = locs[i] + seed_length(seedhit) - 1
        refseqs[i] = subseq(index.genome, extlen, offset, false)
    end

    alns = paralign_score(
        model.submat,
        model.gap_open_penalty,
        model.gap_extend_penalty,
        seq_t(read[seed_stop(seedhit)+1:seed_stop(seedhit)+extlen], false),
        refseqs
    )
    for i in 1:n_hits
        right_scores[i] = alns[i].score
    end

    for i in 1:n_hits
        seedhitext = SeedHitExt(
            seedhit,
            locs[i],
            left_scores[i],
            0,  # TODO: fix
            right_scores[i]
        )
        push!(rs, seedhitext)
    end
end

function align_hit(rs::ReadState, genome::Genome, affinegap)
    seedhit = best_aligned_seed(rs)
    loc = seedhit.location
    read = isforward(seedhit) ? forward_read(rs) : reverse_read(rs)
    rfrag = Vector{DNANucleotide}()
    gfrag = Vector{DNANucleotide}()
    read_anchors = Vector{AlignmentAnchor}()

    # left
    unpack_seq!(rfrag, read, seed_start(seedhit) - 1, seed_start(seedhit) - 1, true)
    unpack_seq!(gfrag, genome, loc - 1, length(rfrag), true)
    aln = pairalign(GlobalAlignment(), rfrag, gfrag, affinegap)
    seqpos = 0
    refpos = loc - length(gfrag) - 1
    push!(read_anchors, AlignmentAnchor(seqpos, refpos, OP_START))
    anchors = Bio.Align.alignment(aln).anchors
    for i in endof(anchors):-1:2
        Δseqpos, Δrefpos = anchorwidth(anchors, i)
        seqpos += Δseqpos
        refpos += Δrefpos
        push!(read_anchors, AlignmentAnchor(seqpos, refpos, anchors[i].op))
    end

    # match
    seqpos += seed_length(seedhit)
    refpos += seed_length(seedhit)
    push!(read_anchors, AlignmentAnchor(seqpos, refpos, OP_SEQ_MATCH))

    # right
    unpack_seq!(rfrag, read, seed_stop(seedhit) + 1, length(read) - seed_stop(seedhit), false)
    unpack_seq!(gfrag, genome, loc + seed_length(seedhit), length(rfrag), false)
    aln = pairalign(GlobalAlignment(), rfrag, gfrag, affinegap)
    anchors = Bio.Align.alignment(aln).anchors
    for i in 2:endof(anchors)
        Δseqpos, Δrefpos = anchorwidth(anchors, i)
        seqpos += Δseqpos
        refpos += Δrefpos
        push!(read_anchors, AlignmentAnchor(seqpos, refpos, anchors[i].op))
    end

    return AlignedSequence(read, compress(read_anchors))
end

function anchorwidth(anchors, i)
    op = anchors[i].op
    if ismatchop(op)
        Δseqpos = anchors[i].seqpos - anchors[i-1].seqpos
        Δrefpos = anchors[i].refpos - anchors[i-1].refpos
        @assert Δseqpos == Δrefpos
    elseif isinsertop(op)
        Δseqpos = anchors[i].seqpos - anchors[i-1].seqpos
        Δrefpos = 0
    elseif isdeleteop(op)
        Δseqpos = 0
        Δrefpos = anchors[i].refpos - anchors[i-1].refpos
    else
        error("invalid op: $op")
    end
    return Δseqpos, Δrefpos
end

function compress(anchors)
    compressed = Vector{AlignmentAnchor}()
    op = anchors[1].op
    for i in 2:endof(anchors)
        if anchors[i].op != op
            push!(compressed, anchors[i-1])
            op = anchors[i].op
        end
    end
    push!(compressed, anchors[end])
    return compressed
end

function unpack_seq!(dst, src, start, len, reversed)
    resize!(dst, len)
    j = start
    for i in 1:len
        dst[i] = src[j]
        j += reversed ? -1 : +1
    end
    return dst
end
