using StatsBase: WeightVec, sample, sample!

function run_alignment(profile::AlignmentProfile, index, read_file)
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
        println(rec.name)
        if isaligned(readstate)
            aln = alignment(readstate)
            chr, loc = locus(index.genome, Bio.Align.first(aln))
            println(chr, ':', loc)
            println(aln)
            println(cigar(aln.aln))
        else
            println("not aligned")
        end
        println()
    end
    info("finished: ", t, " s")
    info(@sprintf("%.1f", n_reads / t), " reads/s")
    # profile
    #reads = open(read_file, format)
    #readstate = ReadState()
    #@profile for rec in reads
    #    setread!(readstate, rec.seq)
    #    align_read!(readstate, index, profile)
    #end
    #Profile.print(STDERR, format=:flat, cols=1000)
end

function align_read!{T,k}(rs::ReadState, index::GenomeIndex{T,k}, profile)
    seedlen = profile.seed_length
    interval = profile.seed_interval

    # sarch seed hits
    search_seed!(rs, true,  index, seedlen, interval)
    search_seed!(rs, false, index, seedlen, interval)

    if !hashit(rs)
        # no clue
        return rs
    end

    # find best alignment from matching seeds
    for seedhit in each_forward_seedhit(rs)
        score_seed!(rs, seedhit, index, profile.score_model)
    end
    for seedhit in each_reverse_seedhit(rs)
        score_seed!(rs, seedhit, index, profile.score_model)
    end

    if !has_aligned_seed(rs)
        return rs
    end

    aln = align_hit(rs, index.genome, profile.score_model)
    set_alignment!(rs, aln)
    return rs
end

function search_seed!{T,k}(rs, forward, index::GenomeIndex{T,k}, seedlen, interval)
    read = forward ? forward_read(rs) : reverse_read(rs)
    println(forward ? "forward" : "backward")
    for s in endof(read):-interval:k
        if hasn(read[s-k+1:s])
            continue
        end
        kmer::DNAKmer{k} = read[s-k+1:s]
        sa_range::UnitRange{Int} = index.sa_table[kmer]
        @show length(sa_range)
        i = s - k
        while i ≥ 1 && read[i] != DNA_N && length(sa_range) > 16
            sa_range = FMIndexes.sa_range(read[i:i], index.fmindex, sa_range)
            i -= 1
        end
        if 1 ≤ length(sa_range) ≤ 16
            seed_range = i+1:s
            push!(rs, SeedHit(seed_range, sa_range, forward))
        end
    end
end

# calculate alignment scores of all hit locations
function score_seed!(rs::ReadState, seedhit::SeedHit, index, model)
    n_hits = count(seedhit)
    indices = 1:n_hits
    offset = 0
    read = isforward(seedhit) ? forward_read(rs) : reverse_read(rs)

    locs = Vector{Int}(n_hits)
    for i in 1:n_hits
        locs[i] = FMIndexes.sa_value(seedhit[indices[i]], index.fmindex) + 1
    end

    refseqs = Vector{seq_t}(n_hits)

    # align left sequences
    for i in 1:n_hits
        len = seed_start(seedhit) - 1
        offset = locs[i] - 2
        refseqs[i] = subseq(index.genome, len, offset, true)
    end

    left_scores = Vector{Score}(n_hits)
    right_scores = Vector{Score}(n_hits)

    alns = paralign_score(
        model.submat,
        model.gap_open_penalty,
        model.gap_extend_penalty,
        read[1:seed_start(seedhit)-1],
        refseqs
    )
    for i in 1:n_hits
        left_scores[i] = alns[i].score
    end

    # align right sequences
    for i in 1:n_hits
        len = length(read) - seed_stop(seedhit)
        offset = locs[i] + seed_length(seedhit) - 1
        refseqs[i] = subseq(index.genome, len, offset, false)
    end

    alns = paralign_score(
        model.submat,
        model.gap_open_penalty,
        model.gap_extend_penalty,
        read[seed_stop(seedhit)+1:end],
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
