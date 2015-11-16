function run_alignment(profile::AlignmentProfile, index, read_file, output)
    format = endswith(read_file, ".fa") ? FASTA :
             endswith(read_file, ".fq") ? FASTQ :
             error("unknown format")
    reads = open(read_file, format)
    readstate = ReadState()
    out = SAMWriter(output, index.genome)
    n_reads = 0
    info("aligning reads")
    writeheader(out)
    #Profile.init(delay=0.01)
    t = @elapsed for rec in reads
        setrecord!(readstate, rec)
        aln = align_read!(readstate, index, profile)
        n_reads += 1
        if isnull(aln)
            write(out, rec)
        else
            write(out, get(aln))
        end
    end
    #Profile.print(STDERR, format=:tree, cols=100000)
    info("finished: ", t, " s")
    info(@sprintf("%.1f", n_reads / t), " reads/s")
end

function align_read!{T,k}(rs::ReadState, index::GenomeIndex{T,k}, profile)
    # sarch seed hits
    seed_interval = profile.seed_interval
    max_seed_hit = profile.max_seed_hit
    search_seed!(rs,  true, index, seed_interval, max_seed_hit)
    search_seed!(rs, false, index, seed_interval, max_seed_hit)

    if !hashit(rs)
        # shift the seed window
        shift = div(seed_interval, 2)
        search_seed!(rs,  true, index, seed_interval, 4max_seed_hit, shift)
        search_seed!(rs, false, index, seed_interval, 4max_seed_hit, shift)
    end

    if !hashit(rs)
        search_seed!(rs,  true, index, 2seed_interval, 8max_seed_hit)
        search_seed!(rs, false, index, 2seed_interval, 8max_seed_hit)
    end

    if !hashit(rs)
        # no clue
        return Nullable()
    end

    # find best alignment from matching seeds
    best = typemin(Int)
    mis = mismatching_score(profile)
    achievable = readlen(rs) * maximum(profile.score_model.submat)
    ntry = 0
    for seedhit in each_prioritized_seedhit(rs, index)
        if best < 3mis
            ntry ≥ 2^3 * profile.max_seed_try && break
        elseif best < 2mis
            ntry ≥ 2^2 * profile.max_seed_try && break
        elseif best < 1mis
            ntry ≥ 2^1 * profile.max_seed_try && break
        else
            ntry ≥ profile.max_seed_try && break
        end
        score = score_seed!(rs, seedhit, index, profile.score_model16)
        if score ≥ achievable
            best = score
            break
        elseif score > best
            best = score
            ntry = 0
        else
            ntry += 1
        end
    end

    if !has_aligned_seed(rs)
        return Nullable()
    end

    bestseed = best_aligned_seed(rs)
    alnseq = align_hit(rs, bestseed, index.genome, profile.score_model)
    return Nullable(AlignedRead(record(rs), alnseq, isforward(bestseed), map_quality(rs, profile)))
end

function search_seed!{T,k}(rs, forward, index::GenomeIndex{T,k}, interval, maxseed, shift=0)
    read = forward ? forward_read(rs) : reverse_read(rs)
    for s in endof(read)-shift:-interval:k
        kmer = read[s-k+1:s]
        if hasn(kmer)
            continue
        end
        sa_range::UnitRange{Int} = index.sa_table[convert(DNAKmer{k}, kmer)]
        i = s - k
        while i > 1 && read[i] != DNA_N && length(sa_range) > maxseed
            byte::UInt8 = read[i]
            sa_range = FMIndexes.sa_range(byte, index.fmindex, sa_range)
            i -= 1
        end
        if 1 ≤ length(sa_range) ≤ maxseed
            seed_range = i+1:s
            push!(rs, SeedHit(seed_range, sa_range, forward))
        end
    end
end

# calculate alignment scores of all hit locations
function score_seed!(rs::ReadState, seedhit::SeedHit, index, model)
    n_hits = count(seedhit)
    read = isforward(seedhit) ? forward_read(rs) : reverse_read(rs)
    refseqs = Vector{seq_t}(n_hits)

    # align left sequences
    for i in 1:n_hits
        len = seed_start(seedhit) - 1
        offset = metadata(seedhit)[i] - 2
        refseqs[i] = subseq(index.genome, len, offset, true)
    end

    left_scores = Vector{Score}(n_hits)
    right_scores = Vector{Score}(n_hits)

    alns = paralign_score(
        model.submat,
        model.gap_open,
        model.gap_extend,
        seq_t(read[1:seed_start(seedhit)-1], true),
        refseqs
    )
    for i in 1:n_hits
        left_scores[i] = alns[i].score
    end

    # align right sequences
    for i in 1:n_hits
        len = length(read) - seed_stop(seedhit)
        offset = metadata(seedhit)[i] + seed_length(seedhit) - 1
        refseqs[i] = subseq(index.genome, len, offset, false)
    end

    alns = paralign_score(
        model.submat,
        model.gap_open,
        model.gap_extend,
        seq_t(read[seed_stop(seedhit)+1:end], false),
        refseqs
    )
    for i in 1:n_hits
        right_scores[i] = alns[i].score
    end

    best::Int = typemin(Int)
    for i in 1:n_hits
        loc = metadata(seedhit)[i]
        seedhitext = SeedHitExt(
            seedhit,
            loc,
            left_scores[i],
            0,  # TODO: fix
            right_scores[i]
        )
        best = max(total_score(seedhitext), best)
        push!(rs, seedhitext)
    end
    return best
end

function align_hit(rs::ReadState, seedhit, genome::Genome, affinegap)
    rfrag = Vector{DNANucleotide}()
    gfrag = Vector{DNANucleotide}()
    read_anchors = Vector{AlignmentAnchor}()

    # left
    loc = seedhit.location
    read = isforward(seedhit) ? forward_read(rs) : reverse_read(rs)
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

# FASTA
function map_quality(rs, profile)
    # base quality
    qual = 3
    L = readlen(rs)
    # infer the number of mismatches with the alignment score
    mis = mismatching_score(profile)
    bestseed = popmax!(rs.seedhit_queue)
    best = 10^(-div(score(bestseed), mis) * qual / 10)
    ∑mis = best + eps()
    while !isempty(rs.seedhit_queue)
        seedhit = popmax!(rs.seedhit_queue)
        if abs(seedhit.location - bestseed.location) < L
            continue
        end
        s = score(seedhit)
        ∑mis += 10^(-div(s, mis) * qual / 10)
        if s < 10mis
            break
        end
    end
    return ceil(Int, -10log10(1 - best/∑mis))
end
