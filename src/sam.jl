type SAMWriter{IOtype<:IO}
    io::IOtype
    genome::Genome
end

# replace all consective spaces with a tab
macro tab_str(s)
    replace(s, r"\s+", '\t')
end

function writeheader(sw::SAMWriter)
    # @HD: the header line
    write(sw.io, "@HD", '\t', "VN:1.5", '\t', "SO:unknown", '\n')

    # @SQ: reference sequences
    for name in eachname(sw.genome)
        len = length(sw.genome, name)
        print(sw.io, "@SQ", '\t', "SN:", name, '\t', "LN:", len, '\n')
    end

    # @PG: program
    write(sw.io, "@PG", '\t', "ID:FMM.jl", '\n')
end

# FLAGs
const MultipleSeguemtns = 0x0001
const ProperlyAligned   = 0x0002
const SegmentUnmapped   = 0x0004
const NextSegment       = 0x0008
const ReverseComplement = 0x0010
const NextReverseComplement = 0x0020
const FirstSegment      = 0x0040
const SecondSegment     = 0x0080
const SecondaryAlignment = 0x0100
const NotPassed         = 0x0200
const Duplicate         = 0x0400
const Supplementary     = 0x0800

# mapped read
function Base.write(sw::SAMWriter, read::AlignedRead)
    # QNAME
    write(sw.io, name(read), '\t')
    # FLAG
    flag = ProperlyAligned
    if !isforward(read)
        flag |= ReverseComplement
    end
    print(sw.io, flag, '\t')  # TODO
    # RNAME and POS
    rname, pos = locus(sw.genome, firstpos(read))
    writetab(sw.io, rname)
    print(sw.io, pos, '\t')
    # MAPQ
    writetab(sw.io, '*')  # TODO
    # CIGAR
    writetab(sw.io, cigar(read))
    # RNEXT
    writetab(sw.io, '*')  # TODO
    # PNEXT
    writetab(sw.io, '0')  # TODO
    # TLEN
    writetab(sw.io, '0')  # TODO
    # SEQ
    writeseq(sw.io, sequence(read))
    # QUAL
    writetab(sw.io, '*')  # TODO
    write(sw.io, '\n')
end

# unmapped read
function Base.write(sw::SAMWriter, read::Bio.Seq.SeqRecord)
    # QNAME
    writetab(sw.io, read.name)
    # FLAG
    flag = SegmentUnmapped
    print(sw.io, flag, '\t')  # TODO
    # RNAME
    writetab(sw.io, '*')
    # POS
    print(sw.io, 0, '\t')
    # MAPQ
    writetab(sw.io, '*')  # TODO
    # CIGAR
    writetab(sw.io, '*')
    # RNEXT
    writetab(sw.io, '*')  # TODO
    # PNEXT
    writetab(sw.io, '0')  # TODO
    # TLEN
    writetab(sw.io, '0')  # TODO
    # SEQ
    writeseq(sw.io, read.seq)
    # QUAL
    writetab(sw.io, '*')  # TODO
    write(sw.io, '\n')
end

function writetab(io::IO, x)
    write(io, x, '\t')
end

function writeseq(io::IO, seq)
    for i in 1:length(seq)
        show(io, seq[i])
    end
    write(io, '\t')
end
