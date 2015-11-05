type AlignedRead
    record::Bio.Seq.SeqRecord
    seq::AlignedSequence
    forward::Bool
    mapquality::Int
end

name(r::AlignedRead) = r.record.name
firstpos(r::AlignedRead) = Bio.Align.first(r.seq)
sequence(r::AlignedRead) = r.seq.seq
Bio.Align.cigar(r::AlignedRead) = cigar(r.seq.aln)
isforward(r::AlignedRead) = r.forward
mapquality(r::AlignedRead) = r.mapquality

# type PairedAlignedRead
# end
