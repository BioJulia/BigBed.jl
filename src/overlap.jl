# BigBed Overlap
# ==============

struct OverlapIterator
    reader::Reader
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32
end

function Base.eltype(::Type{OverlapIterator})
    return Record
end

function Base.IteratorSize(::Type{OverlapIterator})
    return Base.SizeUnknown()
end

function GenomicFeatures.eachoverlap(reader::Reader, interval::GenomicFeatures.Interval)
    if haskey(reader.chroms, seqname(interval))
        id, _ = reader.chroms[seqname(interval)]
    else
        id = typemax(UInt32)
    end
    return OverlapIterator(reader, id, leftposition(interval) - 1, rightposition(interval))
end

mutable struct OverlapIteratorState
    state::BioGenerics.Automa.State
    data::Vector{UInt8}
    done::Bool
    record::Record
    blocks::Vector{BBI.Block}
    current_block::Int
end

function Base.iterate(iter::OverlapIterator)
    data = Vector{UInt8}(undef, iter.reader.header.uncompress_buf_size)
    blocks = BBI.find_overlapping_blocks(iter.reader.index, iter.chromid, iter.chromstart, iter.chromend)
    dummy_stream = NoopStream(IOBuffer())
    state = OverlapIteratorState(
        BioGenerics.Automa.State(dummy_stream, 1, 1, false),
        data,
        isempty(blocks), Record(), blocks, 1)
    return iterate(iter, state)
end

function Base.iterate(iter::OverlapIterator, state::OverlapIteratorState)
    advance!(iter, state)
    if state.done
        return nothing
    end
    return copy(state.record), state
end

function advance!(iter::OverlapIterator, state::OverlapIteratorState)
    while true
        while state.current_block ≤ lastindex(state.blocks) && eof(state.state.stream)
            block = state.blocks[state.current_block]
            seek(iter.reader.stream, block.offset)
            size = BBI.uncompress!(state.data, read(iter.reader.stream, block.size))
            state.state = BioGenerics.Automa.State(
                NoopStream(IOBuffer(view(state.data, 1:size))),
                1, 1, false)
            state.current_block += 1
        end
        if state.done || (state.current_block > lastindex(state.blocks) && eof(state.state.stream))
            state.done = true
            return state
        end

        _read!(state.state.stream, state.state, state.record)
        state.record.reader = iter.reader
        if overlaps(state.record, iter.chromid, iter.chromstart, iter.chromend)
            return state
        end
    end
end

function overlaps(record::Record, chromid::UInt32, chromstart::UInt32, chromend::UInt32)
    return record.chromid == chromid && !(record.chromend ≤ chromstart || record.chromstart ≥ chromend)
end
