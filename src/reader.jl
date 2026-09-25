# BigBed Reader
# =============

import Automa.Stream: @mark, @markpos, @relpos, @abspos
import CodecZlib.TranscodingStreams: TranscodingStreams, TranscodingStream, NoopStream

function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    unsafe_copyto!(dst, dpos, src, spos, n)
    return dst
end

struct Reader <: BioGenerics.IO.AbstractReader
    stream::IO
    header::BBI.Header
    zooms::Vector{BBI.Zoom}
    summary::BBI.Summary
    # chrom name => (ID, length)
    chroms::Dict{String,Tuple{UInt32,Int}}
    # chrom ID => name
    chrom_names::Dict{UInt8,String}
    # data index
    index::BBI.RTree
end

function Base.eltype(::Type{Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.stream
end

"""
    BigBed.Reader(input::IO)

Create a reader for bigBed file format.

Note that `input` must be seekable.
"""
function Reader(input::IO)
    # read header
    header = read(input, BBI.Header)
    if header.magic != BBI.BED_MAGIC
        error("invalid BigBed magic")
    elseif header.version < 3
        error("not a supported version of BigBed")
    end
    # read zoom objects
    zoom_headers = Vector{BBI.ZoomHeader}(undef, header.zoom_levels)
    read!(input, zoom_headers)
    zooms = [BBI.Zoom(input, h, header.uncompress_buf_size) for h in zoom_headers]
    sort!(zooms, by=z->z.header.reduction_level)
    # read summary, B tree, and R tree
    seek(input, header.total_summary_offset)
    summary = read(input, BBI.Summary)
    chromindex = BBI.BTree(input, header.chromosome_tree_offset)
    chroms = Dict(name => (id, Int(len)) for (name, id, len) in BBI.chromlist(chromindex))
    chrom_names = Dict(id => name for (name, (id, len)) in chroms)
    index = BBI.RTree(input, header.full_index_offset)
    return Reader(input, header, zooms, summary, chroms, chrom_names, index)
end

"""
    chromlist(reader::BigBed.Reader)::Vector{Tuple{String,Int}}

Get the `(name, length)` pairs of chromosomes/contigs.
"""
function chromlist(reader::Reader)::Vector{Tuple{String,Int}}
    return sort!([(name, Int(len)) for (name, (id, len)) in reader.chroms], by=x->x[1])
end

const data_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    opt = Automa.RegExp.opt

    record = let
        chromid = re"...."
        chromid.actions[:exit] = [:record_chromid]

        chromstart = re"...."
        chromstart.actions[:exit] = [:record_chromstart]

        chromend = re"...."
        chromend.actions[:exit] = [:record_chromend]

        name = re"[ -~]*"
        name.actions[:enter] = [:mark]
        name.actions[:exit]  = [:record_name]

        score = re"[0-9]+"
        score.actions[:enter] = [:mark]
        score.actions[:exit]  = [:record_score]

        strand = re"[+\-.?]"
        strand.actions[:enter] = [:record_strand]

        thickstart = re"[0-9]+"
        thickstart.actions[:enter] = [:mark]
        thickstart.actions[:exit]  = [:record_thickstart]

        thickend = re"[0-9]+"
        thickend.actions[:enter] = [:mark]
        thickend.actions[:exit]  = [:record_thickend]

        itemrgb = cat(re"[0-9]+", opt(cat(',', re"[0-9]+", ',', re"[0-9]+")))
        itemrgb.actions[:enter] = [:mark]
        itemrgb.actions[:exit]  = [:record_itemrgb]

        blockcount = re"[0-9]+"
        blockcount.actions[:enter] = [:mark]
        blockcount.actions[:exit]  = [:record_blockcount]

        # comma-separated values
        csv(x) = cat(rep(cat(x, ',')), opt(x))

        blocksizes = let
            blocksize = re"[0-9]+"
            blocksize.actions[:enter] = [:mark]
            blocksize.actions[:exit]  = [:record_blocksizes_blocksize]

            csv(blocksize)
        end
        blocksizes.actions[:exit] = [:record_blocksizes]

        blockstarts = let
            blockstart = re"[0-9]+"
            blockstart.actions[:enter] = [:mark]
            blockstart.actions[:exit]  = [:record_blockstarts_blockstart]

            csv(blockstart)
        end
        blockstarts.actions[:exit] = [:record_blockstarts]

        cat(
            chromid, chromstart, chromend,
            opt(cat(      name,
            opt(cat('\t', score,
            opt(cat('\t', strand,
            opt(cat('\t', thickstart,
            opt(cat('\t', thickend,
            opt(cat('\t', itemrgb,
            opt(cat('\t', blockcount,
            opt(cat('\t', blocksizes,
            opt(cat('\t', blockstarts)))))))))))))))))))
    end
    record.actions[:enter] = [:anchor]
    record.actions[:exit]  = [:record]

    data = rep(cat(record, '\0'))

    return Automa.compile(data)
end)()

mutable struct Record
    chromid::UInt32
    chromstart::UInt32
    chromend::UInt32
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # number of columns
    ncols::Int
    # indexes
    name::UnitRange{Int}
    score::UnitRange{Int}
    strand::Int
    thickstart::UnitRange{Int}
    thickend::UnitRange{Int}
    itemrgb::UnitRange{Int}
    blockcount::UnitRange{Int}
    blocksizes::Vector{UnitRange{Int}}
    blockstarts::Vector{UnitRange{Int}}
    # reader
    reader::Reader

    function Record(chromid, chromstart, chromend,
                    data, filled, ncols,
                    name, score, strand, thickstart, thickend,
                    itemrgb, blockcount, blocksizes, blockstarts)
        return new(chromid, chromstart, chromend,
                   data, filled, ncols,
                   name, score, strand, thickstart, thickend,
                   itemrgb, blockcount, blocksizes, blockstarts)
    end
end

const actions = Dict(
    :anchor                        => :(@mark),
    :record_chromid                => :(record.chromid    = unsafe_load(convert(Ptr{UInt32}, pointer(data, p-4))); record.ncols += 1),
    :record_chromstart             => :(record.chromstart = unsafe_load(convert(Ptr{UInt32}, pointer(data, p-4))); record.ncols += 1),
    :record_chromend               => :(record.chromend   = unsafe_load(convert(Ptr{UInt32}, pointer(data, p-4))); record.ncols += 1),
    :mark                          => :(pos = @relpos(p)),
    :record_name                   => :(record.name       = (pos:@relpos(p-1)); record.ncols += 1),
    :record_score                  => :(record.score      = (pos:@relpos(p-1)); record.ncols += 1),
    :record_strand                 => :(record.strand     = @relpos(p); record.ncols += 1),
    :record_thickstart             => :(record.thickstart = (pos:@relpos(p-1)); record.ncols += 1),
    :record_thickend               => :(record.thickend   = (pos:@relpos(p-1)); record.ncols += 1),
    :record_itemrgb                => :(record.itemrgb    = (pos:@relpos(p-1)); record.ncols += 1),
    :record_blockcount             => :(record.blockcount = (pos:@relpos(p-1)); record.ncols += 1),
    :record_blocksizes_blocksize   => :(push!(record.blocksizes, (pos:@relpos(p-1)))),
    :record_blocksizes             => :(record.ncols += 1),
    :record_blockstarts_blockstart => :(push!(record.blockstarts, (pos:@relpos(p-1)))),
    :record_blockstarts            => :(record.ncols += 1),
    :record                        => quote
        appendfrom!(record.data, 1, data, @markpos, p - @markpos)
        record.filled = 1:(p - @markpos)
        found_record = true
        @escape
    end,
    :countrecord                   => :()
)

const initcode = quote
    pos = 0
    linenum = 0
    found_record = false
    cs, linenum = state
end

const loopcode = quote
    if found_record
        @goto __return__
    end
end

Automa.Stream.generate_reader(
    :_readrecord!,
    data_machine,
    arguments = (:(record::Record), :(state::Tuple{Int,Int})),
    actions = actions,
    initcode = initcode,
    loopcode = loopcode,
    returncode = :(return cs, linenum, found_record)
) |> eval

function _read!(stream::TranscodingStream, state::BioGenerics.Automa.State, record::Record)
    initialize!(record)
    cs, ln, found = _readrecord!(stream, record, (state.state, state.linenum))
    state.state = cs
    state.linenum = ln
    state.filled = found
    if found
        return record
    end
    if cs == 0 || eof(stream)
        throw(EOFError())
    end
    throw(ArgumentError("malformed BigBed file"))
end

# Iterator
# --------

mutable struct IteratorState
    state::BioGenerics.Automa.State
    record::Record
    n_records::UInt64
    current_record::UInt64
end

function Base.iterate(reader::Reader)
    seek(reader.stream, reader.header.full_data_offset)
    # this is defined as UInt32 in the spces but actually UInt64
    record_count = read(reader.stream, UInt64)
    datastream = CodecZlib.ZlibDecompressorStream(reader.stream)
    parser_state = BioGenerics.Automa.State(datastream, 1, 1, false)
    state = IteratorState(parser_state, Record(), record_count, 0)
    return iterate(reader, state)
end

function Base.iterate(reader::Reader, state::IteratorState)
    if state.current_record < state.n_records
        _read!(state.state.stream, state.state, state.record)
        state.record.reader = reader
        state.current_record += 1
        return copy(state.record), state
    end
    return nothing
end
