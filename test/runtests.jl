using Test
using BigBed
using Distributions
using GenomicFeatures

import Random
import ColorTypes: RGB
import YAML
import FixedPointNumbers: N0f8

import BioCore:
    leftposition,
    hasleftposition,
    rightposition,
    hasrightposition,
    seqname,
    hasseqname


function get_bio_fmt_specimens(commit="222f58c8ef3e3480f26515d99d3784b8cfcca046")
    path = joinpath(dirname(@__FILE__), "BioFmtSpecimens")
    if !isdir(path)
        run(`git clone https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
    end
    cd(path) do
        #run(`git checkout $(commit)`)
    end
    return path
end

# Generate an array of n random Interval{Int} object.
# With sequence names samples from seqnames, and intervals drawn to lie in [1, maxpos].
function random_intervals(seqnames, maxpos::Int, n::Int)
    seq_dist = Categorical(length(seqnames))
    strand_dist = Categorical(2)
    length_dist = Normal(1000, 1000)
    intervals = Vector{Interval{Int}}(undef, n)
    for i in 1:n
        intlen = maxpos
        while intlen >= maxpos || intlen <= 0
            intlen = ceil(Int, rand(length_dist))
        end
        first = rand(1:maxpos-intlen)
        last = first + intlen - 1
        strand = rand(strand_dist) == 1 ? STRAND_POS : STRAND_NEG
        intervals[i] = Interval{Int}(seqnames[rand(seq_dist)], first, last, strand, i)
    end
    return intervals
end

@testset "BigBed" begin
    @testset "empty" begin
        buffer = IOBuffer()
        data = buffer.data
        writer = BigBed.Writer(buffer, [("chr1", 1000)])
        close(writer)
        reader = BigBed.Reader(IOBuffer(data))
        @test length(collect(reader)) == 0
    end

    @testset "small" begin
        buffer = IOBuffer()
        data = buffer.data
        writer = BigBed.Writer(buffer, [("chr1", 1000)])
        write(writer, ("chr1", 50, 100, "name1"))
        close(writer)
        reader = BigBed.Reader(IOBuffer(data))
        records = collect(reader)
        @test length(records) == 1
        @test BigBed.chrom(records[1]) == "chr1"
        @test BigBed.chromstart(records[1]) === 50
        @test BigBed.chromend(records[1]) === 100
        @test BigBed.name(records[1]) == "name1"
        @test !BigBed.hasscore(records[1])
        @test BigBed.optionals(records[1]) == ["name1"]

        buffer = IOBuffer()
        data = buffer.data
        writer = BigBed.Writer(buffer, [("chr1", 1000)])
        write(writer, ("chr1", 1, 100, "some name", 100, '+', 10, 90, RGB(0.5, 0.1, 0.2), 2, [4, 10], [10, 20]))
        close(writer)
        reader = BigBed.Reader(IOBuffer(data))
        records = collect(reader)
        @test length(records) == 1
        @test BigBed.haschrom(records[1]) === hasseqname(records[1]) === true
        @test BigBed.chrom(records[1]) == seqname(records[1]) == "chr1"
        @test BigBed.haschromstart(records[1]) === hasleftposition(records[1]) === true
        @test BigBed.chromstart(records[1]) === leftposition(records[1]) === 1
        @test BigBed.haschromend(records[1]) === hasrightposition(records[1]) === true
        @test BigBed.chromend(records[1]) === rightposition(records[1]) === 100
        @test BigBed.hasname(records[1])
        @test BigBed.name(records[1]) == "some name"
        @test BigBed.hasscore(records[1])
        @test BigBed.score(records[1]) === 100
        @test BigBed.hasstrand(records[1])
        @test BigBed.strand(records[1]) === STRAND_POS
        @test BigBed.hasthickstart(records[1])
        @test BigBed.thickstart(records[1]) === 10
        @test BigBed.hasthickend(records[1])
        @test BigBed.thickend(records[1]) === 90
        @test BigBed.hasitemrgb(records[1])
        @test BigBed.itemrgb(records[1]) === convert(RGB{N0f8}, RGB(0.5, 0.1, 0.2))
        @test BigBed.hasblockcount(records[1])
        @test BigBed.blockcount(records[1]) === 2
        @test BigBed.hasblocksizes(records[1])
        @test BigBed.blocksizes(records[1]) == [4, 10]
        @test BigBed.hasblockstarts(records[1])
        @test BigBed.blockstarts(records[1]) == [10, 20]
        @test BigBed.optionals(records[1]) == ["some name", "100", "+", "9", "90", "128,26,51", "2", "4,10,", "9,19,"]
    end

    @testset "large" begin
        buffer = IOBuffer()
        data = buffer.data
        binsize = 32
        writer = BigBed.Writer(buffer, [("chr1", 100_000), ("chr2", 90_000)], binsize=binsize)
        for i in 1:10_000
            write(writer, ("chr1", (i-1)*10+1, i*10, string("name", i)))
        end
        n = 0
        p = 1
        while p ≤ 90_000
            sz = min(rand(1:100), 90_000 - p)
            write(writer, ("chr2", p, p+sz, string("name", n + 1)))
            n += 1
            p += sz + 1
        end
        close(writer)
        reader = BigBed.Reader(IOBuffer(data))
        records = collect(reader)
        @test length(records) == 10_000 + n
        records = collect(eachoverlap(reader, Interval("chr1", 50_001, 50_165)))
        @test length(records) == 17
    end

    @testset "round trip" begin
        function test_round_trip(filepath)
            reader = open(BigBed.Reader, filepath)
            buffer = IOBuffer()
            data = buffer.data
            writer = BigBed.Writer(buffer, BigBed.chromlist(reader))
            original = []
            for record in reader
                t = (BigBed.chrom(record), BigBed.chromstart(record), BigBed.chromend(record), BigBed.optionals(record)...)
                write(writer, t)
                push!(original, t)
            end
            close(writer)
            close(reader)

            reader = BigBed.Reader(IOBuffer(data))
            copy = []
            for record in reader
                t = (BigBed.chrom(record), BigBed.chromstart(record), BigBed.chromend(record), BigBed.optionals(record)...)
                push!(copy, t)
            end
            close(reader)

            @test original == copy
        end

        dir = joinpath(get_bio_fmt_specimens(), "BBI")
        for specimen in YAML.load_file(joinpath(dir, "index.yml"))
            valid = get(specimen, "valid", true)
            bigbed = "bigbed" ∈ split(specimen["tags"])
            if valid && bigbed
                test_round_trip(joinpath(dir, specimen["filename"]))
            end
        end
    end

    @testset "overlap" begin
        chromlen = 1_000_000
        Random.seed!(1234)
        chroms = ["one", "two", "three", "four", "five"]
        intervals = IntervalCollection([Interval(i.seqname, i.first, i.last) for i in random_intervals(chroms, chromlen, 10_000)], true)

        buffer = IOBuffer()
        data = buffer.data
        writer = BigBed.Writer(buffer, [(chrom, chromlen) for chrom in chroms])
        for i in intervals
            write(writer, i)
        end
        close(writer)

        reader = BigBed.Reader(IOBuffer(data))
        queries = random_intervals(chroms, chromlen, 1000)
        triplet(x::Interval) = String(x.seqname), x.first, x.last
        triplet(x::BigBed.Record) = BigBed.chrom(x), BigBed.chromstart(x), BigBed.chromend(x)
        @test all(triplet.(collect(eachoverlap(intervals, q))) == triplet.(collect(eachoverlap(reader, q))) for q in queries)
        close(reader)
    end
end
