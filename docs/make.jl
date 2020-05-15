using Pkg
using Documenter, BigBed

makedocs(
    format = Documenter.HTML(
        edit_link = "develop"
    ),
    modules = [BigBed],
    sitename = "BigBed.jl",
    pages = [
        "Home" => "index.md",
        "BigBed" => "man/bigbed.md",
        "API Reference" => "man/api.md"
    ],
    authors = replace(join(Pkg.TOML.parsefile("Project.toml")["authors"], ", "), r" <.*?>" => "" ) * ", The BioJulia Organisation, and other contributors."

)

deploydocs(
    repo = "github.com/BioJulia/BigBed.jl.git",
    devbranch = "develop",
    push_preview = true
)
