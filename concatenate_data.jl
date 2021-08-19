#
#   concatenate_data.jl
#
#   Created by Petter Taule on 13.06.2020
#   Copyright (c) 2020 Petter Taule. All rights reserved.
#

using DelimitedFiles, ArgParse


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--n_files"
            help = "Number of files"
            arg_type = Int
            default = 100
        "file_prefix"
            help = "File prefix"
            required = true
    end

    return parse_args(s)
end

function main()
    parse_args = parse_commandline()

    file_prefix = parse_args["file_prefix"]
    n_files = parse_args["n_files"]

    data = readdlm(file_prefix * "_k_0.dat", comments=true)

    table = zeros(size(data,1), n_files)

    # Index naming of files is zero-based
    for k=0:n_files-1
        data = readdlm(file_prefix * "_k_$k.dat", comments=true)

        table[:, k+1] = data[:,2]
    end

    # Need increasing etaD grid (i.e. decreasing redshift), therefore reverse
    # first dimension of table
    table = reverse(table, dims=1)

    # Write results to file
    writedlm(file_prefix * ".dat", table, "\t\t")
end

main()
