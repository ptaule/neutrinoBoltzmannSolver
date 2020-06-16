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
        "file_prefix"
            help = "File prefix"
            required = true
    end

    return parse_args(s)
end

function main()
    parse_args = parse_commandline()

    file_prefix = parse_args["file_prefix"]

    data = readdlm(file_prefix * "_k_0.dat", comments=true)

    table = data[:,2]

    for k=1:99
        data = readdlm(file_prefix * "_k_$k.dat", comments=true)

        table = [table data[:,2]]
    end
    writedlm(file_prefix * ".dat", table, "\t\t")
end

main()
