#
#   check_error.jl
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

    z_vals = readdlm(file_prefix * "_k_0.dat", comments=true)[:,1]

    for k=0:99
        data = readdlm(file_prefix * "_k_$k.dat", comments=true)

        if data[:,1] != z_vals
            println("z_column for k_index = $k does not equal to that of k_index = 0")
        end

        for z=1:size(data,1)
            rel_err = data[z,4] / data[z,2]
            if rel_err > 1e-2
                println("Relative error $rel_err for k_index = $k and z_index = $z")
            end
        end
    end
end

main()
