#
#   submit.jl
#
#   Created by Petter Taule on 16.08.2021
#   Copyright (c) 2021 Petter Taule. All rights reserved.
#

using ArgParse

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--k_idx_lower"
            help = "lower k index"
            arg_type = Int
            default = 0
        "--k_idx_upper"
            help = "upper k index"
            arg_type = Int
            default = 0
        "--z_a_idx"
            help = "First redshift index"
            arg_type = Int
            default = -1
        "--z_b_idx"
            help = "Last redshift index"
            arg_type = Int
            default = -1
        "--log_dir"
            help = "log directory"
            arg_type = String
            default = "/space/ge52sir/sge_output/"
        "m_nu"
            help = "neutrino mass"
            arg_type = Float64
            required = true
        "class_dir"
            help = "CLASS directory to read metric potential from"
            arg_type = String
            required = true
        "output_dir"
            help = "CLASS directory to read metric potential from"
            arg_type = String
            required = true
    end

    return parse_args(s)
end

function main()
    parse_args = parse_commandline()

    z_a_idx     = parse_args["z_a_idx"]
    z_b_idx     = parse_args["z_b_idx"]
    k_idx_lower = parse_args["k_idx_lower"]
    k_idx_upper = parse_args["k_idx_upper"]
    m_nu        = parse_args["m_nu"]
    class_dir   = parse_args["class_dir"]
    output_dir  = parse_args["output_dir"]
    log_dir     = parse_args["log_dir"]

    if z_a_idx == -1
        z_a_idx = 0
        z_b_idx = 49
    end
    if z_b_idx == -1
        z_b_idx = z_a_idx
    end

    @sync for k_idx = k_idx_lower : k_idx_upper
        job_name = "sk_$(m_nu)_$(k_idx)_$(z_a_idx)"
        log_file = log_dir * job_name * ".log"

        qsub_cmd = `qsub -N $job_name -pe smp 1 -cwd -e $log_dir/error/
        -o $log_dir/output/ run.sh --k_idx=$k_idx --z_a_idx=$z_a_idx
        --z_b_idx=$z_b_idx $m_nu $class_dir $output_dir $log_file`

        @async run(qsub_cmd)
    end
end

main()
