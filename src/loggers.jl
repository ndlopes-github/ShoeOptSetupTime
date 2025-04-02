#= Copyright (C) 2024
Nuno David Lopes.
Created:  2025/03/13
Last changed - N. Lopes:2025/03/13 19:55:07
=#

using DrWatson
@quickactivate "ShoeOptSetupTime"
using LoggingExtras, Dates

function format_logger(prefix, level="info")
    FormatLogger(prefix * "_" * level * ".log.txt"; append=true) do io, args
        println(io, "[", args.level, "] ", args._module, ":[", rpad(args.line, 3, " "), "] | ", args.message)
    end
end

# Define a custom logger with level filtering
# Define a basic formatter to log only time, level, and message

function demux_logger(prefix)
    dt=now()
    prefixdt = prefix*"_"*Dates.format(dt, "yyyy-mm-dd_HH-MM-SS")
    TeeLogger(
        EarlyFilteredLogger(format_logger(prefixdt, "info")) do log_args
            log_args.level == Logging.Info
        end,
        EarlyFilteredLogger(format_logger(prefixdt, "debug")) do log_args
            log_args.level == Logging.Debug
        end,
        EarlyFilteredLogger(format_logger(prefixdt, "warn")) do log_args
            log_args.level == Logging.Warn
        end,
    )
end

# New ConsoleLogger that prints to stderr and accept messages 
#log_level = Logging.Info #Define logging level: Logging.Debug Logging.Info 
#logger = ConsoleLogger(stderr, log_level; show_limited=false)
#global_logger(logger)
