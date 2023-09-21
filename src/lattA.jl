module lattA

import juobs, ADerrors

include("types.jl")
export EnsInfo

include("reader.jl")
export get_corr_wil

include("obs.jl")
export get_m

include("tools.jl")
export model_av

end
