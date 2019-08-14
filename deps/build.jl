using Conda

ENV["PYTHON"] = ""  # Configure PyCall.jl to use Conda.jl's Python

# We need install a specific version of llvmlite which supports the same minor version of
# llvm used with Julia. https://github.com/numba/llvmlite#compatibility
const LLVM_VERSION = Base.libllvm_version
const LLVMLITE_VERSION = if LLVM_VERSION > v"7.0.0"
    v"0.29.0"
elseif LLVM_VERSION >= v"6.0.0"
    v"0.26.0"
elseif LLVM_VERSION >= v"5.0.0"
    v"0.22.0"
elseif LLVM_VERSION >= v"4.0.0"
    v"0.20.0"
elseif LLVM_VERSION >= v"3.9.0"
    v"0.17.0"
elseif LLVM_VERSION >= v"3.8.0"
    v"0.15.0"
elseif LLVM_VERSION >= v"3.7.0"
    v"0.12.1"
elseif LLVM_VERSION >= v"3.6.0"
    v"0.8.0"
elseif LLVM_VERSION >= v"3.5.0"
    v"0.5.1"
else
    error("Can't find a version of llvmlite that supports $LLVM_VERSION")
end

Conda.update()
Conda.add_channel("conda-forge")
Conda.add("numpy==1.17.0")
Conda.add_channel("numba")  # channel containing llvmlite
Conda.add("llvmlite==$LLVMLITE_VERSION")  # https://github.com/numba/llvmlite
Conda.add_channel("invenia")
Conda.add("pandapower==2.1.0")
@info "Verifying pandapower install..."
python_bin = joinpath(Conda.PYTHONDIR, "python")
run(`$python_bin -c "import pandapower, pandapower.networks, pandapower.topology, pandapower.plotting, pandapower.converter, pandapower.estimation"`)
run(`$python_bin -c "from pandapower.networks import case9; case9()"`)