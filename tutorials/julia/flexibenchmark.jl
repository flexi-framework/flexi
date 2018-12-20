using BenchmarkTools
using MPI

const flexilib  =  "/home/iagbole/Codeentwicklung/flexi/build/lib/libflexi.so"
macro fl(m::String,f::String)
  (String("__"*lowercase(m)*"_MOD_"*lowercase(f)) , flexilib )
end

replARGS  = ("parameter_convtest_flexi.ini",)
flexiargs =  collect(replARGS)
nargs     = length(flexiargs)
strout    = ""
for i = 1:nargs
  global strout = String(strout * flexiargs[i] * (" "^(255-length(flexiargs[i]))))
end

# Initialize the MPI stuff
MPI.Finalized()   && error("MPI_Init no longer possible after MPI_Finalize has been called once! Please restart Julia!")
MPI.Initialized() || MPI.Init()

n = 3
nelems  = 64
nsides  = 192
u       = rand(Float64,5,n+1,n+1,n+1,nelems)
ut      = rand(Float64,5,n+1,n+1,n+1,nelems)
u1prim  = rand(Float64,6,n+1,n+1,n+1,nelems)
u2prim  = rand(Float64,6,n+1,n+1,n+1,nelems)
g1x     = rand(Float64,6,n+1,n+1,n+1,nelems)
g1y     = rand(Float64,6,n+1,n+1,n+1,nelems)
g1z     = rand(Float64,6,n+1,n+1,n+1,nelems)
v       = rand(n+1,n+1)

umaster = zeros(Float64,5,n+1,n+1,nsides)
uslave  = zeros(Float64,5,n+1,n+1,nsides)
lminus  = rand( Float64,n+1)
lplus   = rand( Float64,n+1)

function finit()
  global nargs,strout
  ccall( @fl("MOD_Flexi","InitFlexi") , Cvoid, (Ref{Int32},Ptr{UInt8},Ref{Int32}), Int32(nargs), strout, MPI.MPI_COMM_WORLD)
end
function ffinal()
  ccall( @fl("MOD_Flexi","FinalizeFlexi") ,Cvoid, () )
end

function liftcons(dir,uprim,gx,niter)
  dir2  =  Int32(dir)
  ccall( @fl("MOD_Lifting_Volint","Lifting_VolInt_conservative")  ,Cvoid,
             (Ref{Int32},Ptr{Array{Float64,5}},Ptr{Array{Float64,5}}),
             dir2,uprim,gx)
end
function liftnoncons(uprim,gx,gy,gz,niter)
  ccall( @fl("MOD_Lifting_Volint","Lifting_VolInt_Nonconservative")  ,Cvoid,
             (Ptr{Array{Float64,5}},Ptr{Array{Float64,5}},Ptr{Array{Float64,5}},Ptr{Array{Float64,5}}),
             uprim,gx,gy,gz)
end
function dgvolint(ut)
  ccall( @fl("MOD_Volint","VolInt_weakform")  ,Cvoid, (Ptr{Array{Float64,5}},),ut)
end
function prolongtoface(n,u,umaster,uslave,lminus,lplus)
  ccall( @fl("MOD_ProlongToFaceCons","ProlongToFace")  ,Cvoid,
  (Ref{Int32},Ptr{Array{Float64,5}},Ptr{Array{Float64,4}},Ptr{Array{Float64,4}},
  Ptr{Array{Float64,1}},Ptr{Array{Float64,1}},Ref{Bool}),
  n,u,umaster,uslave, lminus, lplus, false)
end

finit()
@btime dgvolint(ut)
@btime prolongtoface(n,u,umaster,uslave,lminus,lplus)
ffinal()
