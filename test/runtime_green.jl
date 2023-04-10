import HarmonicWaterWaves as WW
using BenchmarkTools
using SpecialFunctions

x = rand()
y = rand()
z = x-im*y
@btime expinti($z)

@btime log($x)
