println("Script running.")
start_time = time()

import Pkg
Pkg.activate(dirname(@__DIR__))

using Dates

include("../src/includes.jl")


#run_2d_simulation(EchoConfig(alpha=3000,beta=0,Ny=1,Nz=128,Nd=256,Nt=512,pulses=default_weak_2d_pulses())) #attenuation 1d
#run_2d_simulation(EchoConfig(alpha=3000,beta=0,Ny=1,Nz=128,Nd=256,Nt=512,pulses=default_soliton_2d_pulses())) #sit 1d

#run_2d_simulation(EchoConfig(alpha=3000,beta=0,Ny=1,Nz=128,Nd=256,Nt=512,pulses=default_4pi_2d_pulses())) #1d 4pi
#run_2d_simulation(EchoConfig(alpha=3000,beta=0,Ny=1,Nz=128,Nd=256,Nt=512,pulses=default_6pi_2d_pulses())) #1d 6pi

#run_2d_simulation(EchoConfig(alpha=3000,beta=0,Ny=1,Nz=128,Nd=128,Nt=256,pulses=default_echo_2d_pulses())) #1d echo

#run_2d_simulation(EchoConfig(alpha=0,beta=1,Ny=128,Nz=64,Nd=128,Nt=256,pulses=default_weak_2d_pulses())) #diffraction 2d

run_2d_simulation(EchoConfig(alpha=3000,beta=0,Ny=128,Nz=128,Nd=256,Nt=512,pulses=default_soliton_2d_pulses())) #2d SIT
run_2d_simulation(EchoConfig(alpha=3000,beta=1,Ny=128,Nz=128,Nd=256,Nt=512,pulses=default_soliton_2d_pulses())) #2d SIT diffrac

run_2d_simulation(EchoConfig(alpha=3000,beta=0,Ny=128,Nz=128,Nd=128*2,Nt=256*2,pulses=default_echo_2d_pulses())) #2d echo
run_2d_simulation(EchoConfig(alpha=3000,beta=1,Ny=128,Nz=128,Nd=128*2,Nt=256*2,pulses=default_echo_2d_pulses())) #2d echo diffrac