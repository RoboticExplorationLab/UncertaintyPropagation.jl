# Loading File for Running Structures and Experiments

using PolyChaos
using LinearAlgebra
using Convex, SCS
using Distributions
using Random
using ForwardDiff
using LinearAlgebra
using Mosek, MosekTools
using Parameters
using Plots
using Statistics
using KernelDensity
using TrajectoryOptimization

const TO = TrajectoryOptimization
