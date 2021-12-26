##一些库
using LinearAlgebra
using SparseArrays
using DifferentialEquations.OrdinaryDiffEq
using Graphs
using Plots, GraphPlot, LaTeXStrings

## Graph Generation
Nv, kio_ave, LVseed = (300, 4, 6) # 节点个数, 平均入/出度, 随机数种子
Ne = kio_ave * Nv

# 生成对应的有向图和邻接矩阵（有向）
function random_directed_graph(nv,ne,s)
	g = SimpleDiGraph(nv,ne; seed=s )
	A = adjacency_matrix(g) -transpose(adjacency_matrix(g))
	return g,A
end

g,A = random_directed_graph(Nv,Ne,LVseed)

# 绘制质量流系数对应的稀疏矩阵，并利用稀疏矩阵加速
theme(:juno)
adjMat_plt = heatmap(  reverse(Matrix(A[1:60,1:60]),dims = 1),
	axis = false, xlims=(0,60), ylims=(0,60),
	aspect_ratio=:equal, c = :berlin)
savefig(adjMat_plt,"A_spy.pdf")

## ODE Problem

# Parameters and initial conditions
begin
	x0 = ones(Nv)
	xt = deepcopy(x0)
	u0 = log.(x0)
	M0 = sum(x0)
	tspan = (0.0,30.0)
	ϵ = 0.01
end

function LV_model!(du,u,p,t) # THE Lotka-Volterra Problem
	xt .= exp.(u)
	du .= (A * xt)
	M = sum(xt)
	du .-= (ϵ*(M-M0)) # Mass Conservation Correction
end

# Solve with DifferentialEquations.jl
prob = ODEProblem(LV_model!,u0,tspan)
usol = solve(prob,TsitPap8();
 		dt=0.1,adaptive=false)

ΔM = sum( exp.(usol.u[end]) )-M0 # check Mass Conservation

## 将对数结果转化为原本的质量结果
function usol2xsol(usol)
    xsol = deepcopy(usol)
    @inbounds begin
        for ii=1:length(usol.t)
            xsol.u[ii] = exp.(usol.u[ii])
            xsol.k[ii][1] .*= xsol.u[ii]
        end
		for ii=2:length(usol.t)
	        xsol.k[ii][2] .*= xsol.u[ii]
		end
    end
    return xsol
end

xsol = usol2xsol(usol)

## 绘制质量随时间演化的曲线图
theme(:dark)
massevo_plt=plot(0:0.01:25,[xsol(0:0.01:25)[i,:] for i =1:1:Nv],
	labels = nothing,
	xlabel = L"t",
	ylabel = L"x_i")
savefig(massevo_plt,"massevo.pdf")
