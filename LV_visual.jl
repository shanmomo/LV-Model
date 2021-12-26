include("LV_model.jl")

## 整理拷贝画图数据
range_plt = 1:250 # plot range
tdata = xsol.t[range_plt] # time axis
xdata = xsol.u[range_plt] # mass datas

## 动画布局设置
# generate scatter locations using GraphPlot Pkg
locx,locy = spring_layout(g)
# 坐标中心
mx = sum(locx)/Nv ; my = sum(locx)/Nv;
# 计算各点的相对坐标和距离
related_loc = [[x-mx,y-my] for (x,y) in zip(locx,locy)]
related_dis = [	hypot(loc...) for loc in related_loc]
# 将最远离中心的若干点加权移动到中心附近
nodes_far = sortperm(related_dis)[end-8:end]
@. locx[nodes_far] = ( 0.2 * locx[nodes_far] + 0.8 * mx)
@. locy[nodes_far] = ( 0.2 * locy[nodes_far] + 0.8 * my)

xlimlv = ( minimum(locx) - 0.1 , maximum(locx) + 0.1 )
ylimlv = ( minimum(locy) - 0.10 , maximum(locy) + 0.15 )

gfadj=g.fadjlist # Forward adjacency list

## assign edge width and transparency according to the flow strength.
# w = bias_w + ln(1+f)  ;  α = bias_α + ln(1+f)
function plotedge!(i,j,flow)
	plot!(locx[[i,j]],locy[[i,j]],
		lw = 0.2 + log1p(abs(flow)),
		linealpha = 0.2 + log1p(abs(flow)),
		lc=:white,labels=nothing,)
end

## assign marker size and transparency according to the vertices' mass.
function plotvertecies!(xs)
	scatter!(locx,locy,
		m = (	3.0 .+ sqrt.(xs) .* 5,
				0.8,
			:blues, Plots.stroke(0)),
		labels = nothing )
end

gr() ; theme(:lime)
ani = @animate for tt in range_plt
	plot(xlims = xlimlv, ylims = ylimlv,
		axis = false,
		ticks = false,
		size = [1920,1080])
	xs = xdata[tt]
	for i = 1:Nv
		for j in gfadj[i]
			plotedge!(i,j,xs[i]*xs[j])
		end
	end
	plotvertecies!(xs)
	annotate!(0.85 * xlimlv[1] + 0.15 * xlimlv[2], 0.03*ylimlv[1] + 0.97*ylimlv[2],
		annotationcolor = :white,
		annotationfontsize = 30,
		latexstring("\$ N=$(Nv)\\, ,\\,\\bar{k}_i=\\bar{k}_o=$(kio_ave) \$"))
	annotate!(0.09 * xlimlv[1] + 0.91 * xlimlv[2], 0.03*ylimlv[1] + 0.97*ylimlv[2],
		annotationcolor = :white,
		annotationfontsize = 30,
		latexstring("\$ t= \$"))
	annotate!(0.05 * xlimlv[1] + 0.95 * xlimlv[2], 0.03*ylimlv[1] + 0.97*ylimlv[2],
		annotationcolor = :white,
		annotationfontsize = 30,
		latexstring("\$ {$(round(xsol.t[tt],sigdigits=2))} \$"))
end

gif(ani, "N300K2.gif", fps = 10)
