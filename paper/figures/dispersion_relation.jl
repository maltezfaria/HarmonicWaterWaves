using Plots

gr()

xx = 0:0.01:10

d = 1

# fig = plot(xx,tanh.(d*xx),label="tanh(x)",lw=2,xlabel="k")
# plot!(fig,xx,1 ./ xx, label="1/x",ylim=(0,2),lw=2)
# savefig(fig,"propagative_mode.png")

labeled = false
fig = plot(xlabel="k",lw=2,ylim=(-4,4),xlim=(0,10))
for n in 0:100
    xx = n*(π/(2d))+0.001:0.01:(n+1)*π/(2d)-0.001
    fig = plot!(xx,tan.(d*xx),label= labeled ? "" : "tan(γd)",lc=:blue)
    plot!(fig,xx,-1 ./ xx, label= labeled ? "" : "-1/(γd)",lw=2,lc=:red)
    labeled=true
end
fig
savefig(fig,"dissipative_modes.png")
