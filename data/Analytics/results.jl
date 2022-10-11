    using LaTeXStrings
    using Plots
    using JLD
    using DelimitedFiles
    #Reading variables
    d=load("results/variables.jld");
    N_h2=d["N_h2"];E_h2=d["E_h2"];List_h2=d["List_h2"]
    N_h2_6=d["N_h2_6"];E_h2_6=d["E_h2_6"]; List_h2_6=d["List_h2_6"];
    N_h2O=d["N_h2O"];E_h2O=d["E_h2O"];List_h2O=d["List_h2O"];    #plot1 Dihydrogen in STO-3g basis
    #E vs N
    plot(N_h2,E_h2,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_2 in STO-3G basis set", legendfontsize=4)
    plot!([2],[-1.13728383],marker =:star8,markercolor=:red,markersize=4,label="Ground state energy for H_2 in  STO_3G basis")
    plot!(N_h2_6,E_h2_6,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_2 in 6-31G basis set")
    plot!([2],[-1.15167254],marker =:star8,markercolor=:green,markersize=4,label="Ground state energy for H_2 in  6-31G basis")
    plot!(N_h4,E_h4,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_4 in STO-3G basis set")
    plot!([4],[-1.47624026],marker =:star8,markercolor=:yellow,markersize=4,label="Ground state energy for H_4 in  STO-3G basis")
    plot!(N_h6,E_h6,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_6 in STO-3G basis set")
    p1=plot!([6],[2.21286898],marker =:star8,markercolor=:pink,markersize=4,label="Ground state energy for H_6 in  STO-3G basis")
    #plot!(N_h2O,E_h2O,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_2O in STO-3G basis set")
    #plot!([10],[-75.01053577546307],marker =:star8,markercolor=:blue,markersize=4,label="Ground state energy for H_2O in  STO-3G basis")
    xlabel!("N_\\mu")
    ylabel!("E_\\mu")
    Plots.savefig("Hydrogen_chain.pdf")

    plot(N_h2O[1:end-1],E_h2O[1:end-1],marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_2O in STO-3G basis set", legendfontsize=4)
    plot!([10],[-75.01053577546307],marker =:star8,markercolor=:blue,markersize=4,label="Ground state energy for H_2O in  STO-3G basis")
    xlabel!("N_\\mu")
    ylabel!("E_\\mu")
    Plots.savefig("Water.pdf")

    #E vs mu
    plot(List_h2,abs.(E_h2 .+1.13728383)./abs(-1.13728383),linewidth=0.5,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_2 in STO-3G basis set",yaxis=:log, legendfontsize=3)
    plot!(List_h2_6,abs.(E_h2_6.+1.15167254)./abs(-1.15167254),linewidth=0.5,marker =:circle,markersize=2.5,label ="H_2 in 6-31G basis set")
    plot!(List_h4,abs.(E_h4 .+1.47624026)./abs(-1.47624026),linewidth=0.5,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_4 in STO-3G basis set")
    plot!(List_h6,abs.(E_h6 .-2.21286898)./abs(2.21286898),linewidth=0.5,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_6 in STO-3G basis set")
    plot!(List_h2O[1:end-1],abs.(E_h2O[1:end-1] .+75.01053577546307)./abs(75.01053577546307),linewidth=0.5,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_20 in STO-3G basis set")

    xlabel!("\\mu")
    ylabel!("\\Delta E_\\mu")
    Plots.savefig("energyvsmu.pdf")
    #N vs mu
    plot(List_h2,N_h2,linewidth=0.5,marker =:circle,markersize=2.5,title = "Particle number conservation by adding chemical potential",label ="H_2 in STO-3G basis set", legendfontsize=4)
    plot!(List_h2_6,N_h2_6,linewidth=0.5,marker =:circle,markersize=2.5,label ="H_2 in 6-31G basis set")
    plot!(List_h4,N_h4,linewidth=0.5,marker =:circle,markersize=2.5,label ="H_4 in STO-3G basis set")
    plot!(List_h6,N_h6,linewidth=0.5,marker =:circle,markersize=2.5,label ="H_4 in STO-3G basis set")
    plot!(List_h2O[1:end-1],N_h2O[1:end-1],linewidth=0.5,marker =:circle,markersize=2.5,label ="H_4 in STO-3G basis set")
    xlabel!("\\mu")
    ylabel!("N_\\mu")
    Plots.savefig("Nvsmu.pdf")
##Energy vs Discarded weight with diff tols
#LIH
tol_MPS_all=[1e-2,1e-3,1e-4,1e-5,1e-6,1e-8,1e-10,1e-12,1e-14,1e-16]
E_DW_lih= readdlm("data/E_discarded_lih.output");
plot(E_DW_lih[:,3],abs.(E_DW_lih[:,1].-(-7.88240342)),yaxis=:log,xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="LiH molecule")
E_DW_h2o= readdlm("data/E_discarded_h2o.output");
plot!(E_DW_h2o[:,3],E_DW_h2o[:,1].-(-75.012578266523052),xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="H2O molecule")
E_DW_Be2= readdlm("data/E_discarded_Be21_1_6.output");
plot(E_DW_Be2[:,3],E_DW_Be2[:,1].-(-28.80434531),xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="Be2 molecule")
N2_dW=[0.2545833691599805 0.22990292394715856 0.02807572016952029 0.0032080072409191407 0.0010563765432282345 ]
plot(N2_dW',N2_dmrg'[1:end-1].-(-107.65277153),yaxis=:log,xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="N2 molecule")
xlabel!("Discarded weight")
ylabel!("E_dmrg-E_FCI")
Plots.savefig("Discarded_weight_N2.pdf")

##Disso curves with different tols for LIH
E_lih_true= readdlm("data/LiH/LiH_tols/E_tol_lihh.output");
energy_fci_lih=[-7.06946645    -7.58310908  -7.78983261 -7.86363546 -7.88252884 -7.87663413 -7.86026583 -7.84097366 -7.80861966 -7.79195886 -7.78561546 -7.78343218]# -7.32719572 -7.34578399]
energy_scf_lih=[-7.0478926046  -7.56438711715130 -7.7727824926 -7.8466364488 -7.86338213302325 -7.85304056808981 -7.82972034878235 -7.80059541031344 -7.7386724024 -7.68328537606259 -7.64066210706671 -7.61089630899898]#-7.3667471926 -7.2624180810]
Bond_lengths_lih=[0.506  0.756 1.011 1.261 1.511 1.761 2.011 2.261 2.761 3.261 3.761 4.261 5.261 7.261 10.261]
values_zeros=findall(x->x==0.0, E_lih_true[:,3]);
tol_MPO=E_lih_true[values_zeros,4]
tol_MPS=E_lih_true[values_zeros,6]
r_true=E_lih_true[values_zeros,1]
A_lih=zeros(length(tol_MPO),length(Bond_lengths_lih));

for i=1:length(Bond_lengths_lih)
    E=readdlm("data/LiH/LiH_tols/E_tol_lihh_$(i).output");
    A_lih[:,i]=E[values_zeros,1]
end
plot(Bond_lengths_lih'[1:12],abs.(energy_fci_lih' .-(-7.88252884))./abs(-7.88252884),linewidth=1,marker=:circle,label="Full-CI",title="LiH dissociation curve with different ranks")
#ylims!((-6.5,-8))
plot!(Bond_lengths_lih'[1:12],abs.(energy_scf_lih'.-(-7.88252884))./abs(-7.88252884),linewidth=1,marker =:cross,markersize=2.5,label="SCF")
for i=36:42
    plot!(Bond_lengths_lih',abs.(A_lih[i,:] .-(-7.88252884))./abs(-7.88252884),linewidth=1,marker =:circle,markersize=2.5,legend =true,label="tol_MPS=$(tol_MPS[i]),tol_MPO=$(tol_MPO[i])",legendfontsize=5)
end
xlabel!("Li-H bond length(Å)")
ylabel!("E")
Plots.savefig("Dissos_curve_lih_tols_v3.pdf")

#EvsNmu with diff tol
plot([4],[-7.88240342],markercolor=:green,marker=:cross,markersize=10,label="Ground state energy for LiH ",title="LiH molecule, STO-3G basis,L=6")£
for i in [1e-8]
    for j=1:length(tol_MPS[1:7])
        id1=findall(x->x==i, E_lih_true[:,4]);
        id2=findall(x->x==tol_MPS[j], E_lih_true[:,6]);
        id=intersect(id1,id2)
        #plot!(E_lih_true[id,2],E_lih_true[id,1],linewidth=1,marker =:circle,markersize=2.5,label="$(tol_MPS[j]),$(i)",legend=false)
        plot!(E_lih_true[id,2],E_lih_true[id,1],linewidth=1,marker =:circle,markersize=2.5,label="tol_MPS:$(tol_MPS[j]),tol_MPO:$(i)",legend=false)
    end
end
xlabel!("N_\\mu")
ylabel!("E_\\mu")
Plots.savefig("E_mu_vs_mu_v2.pdf")
#E vs N_mu with diff r
id1=findall(x->x==1e-4, E_lih_true[:,4]);
id2=findall(x->x==1e-8, E_lih_true[:,6]);
id=intersect(id1,id2)
plot(E_lih_true[id,2],E_lih_true[id,1],linewidth=1,marker =:circle,markersize=2.5,label="r=$(Bond_lengths_lih[1])",title="LiH: E_mu vs N_mu ")
for i=2:length(Bond_lengths_lih)
        E=readdlm("data/LiH/LiH_tols/E_tol_lihh_$(i).output");
        id1=findall(x->x==1e-4, E[:,4]);
        id2=findall(x->x==1e-8, E[:,6]);
        id=intersect(id1,id2)
        #plot!(E_lih_true[id,2],E_lih_true[id,1],linewidth=1,marker =:circle,markersize=2.5,label="$(tol_MPS[j]),$(i)",legend=false)
        plot!(E[id,2],E[id,1],linewidth=1,marker =:circle,markersize=2.5,label="r=$(Bond_lengths_lih[i])",legend=:outertopright,legendfontsize=5)
end
xlabel!("N_\\mu")
ylabel!("E_\\mu")
Plots.savefig("E_mu_vs_N_mu_r_badmpo.pdf")

##Pres
#Relative error vs MPS_tol
#MPO 1e-8
tol_MPS=[1e-2,1e-3,1e-4,1e-5,1e-6,1e-8,1e-10,1e-12,1e-14,1e-16]
E_H2=readdlm("Outputs/h2_MPS.output")
Error_H2=E_H2[:,3]
plot(tol_MPS,Error_H2,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H2 molecule in STO-3G basis L=2, tol_MPO=1e-2",legendfontsize=5,legend=:topleft)
E_H4=readdlm("Outputs/h4_MPS.output")
Error_H4=E_H4[:,3]
plot!(tol_MPS,Error_H4,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H4 molecule in STO-3G basis L=4, tol_MPO=1e-8")
E_H6=readdlm("Outputs/h6_MPS.output")
Error_H6=E_H6[:,3]
plot!(tol_MPS,Error_H6,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H6 molecule in STO-3G basis L=6, tol_MPO=1e-8")
E_LiH=readdlm("Outputs/Lih_MPS.output")
Error_LiH=E_LiH[:,3]
plot!(tol_MPS,Error_LiH,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="LiH molecule in STO-3G basis L=6, tol_MPO=1e-8")
E_H2O=readdlm("Outputs/h2o_MPS.output")
Error_H2O=E_H2O[:,3]
plot!(tol_MPS,Error_H2O,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H2O molecule in STO-3G basis L=7, tol_MPO=1e-8")
tol_MPS=[1e-2,1e-3,1e-4,1e-5,1e-6]
Error_Be2=[0.008965665762955877 0.008965665845066932 0.0030723099369860517 0.0030115295586261525 0.0016450423703040942]
plot!(tol_MPS,Error_Be2',yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="Be2 molecule in STO-3G basis L=10,tol_MPO=1e-6")
Error_N2=[0.0028456578774893335 0.0020942851383632666 0.00199425168668599 0.0012269138408158634 0.0008296998268917912]
plot!(tol_MPS,Error_N2',yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="N2 molecule in STO-3G basis L=10,tol_MPO=1e-6")
xlabel!("MPS tolerance")
ylabel!("\\DeltaE")
Plots.savefig("E_vs_MPS.pdf")


#Relative error vs MPS_rank
E_H2=readdlm("h2_MPS.output")
rank_H2=E_H2[:,4]
plot(rank_H2,Error_H2,yaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H2 molecule in STO-3G basis,L=2",)
E_H4=readdlm("h4_MPS.output")
rank_H4=E_H4[:,4]
plot!(rank_H4,Error_H4,yaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H4 molecule in STO-3G basis,L=4")
E_H6=readdlm("h6_MPS.output")
rank_H6=E_H6[:,4]
plot!(rank_H6,Error_H6,yaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H6 molecule in STO-3G basis,L=6")
E_LiH=readdlm("Lih_MPS.output")
rank_LiH=E_LiH[:,4]
plot!(rank_LiH,Error_LiH,yaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="LiH molecule in STO-3G basis,L=6")
xlabel!("Maximum TT-rank (MPS)")
ylabel!("\\DeltaE")
Plots.savefig("E_ranks_MPS.pdf")

#Relat error vs MPO_tol
tol_MPO=[1e-1,1e-2,1e-3,1e-4,1e-6,1e-8];
E_H2=readdlm("h2_MPO.output")
Error_H2=E_H2[:,3]
plot(tol_MPO,Error_H2,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H2 molecule in STO-3G basis,L=2,tol_MPS=1e-12",legendfontsize=4,legend=:topleft)
E_H4=readdlm("h4_MPO.output")
Error_H4=E_H4[:,3]
plot!(tol_MPO,Error_H4,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H4 molecule in STO-3G basis,L=4,tol_MPS=1e-12")
E_lih=readdlm("lih_MPO.output")
Error_lih=E_lih[:,3]
plot!(tol_MPO,Error_lih,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="LiH molecule in STO-3G basis,L=6,tol_MPS=1e-12")
E_H6=readdlm("H6_MPO.output")
Error_H6=E_H6[:,3]
plot!(tol_MPO,Error_H6,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H6 molecule in STO-3G basis,L=6,tol_MPS=1e-12")
E_H2O=readdlm("h2o_MPO2.output")
Error_H2O=E_H2O[:,3]
plot!(tol_MPO,Error_H2O,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H2O molecule in STO-3G basis,L=7,tol_MPS=1e-12")
tol_MPO=[1e-1,1e-2,1e-3,1e-4,1e-6]
Error_Be2=[0.37040775014528343 0.3174050337055484 0.16020895693042714 0.056347322834357914 0.008180762508995193]
plot!(tol_MPO,Error_Be2',yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="Be2 molecule in STO-3G basis,L=10,tol_MPS=1e-6")
xlabel!("MPO tolerance")
ylabel!("\\DeltaE")
Plots.savefig("E_tols_MPO.pdf")

#time vs MPO_tol
Time_H2=E_H2[:,5]
plot(tol_MPO,Time_H2,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H2 molecule in STO-3G basis",legend=:outertopright)
Time_H4=E_H4[:,5]
plot!(tol_MPO,Time_H4,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H4 molecule in STO-3G basis",legend=:outertopright)
Time_lih=E_lih[:,5]
plot!(tol_lih,Time_lih,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="LiH molecule in STO-3G basis",legend=:outertopright)
Time_H6=E_H6[:,5]
plot!(tol_H6,Time_H6,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H6 molecule in STO-3G basis",legend=:outertopright)

xlabel!("MPO tolerance")
ylabel!("Time(s)")

#Relat error vs MPO_rank
rank_H2=E_H2[:,4]
plot(rank_H2,Error_H2,yaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H2 molecule in STO-3G basis,L=2,tol_MPS=1e-12",legendfontsize=3)
rank_H4=E_H4[:,4]
plot!(rank_H4,Error_H4,yaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H4 molecule in STO-3G basis,L=4,tol_MPS=1e-12")
rank_LiH=E_lih[:,4]
plot!(rank_LiH,Error_H6,yaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="LiH molecule in STO-3G basis,L=6,tol_MPS=1e-12")
rank_H6=E_H6[:,4]
plot!(rank_H6,Error_H6,yaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H6 molecule in STO-3G basis,L=6,tol_MPS=1e-12")
rank_H2O=E_H2O[:,4]
plot!(rank_H2O,Error_H2O,yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="H2O molecule in STO-3G basis,L=7,tol_MPS=1e-12")
rank_Be2=[1.0,2.0,10.0,80.0,98.0]
plot!(rank_Be2,Error_Be2',yaxis=:log,xaxis=:log,linewidth=1,marker =:circle,markersize=2.5,label="Be2 molecule in STO-3G basis,L=10,tol_MPS=1e-12")

xlabel!("Maximum MPO rank")
ylabel!("\\DeltaE")
Plots.savefig("E_ranks_MPO.pdf")
##Discarded_weight_C2
#C2
Energy_C2=[-74.69021094]
tol_MPS_all=[1e-2,1e-3,1e-4,1e-5,1e-6]
C2_dmrg=[-74.28782072631164 -74.46888116687651 -74.50379520732703  -74.58147682423747 -74.61514515877109]
C2_dW=[0.46246039034552827 0.04324855060240912  0.011900971579808337 0.006266968497053382 0.0013695896361487713]
plot(C2_dW',C2_dmrg'[1:end].-(-74.69021094),xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="N2 molecule")

Be2_dmrg=[-28.48501754403231  -28.546095167674288  -28.708385916266717 -28.73478902468306 -28.72760849965315 -28.80652296360926 ]
Be2_dW=[0.14924605460794588  0.10318905523941185 0.01813624573474323 0.0033344092517695466 0.0013938419818108842]
plot(Be2_dW',Be2_dmrg'[1:end-1].-(-28.80434531),xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="N2 molecule")

N2_dmrg=[-107.34671950621627 -107.42799237910548 -107.43878510711977 -107.5223791511034 -107.55843364350218 -107.59045097357283]
N2_dW=[0.2545833691599805 0.22990292394715856 0.02807572016952029 0.0032080072409191407 0.0010563765432282345 ]
plot(N2_dW',N2_dmrg'[1:end-1].-(-107.65277153),xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="N2 molecule")


##Totalspin
#Lih
Bond_lengths_lih=[0.506  0.756 1.011 1.261 1.511 1.761 2.011 2.261 2.761 3.261 3.761 4.261 5.261 7.261 10.261]
Energy3=[-7.069466454629781 -7.58310908456387 -7.789705342926 -7.863635333395484 -7.8825286901531495 -7.876634127326262 -7.860265796202383 -7.840973635325725 -7.808618390149074 -7.7919561911536075 -7.7856146936559405 -7.783432184429271 -7.782501479366335 -7.7824189319850365 -7.782418427943522]
Lih=readdlm("Weight_lih.output")
plot(Bond_lengths_lih',[Lih' Energy3'],linewidth=1,marker =:circle,markersize=2.5,label="LiH molecule",title=["Total spin operator vs bond lengths" "Energy vs bond lengths"],xlabel="Bond length(Å)",ylabel=["<S^2>" "E"],layout = 2)
Plots.savefig("total_spin_lih.pdf")

Bond_lengths_h2o=[0.1144  0.3 0.4 0.5144 0.6  0.7 0.8 0.9 1 1.5144  1.9856 2.4856 2.9856 3.4856 4.4856 6.4856 8.4856 9.4856]
Energy4=[-70.24210728900874  -73.13666784758006 -73.62882057047491 -74.00186046630948 -74.21262183547283 -74.41056142037499 -74.56866467476631 -74.69464955617035 -74.79390922848496 -75.0125519782931 -74.98206814432595 -74.89092352158461 -74.81716018660346 -74.77370073485763 -74.74357726913104 -74.73743042643378  -74.73732073139882 -74.73731395992603]
4.5455611157808615
H2O=readdlm("W_H2O.output")
plot(Bond_lengths_h2o',[H2O' Energy4'],linewidth=1,marker =:circle,markersize=2.5,label="H2O molecule",layout = 2,title=["Total spin operator vs bond lengths" "Energy vs bond lengths"],xlabel="Bond length(Å)",ylabel=["<S^2>" "E"])
Plots.savefig("total_spin_h2o.pdf")

Bond_lengths_h2=[0.25 0.5 0.6 0.7414 0.8 0.9914 1.2414 1.4914 1.7414 1.9914 2.4914 2.9914 3.4914 3.9914]#,4.9914,5.9914]
energy_scf_h2=[-0.305000531255344,-1.04299627382273,-1.10112824196162,-1.11668438720084,-1.11085039771967,-1.06850935949785,-0.991868551365561,-0.913453211828911,-0.843548158039365,-0.785568450214639,-0.704003871693402, -0.656647628792305,-0.630154538852619,-0.615068165498246,-0.599121806441454,-0.590035079802113]
energy_fci_h2=[-0.31226990,-1.05515979,-1.11628601,-1.13727017,-1.13414767,-1.10295561,-1.04764544,-0.99953646,-0.96716840,-0.94907006,-0.93614316,-0.93364737,-0.93323074,-0.93317166,-0.93316377,-0.93316370]
#tol_MPO=1e-2,bond=8
energy_dmrg_h2=[-0.3122699002638994 -1.0551597938192845 -1.1162860066241516 -1.1372701748287717 -1.1341476669708221 -1.102955605430223 -1.0476454361534762 -0.9995364595556123 -0.9671684029721406 -0.9490700634305685 -0.9361431637435632 -0.9336473680926343 -0.9332307372506189   -0.93317166]#,-0.9331637664803448,-0.9331637009473837]
H2=[3.964783258150895e-30 7.299213014497305e-31 9.391427491377852e-31 7.4131133674860774e-31 5.775435613754226e-30 1.7461328558872774e-30 3.1693511356094015e-30 4.537034821977264e-31 3.680035920706221e-32 3.67867196562087e-31 7.2646336073914615e-31 5.053613076487893e-31 4.384877365115662e-31 1.3590771332753966e-27]
plot(Bond_lengths_h2',[H2' energy_dmrg_h2'],linewidth=1,marker =:circle,markersize=2.5,label="H2 molecule",layout = 2,title=["Total spin operator vs bond lengths" "Energy vs bond lengths"],xlabel="Bond length(Å)",ylabel=["<S^2>" "E"])
Plots.savefig("total_spin_h2.pdf")
##Change order
#1e-3
Energy2_default=[-7.0646008510134175	-7.444378369621541 -7.772782492586138 -7.860317421169011 -7.878426541414133 -7.735125448313238 -7.829720348782118 -7.776881554968909 -7.77850539867692 -7.780323878926449 -7.781445348623647 -7.781901057983667 -7.782082486855148 -7.781037160720833 -7.782099703266962]
Energy_o1_327=[-7.060965789202651 -7.57967930070199 -7.772782492588609 -7.846636448813843 -7.877800350233734 -7.745673305250666 -7.715875401541436 -7.776881559497467 -7.778505395646366 -7.780323877926042 -7.781445347814225 -7.781901063436323 -7.782082505358305 -7.781037122222708 -7.782099702207177]
Energy_o2_334=[-7.064035354846116 -7.580808635691555 -7.773939222167618 -7.846636448841723 -7.877816400122822 -7.853727487032359 -7.775769751858088 -7.77699592170214 -7.778579255324616 -7.780365230323522 -7.781467067674727 -7.7819120229843275 -7.782085177587835 -7.78207161703539 -7.782099702209134]
Energy_o2_322=[-6.985816585358591 -7.580891131508159 -7.77278249258898 -7.846636448843054 -7.863382133021259 -7.870351005240117 -7.775769740095136 -7.776995935971711 -7.778579254408537 -7.780365227069305 -7.77586135391634 -7.7819120282711145 -7.782085178730204 -7.696568738170398 -7.782099702208066]
Energy_o2_17=[-6.985779997579135 -7.444378345878289 -7.77278249259083 -7.860317421235563 -7.881073023544048 -7.85100578556124 -7.77564186286731 -7.7768815535264 -7.778505392215829 -7.780323877587781 -7.7814453484164146 -7.781901061582701 -7.696656436813727 -7.697605870844446 -7.696686530265212]
Energy_dmrg=[-7.069466454634794 -7.583109084568589 -7.789832609809361 -7.863635518658474 -7.882528846313271 -7.876634127330327 -7.860265796206015 -7.840973635329235 -7.808619663019514  -7.7919588556424335 -7.7856154573905005 -7.7834321844342 -7.782501479672469 -7.782418615404309 -7.782418427263477]
#Energy_o2_260=[-7.052388231169452 -7.564387117152868 -7.772782492589005 -7.846636448845388 -7.8633821330246985 -7.825715934104483 -7.829720348782107 -7.800595410308519 -7.745329107538503 -7.717889711772538 -7.737397393495885 -7.75975717985268 -7.7771635561561 -7.7820715356489165 -7.782099702209313]
Energy_o2_260=[-7.053008157977529 -7.566099958642702 -7.774012495285453 -7.846636448842325 -7.864403445519518 -7.840085082649438	-7.8297203487586025 -7.801775201916329 -7.746014459501261 -7.71788971766757 -7.737397400470715 -7.7597571768931894 -7.777163557214827 -7.781037296573158 -7.7820997020802025]
plot(Bond_lengths_lih',Energy_dmrg',linewidth=1,marker=:circle,label="Full-CI",title="LiH dissociation curve for tol-MPO=1e-8,tol-MPS=1e-3")
plot!(Bond_lengths_lih',Energy2_default',linewidth=1,marker=:circle,label="Canonical order")
plot!(Bond_lengths_lih',Energy_o1_327',linewidth=1,marker=:circle,label="Order 1")
plot!(Bond_lengths_lih',Energy_o2_334',linewidth=1,marker=:circle,label="Order 2")
plot!(Bond_lengths_lih',Energy_o2_322',linewidth=1,marker=:circle,label="Order 3")
#plot!(Bond_lengths_lih',Energy_o2_17',linewidth=1,marker=:circle,label="Order 4")
plot!(Bond_lengths_lih',Energy_o2_260',linewidth=1,marker=:circle,label="Order 5")

xlabel!("Bond length(Å)")
ylabel!("E")
Plots.savefig("order_dissos.pdf")
##
#dissoc LIH
#1e-2
Bond_lengths_lih=[0.506  0.756 1.011 1.261 1.511 1.761 2.011 2.261 2.761 3.261 3.761 4.261 5.261 7.261 10.261]
Energy1=[-6.979769540652963 -7.425286834991513 -7.597788510326041 -7.661739590632109 -7.6815133032188845 -7.6227486849241295 -7.669325887126533 -7.653335625755423 -7.703507018604413 -7.699328624904352 -7.697596214736528 -7.696952025515815 -7.696535111919673 -7.697605870490759 -7.696682200110796]
plot(Bond_lengths_lih',Energy1',linewidth=1,marker=:circle,label="Full-CI",title="LiH dissociation curve with different ranks")
#1e-3
Energy2=[-7.0646008510134175	-7.444378369621541 -7.772782492586138 -7.860317421169011 -7.878426541414133 -7.735125448313238 -7.829720348782118 -7.776881554968909 -7.77850539867692 -7.780323878926449 -7.781445348623647 -7.781901057983667 -7.782082486855148 -7.781037160720833 -7.782099703266962]
plot!(Bond_lengths_lih',Energy2',linewidth=1,marker=:circle,label="Full-CI",title="LiH dissociation curve with different ranks")
#1e-10
Energy3=[-7.069466454629781 -7.58310908456387 -7.789705342926 -7.863635333395484 -7.8825286901531495 -7.876634127326262 -7.860265796202383 -7.840973635325725 -7.808618390149074 -7.7919561911536075 -7.7856146936559405 -7.783432184429271 -7.782501479366335 -7.7824189319850365 -7.782418427943522]
plot!(Bond_lengths_lih',Energy3',linewidth=1,marker=:circle,label="Full-CI",title="LiH dissociation curve with different ranks")
energy_dmrg_h2o=[-70.24210921585846 -72.22793929 -73.1367043703099  -73.62882058793929 -74.0018605221657 -74.21262191183824 -74.41056146652681 -74.56866473302975 -74.69464964343852 -74.79390926152588 -75.01257819471171  -74.9820681971939 -74.89092363264108 -74.81716019045031 -74.77370076484597  -74.74357754783232 -74.73743147745849 -74.79390926164747 -74.73732225411209 -74.73731753664063]
##Compare 1-site dmrg with 2-sites dmrg
#LIH deltaE vs MPS tol
one_site=readdlm("Lih_data_1_site.output");
two_site=readdlm("Lih_data.output");

plot(one_site[:,4],abs.(one_site[:,1].- (-7.8824034242603345))./abs( -7.8824034242603345),xaxis=:log,linewidth=1,marker=:circle,label="1-site DMRG",title="Comparison between 1-site and 2-site DMRG")
plot!(two_site[:,4],abs.(two_site[:,1].- (-7.8824034242603345))./abs( -7.8824034242603345),linewidth=1,marker=:circle,label="2-site DMRG")
xlabel!("MPS accuracies")
ylabel!("\\DeltaE")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/energy.pdf")

#LIH max_rnk accuracies
plot(one_site[:,4],one_site[:,5],xaxis=:log,linewidth=1,marker=:circle,label="1-site DMRG",title="Comparison between 1-site and 2-site DMRG")
plot!(two_site[:,4],two_site[:,5],linewidth=1,marker=:circle,label="2-site DMRG")
xlabel!("MPS accuracies")
ylabel!("maximum ranks")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/maximum_ranks.pdf")
#LIH sweeps
plot(one_site[:,4],one_site[:,6],xaxis=:log,linewidth=1,marker=:circle,label="1-site DMRG",title="Comparison between 1-site and 2-site DMRG")
plot!(two_site[:,4],two_site[:,6],linewidth=1,marker=:circle,label="2-site DMRG")
xlabel!("MPS accuracies")
ylabel!("Sweeps")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/sweeps.pdf")
#Dissociatin curve fo 2-sites
Dissos=readdlm("Lih_data.output");
Dissos=reshape(Dissos,9,15,6);
Bond_lengths_lih=[0.506  0.756 1.011 1.261 1.511 1.761 2.011 2.261 2.761 3.261 3.761 4.261 5.261 7.261 10.261]
Energy3=[-7.069466454629781 -7.58310908456387 -7.789705342926 -7.863635333395484 -7.8825286901531495 -7.876634127326262 -7.860265796202383 -7.840973635325725 -7.808618390149074 -7.7919561911536075 -7.7856146936559405 -7.783432184429271 -7.782501479366335 -7.7824189319850365 -7.782418427943522]
plot(Bond_lengths_lih',Energy3',linewidth=1,marker=:circle,label="Full-CI",title="Dissociation curve for LiH with 2-sites DMRG")

for i=1:size(Dissos,1)
    plot!(Bond_lengths_lih',Dissos[i,:,1],linewidth=1,marker=:circle,label="tol_MPS=$(Dissos[i,1,4])")
end
xlabel!("Bond length(Å)")
ylabel!("E")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/dissocition_curve.pdf")

#MPO_tols
MPS_=1e-8
tol_MPO_all=[1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8]
Dissos_MPO=readdlm("Lih_data_MPO.output");
Dissos_MPO=reshape(Dissos_MPO,7,15,6);
Bond_lengths_lih=[0.506  0.756 1.011 1.261 1.511 1.761 2.011 2.261 2.761 3.261 3.761 4.261 5.261 7.261 10.261]
Energy3=[-7.069466454629781 -7.58310908456387 -7.789705342926 -7.863635333395484 -7.8825286901531495 -7.876634127326262 -7.860265796202383 -7.840973635325725 -7.808618390149074 -7.7919561911536075 -7.7856146936559405 -7.783432184429271 -7.782501479366335 -7.7824189319850365 -7.782418427943522]
plot(Bond_lengths_lih',Energy3',linewidth=1,marker=:star,label="Full-CI",title="Dissociation curve for LiH with 2-sites DMRG")

for i=1:size(Dissos_MPO,1)
    plot!(Bond_lengths_lih',Dissos_MPO[i,:,1],linewidth=1,marker=:circle,label="tol_MPO=$(Dissos_MPO[i,1,4])",legend=:outertopright)
end
xlabel!("Bond length(Å)")
ylabel!("E")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/dissocition_curve_MPO.pdf")

#E_mu vs N_mu for diff MPS tols
Curve=readdlm("examples/Lih_data_Nmu.output");
Curve=reshape(Curve,(9,25,:))
plot(Curve[1,:,2],Curve[1,:,1],linewidth=1,marker =:circle,markersize=2.5,label="tol_MPS:$(Curve[1,1,4])",title="2-sites DMRG LiH: E_mu vs N_mu ",legend=:topleft)
for i=2:9
    plot!(Curve[i,:,2],Curve[i,:,1],linewidth=1,marker =:circle,markersize=2.5,label="tol_MPS:$(Curve[i,1,4])")
end
xlabel!("N_\\mu")
ylabel!("E_\\mu")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/E_MU_vs_mu_NU.pdf")

#Diff simulations
Curve=readdlm("examples/Lih_data_Nmu_1_simul.output");
Curve=reshape(Curve,(100,25,:))
plot(Curve[1,:,2],Curve[1,:,1],linewidth=1,marker =:circle,markersize=2.5,label="tol_MPS:$(Curve[1,1,4])",title="2-sites DMRG with different initial guesses LiH: E_mu vs N_mu ",legend=false)
for i=2:80
    plot!(Curve[i,:,2],Curve[i,:,1],linewidth=1,marker =:circle,markersize=2.5,label="tol_MPS:$(Curve[i,1,4])")
end
xlabel!("N_\\mu")
ylabel!("E_\\mu")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/E_MU_vs_mu_NU_simul.pdf")


#Diff orders
Curve=readdlm("examples/Lih_data_Nmu_1_orders.output");
Curve=reshape(Curve,(25,37,:))
plot(Curve[:,1,2],Curve[:,1,1],linewidth=1,marker =:circle,markersize=2.5,label="tol_MPS:$(Curve[1,1,4])",title="2-sites DMRG with different orders LiH: E_mu vs N_mu ",legend=false)
for i=2:37
    plot!(Curve[:,i,2],Curve[:,i,1],linewidth=1,marker =:circle,markersize=2.5)
end
xlabel!("N_\\mu")
ylabel!("E_\\mu")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/E_MU_vs_mu_NU_orders.pdf")

#search of local minima in a specific point
true_energy=[ -7.876634127327094]
energies=readdlm("energies_order.output");
plot(abs.(S .-(-7.876634127327094))./abs(-7.876634127327094),yaxis=:log,linewidth=1,marker=:circle,label="Error vs number of simulations",legend=:outertopright)
plot!(abs.(energies .-(-7.876634127327094))./abs(-7.876634127327094),yaxis=:log,linewidth=1,marker=:circle,label="Error vs number of permutations",legend=:outertopright)
plot(S,linewidth=1,marker=:circle,label="Error vs number of simulations",legend=:outertopright)
plot!(energies,linewidth=1,marker=:circle,label="Error vs number of permutations",legend=:outertopright)
plot!(ones(length(S)).*true_energy,linewidth=2,label="Energy",legend=:outertopright)
xlabel!("Number of simulations/permuations")
ylabel!("E")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/local_minimas_curve.pdf")
plot(N_,linewidth=1,marker=:circle,label="Error vs number of simulations",legend=:outertopright)
plot!(electrons,linewidth=1,marker=:circle,label="Error vs number of permutations")
xlabel!("Number of simulations/permuations")
ylabel!("Number of electrons")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/local_minimas_curve_elec.pdf")

#Continuity test
Test_c=readdlm("examples/continuity_lih.output");
λ=[1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 0.5 1 5 10]
plot(λ',abs.(Test_c[2:end,1] .- (-7.88240342))./abs(-7.88240342),xaxis=:log,linewidth=1,marker=:circle)
plot(λ',Test_c[2:end,1],linewidth=1,marker=:circle)
λ=[0 1e-10 1e-9  1e-5 1e-2 1e-1 0.5 1 3 10 100]
E=[ -20.147797077859526  -20.14779707603147  -20.147797059405452  -20.147578203160258  -19.93012503196069 -17.971015381779555 -10.876206802081846  -7.882403424260322 -4.265094496346308  -3.744428615844959   -3.7444286154324806]
rks=[1 1 1 1 1 1 12 17 6 2 2]
plot(λ'[2:end],E'[2:end],xaxis=:log,marker =:circle,markersize=2.5,label="LIH molecule",legend=:topleft)
xlabel!("\\lambda")
ylabel!("E")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/Continuity_1.pdf")

plot(λ'[2:end],rks'[2:end],xaxis=:log,yaxis=:log,marker =:circle,markersize=2.5,label="LIH molecule",legend=:topleft)
xlabel!("\\lambda")
ylabel!("Maximum ranks")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/Continuity_2.pdf")

Test_c=readdlm("Continuity2.output");
l = @layout [a ; b c]
p1 = plot(Test_c[2:end,5],Test_c[2:end,1],xaxis=:log,linewidth=1,legend=false,marker =:circle,markersize=2.5,xlabel="\\lambda",ylabel="E")
p2 = plot(Test_c[2:end,5],Test_c[2:end,2],xaxis=:log,linewidth=1,legend=false,marker =:circle,markersize=2.5,xlabel="\\lambda",ylabel= "Maximum ranks")
p3 = plot(Test_c[2:end,5],Test_c[2:end,4],xaxis=:log,linewidth=1,legend=false,marker =:circle,markersize=2.5,xlabel="\\lambda",ylabel="\\mu")
plot(p1, p2, p3, layout = l)
plot(Test_c[2:end,6],[Test_c[2:end,1] Test_c[2:end,3] Test_c[2:end,5]],xaxis=:log,linewidth=1,legend=false,marker =:circle,markersize=2.5,xlabel="\\lambda",ylabel=["E" "Maximum ranks" "Chemical potential \\mu"],layout = 3)
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/Modified_hamiltonian.pdf")

#Disso with parametrized hamiltonian
Test_c=readdlm("Continuity2_disso.output");
E=[]
for i=1:size(Test_c,1)
    if isapprox(Test_c[i,5],1.0,atol=10e-6)
        push!(E,Test_c[i,1]);
    end
end
Dissos=readdlm("Lih_data.output");
Dissos=reshape(Dissos,9,15,6);
Bond_lengths_lih=[0.506  0.756 1.011 1.261 1.511 1.761 2.011 2.261 2.761 3.261 3.761 4.261 5.261 7.261 10.261]
Energy3=[-7.069466454629781 -7.58310908456387 -7.789705342926 -7.863635333395484 -7.8825286901531495 -7.876634127326262 -7.860265796202383 -7.840973635325725 -7.808618390149074 -7.7919561911536075 -7.7856146936559405 -7.783432184429271 -7.782501479366335 -7.7824189319850365 -7.782418427943522]
plot(Bond_lengths_lih',Energy3',linewidth=1,marker=:circle,label="Full-CI",title="Dissociation curve for LiH with 2-sites DMRG")
plot!(Bond_lengths_lih',Dissos[1,:,1],linewidth=1,marker=:circle,label="Without correction",title="Dissociation curve for LiH with 2-sites DMRG")
plot!(Bond_lengths_lih',E,linewidth=1,marker=:circle,label="With correction",title="Dissociation curve for LiH with 2-sites DMRG")
xlabel!("Bond length(Å)")
ylabel!("E")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/Correction.pdf")
Dissos_FCI=readdlm("Dissos_FCI2.output");
Bond_lengths_lih=[0.506  0.756 1.011 1.261 1.511 1.761 2.011 2.261 2.761 3.261 3.761 4.261 5.261 7.261 10.261]
Energy3=[-7.069466454629781 -7.58310908456387 -7.789705342926 -7.863635333395484 -7.8825286901531495 -7.876634127326262 -7.860265796202383 -7.840973635325725 -7.808618390149074 -7.7919561911536075 -7.7856146936559405 -7.783432184429271 -7.782501479366335 -7.7824189319850365 -7.782418427943522]
plot(Bond_lengths_lih'[1:end-1],Dissos_FCI[1:end-1,1],linewidth=1,marker=:circle,label="Exact solution for tol=1e-2",title="Dissociation curve for LiH with 2-sites DMRG")
plot!(Bond_lengths_lih',Dissos[1,:,1],linewidth=1,marker=:circle,label="Without correction tol=1e-2",title="Dissociation curve for LiH with 2-sites DMRG")
plot!(Bond_lengths_lih',Energy3',linewidth=1,marker=:circle,label="Full-CI",title="Dissociation curve for LiH with 2-sites DMRG")
xlabel!("Bond length(Å)")
ylabel!("E")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/dissos.pdf")

##

guess1=readdlm("intial_guess1.output");
guess2=readdlm("intial_guess2.output");
guess3=readdlm("intial_guess3.output");
guess4=readdlm("guess_1_site.output");
plot(guess2[:,1],marker =:circle,markersize=2.5,title="DMRG ground state energy for random intial guesses",label="Tol_MPS=1e-2")

plot(guess4[:,1],marker =:circle,markersize=2.5,title="DMRG ground state energy for random intial guesses",label="1-site DMRG Tol_MPS=1e-3")
plot!(guess3[:,1],marker =:circle,markersize=2.5,title="DMRG ground state energy for random intial guesses",label="2-site DMRG Tol_MPS=1e-3")
plot!(ones(length(guess3[:,1])).*-7.876439010719829,markersize=2.5,title="DMRG ground state energy for random intial guesses",label="Exact solution Tol_MPS=1e-3")
ylabel!("E")
xlabel!("Number of runs")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/random_simulations_sites.pdf")

guess5=readdlm("examples/Lih_orders_low2.output");
plot(guess5[:,1],marker =:circle,markersize=2.5,title="DMRG ground state energy for random intial guesses",label="Energy vs number of permutations")
plot!(guess2[:,1],marker =:circle,markersize=2.5,title="DMRG ground state energy for random intial guesses",label="Energy vs number of runs")
ylabel!("E")
xlabel!("Number of runs/permutations")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/orders.pdf")

##DW
tol_MPS_all=[5e-2 4e-2 3e-2 2e-2 1e-2 5e-3 4e-3 3e-3 2e-3 1e-3 5e-4 2e-4 1e-4 8e-5 5e-5 1e-6]
DW_lih=readdlm("Discarded_lih.output");
plot!(DW_lih[1:end,4],abs.(DW_lih[1:end,1].- (-7.88240342)),yaxis=:log,xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="LiH molecule 2-sites DMRG",legend=:topleft)
E_DW_lih= readdlm("Lih_discarded_weight_1_site.output");
plot!(E_DW_lih[1:end,3],abs.(E_DW_lih[1:end,1].-(-7.88240342)),yaxis=:log,xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="LiH molecule 1-site DMRG")
xlabel!("Discarded weight")
ylabel!("\\Delta E")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/disc_w_lih.pdf")

DW_Be2=readdlm("examples/Discarded_Be2.output");
plot!(DW_Be2[:,4],abs.(DW_Be2[:,1].- (-28.80434531)),xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="Be2 molecule")
DW_C2=readdlm("examples/Discarded_C2.output");
plot!(DW_C2[:,4],abs.(DW_C2[:,1].- (-75.012578266523052)),xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="C2 molecule")
DW_N2=readdlm("examples/Discarded_N2.output");
plot!(DW_N2[:,4],abs.(DW_N2[:,1].- (-107.65277153)),yaxis=:log,xaxis=:log,marker =:circle,markersize=2.5,title="DMRG ground state energy vs Discarded weight",label="N2 molecule")
##Disso for mPO tols
Dissos=readdlm("examples/Full_CI_MPO.output");
Dissos=reshape(Dissos,7,15,5);
Bond_lengths_lih=[0.506  0.756 1.011 1.261 1.511 1.761 2.011 2.261 2.761 3.261 3.761 4.261 5.261 7.261 10.261]
Energy3=[-7.069466454629781 -7.58310908456387 -7.789705342926 -7.863635333395484 -7.8825286901531495 -7.876634127326262 -7.860265796202383 -7.840973635325725 -7.808618390149074 -7.7919561911536075 -7.7856146936559405 -7.783432184429271 -7.782501479366335 -7.7824189319850365 -7.782418427943522]
plot(Bond_lengths_lih',Energy3',linewidth=1,marker=:circle,label="Full-CI",title="Dissociation curve for LiH with 2-sites DMRG with fixed MPS accuracy")
for i=1:7
    plot!(Bond_lengths_lih', Dissos[i,:,1],linewidth=1,marker=:circle,label="tol_MPO=$(Dissos[i,1,5])",title="Dissociation curve for LiH with 2-sites DMRG")
end
xlabel!("Bond length(Å)")
ylabel!("E")
Plots.savefig("DMRG_preliminary_test/figures/pres/sites/disso_MPO.pdf")
