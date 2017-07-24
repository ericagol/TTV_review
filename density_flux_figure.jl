# shows correlation between planet density vs. flux with an emphasis on TTV:
using CGS
using PyPlot

data = readdlm("exoplanet_eu.dat",',')
#  1          2      3     4     5      6      7      8           9        10      11 12      13   14    15      16    17    18    19      20 21 
# Name       Mass   sM1   sM2   Rad    sR1    sR2     P          sP1      sP2    Meth Ms     sMs1 sMs2   Rs     sRs1  sRs2  Teff   sT1    sT2 mult
mass   = vec(data[:,2])
radius = vec(data[:,5])
period = vec(data[:,8])
Teff = vec(data[:,18])
method = vec(data[:,11])
mult = vec(data[:,21])

density = 3./4./pi*mass.*MSUN./(radius.*RSUN).^3
mstar = convert(Vector{Float64},data[:,12])
flux = (Teff./5777.).^4./(period./365.25).^(4//3).*mstar.^(2//3)

nplanet = length(mass)
irv =[]
irv1 =[]
ittv = []
for i=1:nplanet
  if mult[i] > 1
    if method[i] == "RV"
      irv =[irv;i]
      if period[i] > 30
        println(data[i,:])
      end
    else
      ittv = [ittv;i]
    end
  else
    irv1 = [irv1;i]
  end
end


fig,ax = subplots()

#ax[:loglog](flux[irv],density[irv],"ro",alpha=0.5,label = "RV",s=(Teff[irv]./5777.).^2.*10.)
#ax[:scatter](flux[irv],density[irv],c="r",alpha=0.25,label = "RV",s=(Teff[irv]./5777.).^2*300.)
#ax[:scatter](flux[irv],density[irv],c="r",alpha=0.25,label = "RV",s=mass[irv].*50.)
#ax[:scatter](flux[irv],density[irv],c="r",alpha=0.25,label = "RV",s=mult[irv]*50.)
ax[:scatter](period[irv],density[irv],c="r",alpha=0.25,label = "RV",s=mult[irv]*50.)
ax[:scatter](period[irv1],density[irv1],c="g",alpha=0.25,label = "RV single",s=mult[irv1]*50.)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
#ax[:loglog](flux[ittv],density[ittv],"bo",alpha=0.5,label="TTV",s=(Teff[ittv]./5777.).^2)
#ax[:scatter](flux[ittv],density[ittv],c="b",alpha=0.25,label="TTV",s=(Teff[ittv]./5777.).^2.*300.)
#ax[:scatter](flux[ittv],density[ittv],c="b",alpha=0.25,label="TTV",s=mass[ittv].*50.)
#ax[:scatter](flux[ittv],density[ittv],c="b",alpha=0.25,label="TTV",s=mult[ittv].*50.)
ax[:scatter](period[ittv],density[ittv],c="b",alpha=0.25,label="TTV",s=mult[ittv].*50.)
ax[:set_title]("Planet density versus period")
ax[:axis]([1,2e2,5e-3,10])
#ax[:set_xlabel]("Flux/(Solar Constant)")
ax[:set_xlabel]("Period [d]")
ax[:set_ylabel]("Density [g/cc]")
ax[:legend](loc="lower left")
