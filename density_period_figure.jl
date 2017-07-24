# shows correlation between planet density vs. flux with an emphasis on TTV:
using CGS
using PyPlot

data = readdlm("exoplanet_eu.dat",',')
#  1          2      3     4     5      6      7      8           9        10      11 12      13   14    15      16    17    18    19      20 21 
# Name       Mass   sM1   sM2   Rad    sR1    sR2     P          sP1      sP2    Meth Ms     sMs1 sMs2   Rs     sRs1  sRs2  Teff   sT1    sT2 mult
mass   = vec(data[:,2])
sM1 = vec(data[:,3])
sM2 = vec(data[:,4])
radius = vec(data[:,5])
sR1 = vec(data[:,6])
sR2 = vec(data[:,7])
period = vec(data[:,8])
Teff = vec(data[:,18])
method = vec(data[:,11])
mult = vec(data[:,21])

nplanet = length(mass)
#density = 3./4./pi*mass.*MSUN./(radius.*RSUN).^3
density = zeros(nplanet)
sdensity = zeros(nplanet)
# this should have mstar, not mass (of planet):
mstar = convert(Vector{Float64},data[:,12])
flux = (Teff./5777.).^4./(period./365.25).^(4//3).*mstar.^(2//3)

irv =[]
irv1 =[]
ittv = []
fig,ax = subplots()
nrand = 1000
for i=1:nplanet
  if mass[i] > (abs(sM1[i])*3.0)
# Monte Carlo density (not quite correct...  M/R are correlated in TTV case)
  r1 = randn(nrand)
  mr = zeros(nrand)
  rr = zeros(nrand)
  for j=1:nrand
    if r1[j] < 0
     mr[j]=mass[i]-sM1[i]*r1[j]
    else 
     mr[j]=mass[i]+sM2[i]*r1[j]
    end
  end
  r2 = randn(nrand)
  for j=1:nrand
    if r2[j] < 0
     rr[j]=radius[i]-sR1[i]*r2[j]
    else 
     rr[j]=radius[i]+sR2[i]*r2[j]
    end
  end
  density[i] = mean(3./4./pi*mr.*MSUN./(rr.*RSUN).^3)
  sdensity[i] = std(3./4./pi*mr.*MSUN./(rr.*RSUN).^3)
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
      if period[i] > 30
        println(data[i,:])
      end
  end
  end
end



#ax[:loglog](flux[irv],density[irv],"ro",alpha=0.5,label = "RV",s=(Teff[irv]./5777.).^2.*10.)
#ax[:scatter](flux[irv],density[irv],c="r",alpha=0.25,label = "RV",s=(Teff[irv]./5777.).^2*300.)
#ax[:scatter](flux[irv],density[irv],c="r",alpha=0.25,label = "RV",s=mass[irv].*50.)
#ax[:scatter](flux[irv],density[irv],c="r",alpha=0.25,label = "RV",s=mult[irv]*50.)
#ax[:scatter](period[irv],density[irv],c="r",alpha=0.25,label = "RV",s=mult[irv]*50.)
ax[:errorbar](period[irv],density[irv],sdensity[irv],c="r",alpha=0.75,fmt="o",label="RV")
#ax[:scatter](period[irv1],density[irv1],c="g",alpha=0.25,label = "RV single",s=mult[irv1]*50.)
ax[:errorbar](period[irv1],density[irv1],sdensity[irv1],c="g",alpha=0.75,label = "RV single",fmt="o")
ax[:set_xscale]("log")
ax[:set_yscale]("log")
#ax[:loglog](flux[ittv],density[ittv],"bo",alpha=0.5,label="TTV",s=(Teff[ittv]./5777.).^2)
#ax[:scatter](flux[ittv],density[ittv],c="b",alpha=0.25,label="TTV",s=(Teff[ittv]./5777.).^2.*300.)
#ax[:scatter](flux[ittv],density[ittv],c="b",alpha=0.25,label="TTV",s=mass[ittv].*50.)
#ax[:scatter](flux[ittv],density[ittv],c="b",alpha=0.25,label="TTV",s=mult[ittv].*50.)
#ax[:scatter](period[ittv],density[ittv],c="b",alpha=0.5,label="TTV",s=mult[ittv].*50.)
ax[:errorbar](period[ittv],density[ittv],sdensity[ittv],c="b",alpha=0.75,label="TTV",fmt="o")
#ax[:set_title]("Planet density versus period")
ax[:axis]([1,2e2,5e-3,10])
#ax[:set_xlabel]("Flux/(Solar Constant)")
ax[:set_xlabel]("Period [d]")
ax[:set_ylabel]("Density [g/cc]")
ax[:legend](loc="lower left")
