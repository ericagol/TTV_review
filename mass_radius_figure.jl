# shows correlation between planet density vs. flux with an emphasis on TTV:
using CGS
using PyPlot

data = readdlm("exoplanet_eu.dat",',')
#  1          2      3     4     5      6      7      8           9        10      11 12      13   14    15      16    17    18    19      20 21 
# Name       Mass   sM1   sM2   Rad    sR1    sR2     P          sP1      sP2    Meth Ms     sMs1 sMs2   Rs     sRs1  sRs2  Teff   sT1    sT2 mult
mass   = convert(Vector{Float64},data[:,2])
sM1 = convert(Vector{Float64},data[:,3])
sM2 = convert(Vector{Float64},data[:,4])
radius = convert(Vector{Float64},data[:,5])
sR1 = convert(Vector{Float64},data[:,6])
sR2 = convert(Vector{Float64},data[:,7])

period = convert(Vector{Float64},data[:,8])
Teff = convert(Vector{Float64},data[:,18])
method = vec(data[:,11])
mult = vec(data[:,21])

density = 3./4./pi*mass.*MEARTH./(radius.*REARTH).^3
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
    else
      ittv = [ittv;i]
    end
  else
    irv1 = [irv1;i]
  end
end


fig,ax = subplots()
#fig,ax = figure()

logf = log10.(flux)
fnorm = (logf-1)/3
colormapname = "RdYlBu_r"
ax[:errorbar](mass[irv],radius[irv],xerr=sM1[irv],yerr=sR1[irv],c="b",alpha=1.00,fmt=".")
sc = ax[:scatter](mass,radius,c=logf,alpha=1,cmap=colormapname)
fig[:colorbar](sc,label="log_10(Flux [F_earth])")
ax[:scatter](mass[irv],radius[irv],c=fnorm[irv],alpha=1.00,label = "RV",s=200,cmap=colormapname)
ax[:errorbar](mass[irv1],radius[irv1],xerr=sM1[irv1],yerr=sR1[irv1],c="b",alpha=1.00,fmt=".")
#sc = ax[:scatter](mass[irv1],radius[irv1],c=fnorm[irv1],alpha=1.00,marker="^",label = "RV single",s=200,cmap=colormapname) #,s=mult[irv1]*50.)
ax[:scatter](mass[irv1],radius[irv1],c=logf[irv1],alpha=1.00,marker="^",label = "RV single",s=200,cmap=colormapname) #,s=mult[irv1]*50.)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
ax[:errorbar](mass[ittv],radius[ittv],xerr=sM1[ittv],yerr=sR1[ittv],c="b",alpha=1.,fmt=".")
ax[:scatter](mass[ittv],radius[ittv],c=fnorm[ittv],alpha=1.0,marker = "s",label="TTV",s=200,cmap=colormapname)

ax[:set_title]("Planet radius versus mass ")
ax[:axis]([0.2,25.,0.5,12])
#ax[:set_xlabel]("Flux/(Solar Constant)")
ax[:set_xlabel]("Mass [M_earth]")
ax[:set_ylabel]("Radius [R_earth]")
ax[:legend](loc="upper left")
