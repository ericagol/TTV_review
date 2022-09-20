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

# Planet density:
density = 3./4./pi*mass.*MEARTH./(radius.*REARTH).^3
rstar = convert(Vector{Float64},data[:,15]) 
mstar = convert(Vector{Float64},data[:,12])
flux = (Teff./5777.).^4./(period./365.25).^(4//3).*mstar.^(2//3)

semi = (.25*(period*24.*3600.).^2.*GRAV.*mstar.*MSUN/pi^2).^(1./3.) # semi-major axis in cm
albedo = 0.0   # set planet albedo
# compute the "equilibrium" temperature of planet (assuming perfect recirculation):
tplanet = Teff.*(1.0-albedo)^.25.*sqrt(0.5*rstar.*RSUN./semi)
# compute sound speed:
cs = sqrt(5./3.*BOLTZMANN*tplanet./(MHYDROGEN*2.3))
# escape speed:
vesc = sqrt(GRAV*MEARTH*mass./(radius*REARTH))
# ratio
ratio = vesc./cs

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

logf = log10(flux)
fnorm = (logf-1)/3
colormapname = "RdYlBu_r"
#ax[:errorbar](ratio[irv],density[irv],xerr=sM1[irv],yerr=sR1[irv],c="b",alpha=1.00,fmt=".")
sc = ax[:scatter](ratio,density,c=logf,alpha=1,cmap=colormapname)
fig[:colorbar](sc,label="log_10(Flux [F_earth])")
ax[:scatter](ratio[irv],density[irv],c=fnorm[irv],alpha=1.00,label = "RV",s=200,cmap=colormapname)
#ax[:errorbar](ratio[irv1],density[irv1],xerr=sM1[irv1],yerr=sR1[irv1],c="b",alpha=1.00,fmt=".")
#sc = ax[:scatter](mass[irv1],radius[irv1],c=fnorm[irv1],alpha=1.00,marker="^",label = "RV single",s=200,cmap=colormapname) #,s=mult[irv1]*50.)
ax[:scatter](ratio[irv1],density[irv1],c=logf[irv1],alpha=1.00,marker="^",label = "RV single",s=200,cmap=colormapname) #,s=mult[irv1]*50.)
ax[:set_xscale]("log")
ax[:set_yscale]("log")
#ax[:errorbar](ratio[ittv],density[ittv],xerr=sM1[ittv],yerr=sR1[ittv],c="b",alpha=1.,fmt=".")
ax[:scatter](ratio[ittv],density[ittv],c=fnorm[ittv],alpha=1.0,marker = "s",label="TTV",s=200,cmap=colormapname)
ax[:plot]([6.,6.],[1e-3,1e3])
ax[:set_title]("Density vs. v_esc/c_s")
#ax[:axis]([0.5,25.,1,12])
#ax[:set_xlabel]("Flux/(Solar Constant)")
ax[:set_xlabel]("v_esc/c_s")
ax[:set_ylabel]("density [g/cc]")
ax[:legend](loc="upper left")
