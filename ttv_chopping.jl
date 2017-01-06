# Create an example of two systems with same TTV amplitude,
# but different chopping amplitudes/eccentricities.
using PyPlot
include("regress.jl")

path_to_ttvfast = "/Users/ericagol/Software/TTVFast/jl_version/"

# Demo how to call TTV fast from Julia
if !isdefined(:TTVFAST_INITIALIZED)
  include(string(path_to_ttvfast,"TTVFast.jl"))
  using TTVFast
end

# Integration parameters
inflag = 0
t_start=0.0
#t_stop=4*365.2425
t_stop=1000.0
#t_stop=180.0
dt = 0.5

# initializing array of input parameters to TTVFast
# should make into nice function
#p2 = 1.48
p1 = 10.0
p2 = 1.52*p1
m10 = 1e-6
m20 = 1e-6
nplanets=2;  # hardwired for now
  duration = t_stop-t_start
  p = Array(Cdouble,2+nplanets*7);
#  p[1] = 4pi/365.2425^2; # G
  p[1] = 0.000295994511; # G
  p[2] = 1.0; # Mstar
  num_events=0;
  i=1
  p[2+7*(i-1)+1] = m10 ;   # Planet Mass
  p[2+7*(i-1)+2] = p1 ; # Period
  p[2+7*(i-1)+3] = 0.00 ;  # Eccentricity
  p[2+7*(i-1)+4] = 90.0;    # Inclination
  p[2+7*(i-1)+5] = 0.0 # 2pi*rand(); # Longitude of Ascending Node
  p[2+7*(i-1)+6] = 0.0; # Argument of Pericenter
  p[2+7*(i-1)+7] = 360*rand(); # Mean Anomaly
  num_events += ceil(Integer, duration/p[2+7*(i-1)+2] +1); # /* large enough to fit all the transits calculated by the code*/
  i=2
  p[2+7*(i-1)+1] = m20 ;   # Planet Mass
  p[2+7*(i-1)+2] = p2 ; # Period
  p[2+7*(i-1)+3] = 0.00 ;  # Eccentricity
  p[2+7*(i-1)+4] = 90.0;    # Inclination
  p[2+7*(i-1)+5] = 0.0 # 2pi*rand(); # Longitude of Ascending Node
  p[2+7*(i-1)+6] = 0.0; # Argument of Pericenter
  p[2+7*(i-1)+7] = 360*rand(); # Mean Anomaly
  num_events += ceil(Integer, duration/p[2+7*(i-1)+2] +1); # /* large enough to fit all the transits calculated by the code*/

ttvfast_input = ttvfast_inputs_type(p, t_start=t_start, t_stop=t_stop, dt=dt)

incl_rvs = false # Make sure first rv_time is _after_ t_start+dt/2 or else won't get any RV outputs
if incl_rvs
  num_rvs = 100
  rv_times = collect(linspace(t_start+0.501*dt,t_stop, num_rvs))
  ttvfast_output = ttvfast_outputs_type(num_events ,rv_times)
else
  ttvfast_output = ttvfast_outputs_type(num_events)
end

println(STDERR, "# About to call TTVFast")
tic();
ttvfast!(ttvfast_input,ttvfast_output)
toc();

ntrans_max = ttvfast_output.max_num_events
ntrans_max -= 4
ntrans = [0,0]
i1=Int64[]
i2=Int64[]
for i=1:ntrans_max
  tmp = get_event(ttvfast_output,i)
  if tmp.planet  == 0 && tmp.time != 0.0
     ntrans[1] +=1
     i1 = [i1;i]
  end
  if tmp.planet == 1 && tmp.time != 0.0
     ntrans[2] +=1
     i2 = [i2;i]
  end
end

tt1_0 = zeros(ntrans[1])
epoch1 = zeros(ntrans[1])
for i=1:ntrans[1]
  tmp = get_event(ttvfast_output,i1[i])
  tt1_0[i]=tmp.time
  epoch1[i]=tmp.epoch
end
tt2_0 = zeros(ntrans[2])
epoch2 = zeros(ntrans[2])
for i=1:ntrans[2]
  tmp = get_event(ttvfast_output,i2[i])
  tt2_0[i]=tmp.time
  epoch2[i]=tmp.epoch
end

# Now, carry out a regression to plot TTVs:

y = tt1_0
x = zeros(2,ntrans[1])
x[1,:]=1.0
x[2,:]=epoch1
sigy = ones(ntrans[1])
coeff1,cov =regress(x,y,sigy)

y = tt2_0
x = zeros(2,ntrans[2])
x[1,:]=1.0
x[2,:]=epoch2
sigy = ones(ntrans[2])
coeff2,cov =regress(x,y,sigy)

#fig,axes = subplots(2,1)
#ax=axes[1,1]
tt10_0 =coeff1[1]+coeff1[2]*epoch1
ttv1_0  = (tt1_0-tt10_0)*24.*60.
#ax[:plot](tt10,ttv1,"ro")
#ax[:axis]([0,80,-10,10])
#ax[:set_ylabel]("TTV [min]")
#ax=axes[2,1]
tt20_0 = coeff2[1]+coeff2[2]*epoch2
ttv2_0  = (tt2_0-tt20_0)*24.*60.
#ax[:plot](tt20,ttv2,"bo")
#ax[:axis]([0,80,-10,10])
#ax[:set_xlabel]("Time [day]")
#ax[:set_ylabel]("TTV [min]")

# Next, write code to re-run with a different chopping amplitude:

# Let's use Lithwick et al.'s notation:

# We want m_2_0 * (-f)  = m_2 * (-f - 3/2 * (f* e_1*cos(om_1) + g * e_2 * cos(om_2))/Delta)
# m_1_0 * (-g)  = m_2 * (-g + 3/2 * (f* e_1*cos(om_1) + g * e_2 * cos(om_2))/Delta)

j0 = 3
Delta =  p2/p1*(j0-1)/j0-1.0
# for j0=3, Table 3 in Lithwick et al.:
f = -2.025 + 6.21*Delta
g = 2.484-5.99*Delta
e1 = 0.04
e2 = 0.04
om1 = 0.0
om2 = pi
zreal = 1.5*(f*e1*cos(om1)+g*e2*cos(om2))/Delta
zimag = 1.5*(f*e1*sin(om1)+g*e2*sin(om2))/Delta
m2 = f/(f+zreal)*m20
m1 = g/(g-zreal)*m10
println("M1: ",m1,"; M2: ",m2)

p0 = zeros(length(p))
for i=1:length(p)
  p0[i]=p[i]
end
i=1
p[2+7*(i-1)+1] = m1 ;   # Planet Mass
p[2+7*(i-1)+3] = e1 ;  # Eccentricity
p[2+7*(i-1)+6] = om1*180/pi + 90;  # Argument of peri
p[2+7*(i-1)+7] -=  p[2+7*(i-1)+6] ;  # Mean longitude - offset by longitude of periastron
i=2
p[2+7*(i-1)+1] = m2 ;   # Planet Mass
p[2+7*(i-1)+3] = e2 ;  # Eccentricity
p[2+7*(i-1)+6] = om2*180/pi + 90;  # Argument of peri
p[2+7*(i-1)+7] -= p[2+7*(i-1)+6] ;  # Mean longitude - offset by longitude of periastron

ttvfast_input = ttvfast_inputs_type(p, t_start=t_start, t_stop=t_stop, dt=dt)

incl_rvs = false # Make sure first rv_time is _after_ t_start+dt/2 or else won't get any RV outputs
if incl_rvs
  num_rvs = 100
  rv_times = collect(linspace(t_start+0.501*dt,t_stop, num_rvs))
  ttvfast_output = ttvfast_outputs_type(num_events ,rv_times)
else
  ttvfast_output = ttvfast_outputs_type(num_events)
end

println(STDERR, "# About to call TTVFast")
tic();
ttvfast!(ttvfast_input,ttvfast_output)
toc();

ntrans_max = ttvfast_output.max_num_events
ntrans = [0,0]
i1=Int64[]
i2=Int64[]
for i=1:ntrans_max
  tmp = get_event(ttvfast_output,i)
#    println(tmp.planet," ",tmp.time," ",tmp.rsky," ",tmp.vsky)
end
for i=1:ntrans_max
  tmp = get_event(ttvfast_output,i)
  if tmp.planet  == 0 && tmp.time != 0.0
     ntrans[1] +=1
     i1 = [i1;i]
  end
  if tmp.planet == 1 && tmp.time != 0.0
     ntrans[2] +=1
     i2 = [i2;i]
  end
end

tt1 = zeros(ntrans[1])
epoch1 = zeros(ntrans[1])
for i=1:ntrans[1]
  tmp = get_event(ttvfast_output,i1[i])
  tt1[i]=tmp.time
  epoch1[i]=tmp.epoch
end
tt2 = zeros(ntrans[2])
epoch2 = zeros(ntrans[2])
for i=1:ntrans[2]
  tmp = get_event(ttvfast_output,i2[i])
  tt2[i]=tmp.time
  epoch2[i]=tmp.epoch
end



# Now, carry out a regression to plot TTVs:

y = tt1
x = zeros(2,ntrans[1])
x[1,:]=1.0
x[2,:]=epoch1
sigy = ones(ntrans[1])
coeff1,cov =regress(x,y,sigy)

y = tt2
x = zeros(2,ntrans[2])
x[1,:]=1.0
x[2,:]=epoch2
sigy = ones(ntrans[2])
coeff2,cov =regress(x,y,sigy)

fig,axes = subplots(2,1)
ax=axes[1,1]
tt10 =coeff1[1]+coeff1[2]*epoch1
ttv1  = (tt1-tt10)*24.*60.
ax[:plot](tt10_0,ttv1_0,"go")
ax[:set_title](L"Planet 1, $P_1= 10$ days")
ax[:plot](tt10,ttv1,"ro")
ax[:plot](tt10_0,ttv1_0,"g-")
ax[:plot](tt10,ttv1,"r-")
ax[:axis]([0,1000,-0.5,0.5])
ax[:set_ylabel](L"O-C$_{TTV}$ [min]")
ax=axes[2,1]
tt20 = coeff2[1]+coeff2[2]*epoch2
ttv2  = (tt2-tt20)*24.*60.
ax[:plot](tt20_0,ttv2_0,"go")
ax[:set_title](L"Planet 2, $P_2=15.2$ days")
ax[:plot](tt20,ttv2,"ro")
ax[:plot](tt20_0,ttv2_0,"g-")
ax[:plot](tt20,ttv2,"r-")
ax[:axis]([0,1000,-0.5,0.5])
ax[:set_xlabel]("Time [day]")
ax[:set_ylabel](L"O-C$_{TTV}$[min]")
savefig("ttv_chopping.pdf", bbox_inches="tight")

