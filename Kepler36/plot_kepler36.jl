using PyPlot

data = readdlm("kepler36_lc.dat")

time  = vec(data[:,1])
fflat = vec(data[:,2])
ntphot = 65261
for i=1:ntphot
  if fflat[i] > 1.0005
     fflat[i] = 1.0
  end
end 

# Make a river plot:

#t0 = [54960.9753,54955.9132]
t0 = [54960.8553,54955.9132]
#p  = [13.83989,16.23855]
p  = [13.84889,16.23255]

i0 = 1

i11 = round(Int64,ceil((time[1]-t0[1])/p[1]))
i21 = round(Int64,floor((time[ntphot]-t0[1])/p[1]))

n1 = round(Int64,i21-i11+1)
println("Number of transits of inner planet: ",n1)

i12 = round(Int64,ceil((time[1]-t0[2])/p[2]))
i22 = round(Int64,floor((time[ntphot]-t0[2])/p[2]))
n2 = round(Int64,i22-i12+1)
println("Number of transits of outer planet: ",n2)

scatter(time,fflat)

for i=i11:i21
#  plot((t0[1]+p[1]*i)*[1,1],[0.9995,1.0001])
end
for i=i12:i22
#  plot((t0[2]+p[2]*i)*[1,1],[0.9995,1.0001])
end


nwidth = 50
twidth = 25.0  # Width of river plot in hours
im1 = zeros(n1,nwidth)
im2 = zeros(n2,nwidth)

# Loop over time array until we come close to 
i1 = 1
for i=i11:i21
  tc = t0[1]+p[1]*i
  for j=1:nwidth
# Compute the central time of the pixel:
    tpix = tc+(j-0.5-nwidth/2.)*twidth/nwidth/24.
# Loop over i1 until we pass it:
    while time[i1] < tpix && i1 < ntphot
      i1+=1
    end
# Add to the pixel:
    im1[i-i11+1,j]=(fflat[i1-1]*(time[i1]-tpix)/(time[i1]-time[i1-1])
       +fflat[i1]*(tpix-time[i1-1])/(time[i1]-time[i1-1]))
    if im1[i-i11+1,j] < 0.9997
      im1[i-i11+1,j] = 1.0
    end
  end
  im1[i-i11+1,:]/=median(im1[i-i11+1,:])
# Move on to the next transit
end

# Now do the same for the other planet:
# Loop over time array until we come close to
i1 = 1
dt = (collect(linspace(0,nwidth-1,nwidth)+0.5)/nwidth-0.5)*twidth
for i=i12:i22
  tc = t0[2]+p[2]*i
  for j=1:nwidth
# Compute the central time of the pixel:
    tpix = tc+(j-0.5-nwidth/2.)*twidth/nwidth/24.
# Loop over i1 until we pass it:
    while time[i1] < tpix && i1 < ntphot
      i1+=1
    end
# Add to the pixel:
#    println(time[i1-1]," ",tpix," ",time[i1])
    im2[i-i12+1,j]=(fflat[i1-1]*(time[i1]-tpix)+fflat[i1]*(tpix-time[i1-1]))/(time[i1]-time[i1-1])
  end
  im2[i-i12+1,:]/=median(im2[i-i12+1,:])
#  clf()
#  plot(dt,vec(im2[i-i12+1,:]))
#  plot((time-tc)*24.,fflat)
#  axis([-twidth/2,twidth/2,0.995,1.005])
#  read(STDIN,Char)
# Move on to the next transit
end


fig,axes = subplots(1,2)

ax = axes[1]

ax[:imshow](im1,interpolation="nearest")

ax = axes[2]
#fig,axes = subplots()

ax[:imshow](im2,interpolation="nearest")
#axes[:imshow](im2,interpolation="nearest")
