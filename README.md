# TTV_review
Chapter for Springer's Exoplanet Handbook on TTV/TDV

exoplanet_eu.dat:  
Downloaded m < 25 M_earth planets with radius measurements from exoplanet.eu on 10/7/2016
Total of 157 planets.  Discarded planets for which lower limit is < 2-sigma from zero.
Left with 76 planets.  Had to look up parameters & uncertainties for some - used NeXScI
(http://exoplanetarchive.ipac.caltech.edu/index.html) & reference papers (found via
exoplanet.eu).  Main mass refs. are Marcy et al., Hadden & Lithwick, etc.

This matches the number in Jontof-Hutter's table as well, which are in the file planet_mr.txt

J-F misses GJ 1132, BD+20 594 b, -20d, -21b, -23b, -56b, -80bcde, -102d -113 -114
eu misses 289-b [ ], K2-38 (aka EPIC 204221263) [ ], Kepler-4b [ ], Kepler-10c [ ], -11b,c (<2-sig), -18b (<2-sig), -33def [ ]
  -89c [ ],  -98b [ ], -99b [ ], -106c [ ], K2-19c [ ], Corot-22b [ ]

So there should be a total of ~89 planets!

10/9/2016

Both the HARPS & Keck teams clamin that K2-3d is affected by stellar variability.  Density is
abnormally high, so we'll eliminate this one.

So, now I need to figure out how to make Dan's plot. [ ]

10/18/2016
To get the Julia call of TTVFast working (see https://github.com/kdeck/TTVFast/tree/master/jl_version), 
I had to make sure the path  is added for the TTVFast.jl function, and I also had to 
add an absolute path to the TTVFast.jl code since a relative path doesn't work when 
one moves to a different directory.   There probably is a better way to code this up
so that a path can be passed to the library in TTVFast.jl, but at the moment I'm not 
sure how to do this.

Okay, got the Julia version of TTVFast working, and managed to repeat runs with different
parameters.  Unfortunately TTVFast seems to be missing transits on occasion...  not sure
why. [ ]  I should probably write parameters to a file, call TTVFast directly, and see if the
problem persists.

10/19/2016

Okay, I successfully made a chopping comparison plot with masses differing by a factor of ~10:

Here are the zero-eccentricity & e_1=e_2 = 0.04 eccentricity parameters:
julia> p0
16-element Array{Float64,1}:
   0.000295995
   1.0        
   1.0e-6     
  10.0        
   0.0        
  90.0        
   0.0        
   0.0        
 293.193      
   1.0e-6     
  15.2        
   0.0        
  90.0        
   0.0        
   0.0        
 183.256      

julia> p
16-element Array{Float64,1}:
   0.000295995
   1.0        
   1.09465e-7 
  10.0        
   0.04       
  90.0        
   0.0        
  90.0        
 203.193      
   9.03319e-8 
  15.2        
   0.04       
  90.0        
   0.0        
 270.0        
 -86.7444     

Skewness?  Nesvorny - perturbing planet may have to be inclined - would give TDVs.

There is a scatter in the densities of super-Earths.  5 M_earth bodies aren't necessarily
rocky.

Maybe cutoff at 15 M_earth.

Add composition lines to mass-radius plot. [ ]  Puffy planets have low insolation.

What is high-density massive planet? M ~ 15 M_earth; R < 2 R_earth.  Make sure there
aren't outliers with poor data.  [ ]  Make error bars black. [ ]  Check on lowest mass
planet - why are there two tickmarks at end of error bars? [ ]  Plot dots on top of
error bars. [ ]  Don't need a title - describe in caption. [ ]  Add M_earth symbol. [ ]
Add in more tickmarks. [ ]  Make axes & tickmarks heavier. [ ]
