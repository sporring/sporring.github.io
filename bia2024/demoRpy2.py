# This script gives an example on how rpy2 can be used to work with R and the spatstat and
# friends package from python. The following has been tested on an intel and m1 mac. You may
# need to customize for your setup
# 
# Jon Sporring
# Department of Computer Science
# University of Copenhagen
# August 11, 2024
#
# 1. Install R, which to my experience works best directly from  https://cran.r-project.org/
# then start R and install some packes:
#   install.packages("spatstat")
#   install.packages("lazyeval")
#   install.packages("GET")
# In R you can get access to lots of manuals, e.g.,
#   library("spatstat")
#   ?ppp
# and there are also more general descriptions
#   library("GET")
#   vignette("pointpatterns")
#
# 2. Install rpy2, of which to my experience only the pip version is stable. I.e., create
# an environment and pip install: 
#   conda create -n rpy2 python=3.12 matplotlib
#   conda activate rpy2
#   pip3 install rpy2==3.5.12 # pcre2-8 library problems on some systems for later versions
# now you should be able to check the basic settings of rpy2
#   python -m rpy2.situation 
#
# 3. It is not necessary, and frankly, I find it not very stable, but if you wish to produce
# plots with R from Python on a Mac, then on a mac you also need the X11 graphics driver,
# which is called XQuartz, https://www.xquartz.org/
#
# 4. Now you can start Python, and I find the introduction on
#    https://rpy2.github.io/doc/latest/html/introduction.html
# instructive

# The following contains 2 demonstrations of working with spatstat, GET, numpy, and matplotlib.
# 
# first we import a bunch of basic stuff and prepare working with rpy2
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
import numpy as np
from rpy2.robjects import FloatVector
base = importr('base')
spatstat = importr('spatstat')
# shortcut
ro = robjects.r

# The following is an example of how to create a random 2-dimensional point pattern and make a
# python wrapper for the Ripley's K function, Kest:
# 
# First we generate a random set of coordinates
x = FloatVector(np.random.rand(100,1))
y = FloatVector(np.random.rand(100,1))
# then we create a point pattern data structure
bounds = FloatVector([0,1])
pp = ro.ppp(x, y, window = ro.owin(bounds,bounds))
# Finally, we create the Ripley's K-function data structure, extract its function part and link it
# to Python as the function f
f = base.as_function(ro.Kest(pp))
# This function can now be called, e.g., to be plottet
n=100
t = [i/float(n) for i in range(n+1)]
ft = [f(i)[0] for i in t] # the return value of f(i) is a vector, so we need to extract its element.
plt.plot(t,ft)
plt.show()

# The next example demonstrates how we can setup the envelope function, which relies on an R
# expression, i.e., a function to be called in the simulation loop of the envelope function. This
# is similar to a Python lambda function. The script also demonstrates how to extract 2 types
# of attributes from an R object.
#
# The original R code is as follows:
#   library("GET")
#   library("spatstat.model")
#   library("ggplot2")
#   X <- spruces
#   X
#   env <- envelope(X, nsim=1999, savefuns=TRUE, simulate=expression(runifpoint(ex=X)), verbose=FALSE)
#   res <- global_envelope_test(env)
#   plot(res)
# My Rpy2 version is
GET = importr("GET") # We don't access GET but it is still needed to load the library
model = importr("spatstat.model")
X = ro.spruces
print(X)
# The spatstat envolope function requires an unevaluated expression in order to perform repeated
# experiments. To do this from Python, we must move the relevant data to the R environment and use
# the rpy2.robjects.language.LangVector.from_string functionality 
robjects.globalenv["X"] = X
sim = rpy2.robjects.language.LangVector.from_string('expression(runifpoint(ex=X))')
env = ro.envelope(X, fun=ro.Lest, nsim=1999, savefuns=True, simulate=sim, verbose=False)
res = ro.global_envelope_test(env)
print(res) # Two type of data: attributes and elements in a double list
print(base.attributes(res)) # the full list
# One can plot it directly with
#   ro.plot(res)
# However, I find the Xquartz/X11 integration poor on my mac, so I prefer to extract data and plot
# in python
p = ro.attr(res,'p')[0] # some needs to be indexed liked this
lo = res.rx2('lo')
hi = res.rx2('hi')
r = res.rx2('r')
obs = res.rx2('obs')
central = res.rx2('central')
plt.plot(r,lo, 'k--')
plt.plot(r,hi, 'k--')
plt.plot(r,central, 'k-')
plt.plot(r,obs, 'k-')
plt.show()
