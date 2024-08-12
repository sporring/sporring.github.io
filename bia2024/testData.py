# A solution to the exercise: 
#    Use GET's global_envelope_test to test whether x and/or z are likely to be random
# posed at the DTU, UCPH, AAU summerschool 2024
#
# Jon Sporring
# Department of Computer Science
# University of Copenhagen
# August 11, 2024

# load the two datasets
import numpy as np
x = np.loadtxt("x.csv",delimiter=",", dtype=float)
z = np.loadtxt("z.csv",delimiter=",", dtype=float)

import matplotlib.pyplot as plt
ax1 = plt.subplot(3,2,1)
ax1.plot(x[:,0],x[:,1],'k.')
ax1.set_title('x')
ax2 = plt.subplot(3,2,2)
ax2.plot(z[:,0],z[:,1],'k.')
ax2.set_title('z')

# setup rpy2 to interface with the spatstat and GET packages from R
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
import numpy as np
from rpy2.robjects import FloatVector
base = importr('base')
spatstat = importr('spatstat')
ro = robjects.r

# Estimate L from samples
def sampleL(x, n=100, bounds = FloatVector([0,1])):
    t = [i/float(n) for i in range(n+1)]
    ppx = ro.ppp(FloatVector(x[:,0]), FloatVector(x[:,1]), window = ro.owin(bounds,bounds))
    fx = base.as_function(ro.Lest(ppx))
    fxt = [fx(i)[0] for i in t] 
    return (t,fxt)

(t,fxt) = sampleL(x)
(_,fzt) = sampleL(z)
ax3 = plt.subplot(3,2,3)
ax3.plot(t,fxt,'r-')
ax3.set_title('Lest')
ax4 = plt.subplot(3,2,4)
ax4.plot(t,fzt,'b-')
ax4.set_title('Lest')

# Simulate random point sets and perform an envelope test
GET = importr("GET")
model = importr("spatstat.model")
def envelopeTest(x, nsim=1999, bounds = FloatVector([0,1]), fun=ro.Kest):
    X = ro.ppp(FloatVector(x[:,0]), FloatVector(x[:,1]), window = ro.owin(bounds,bounds))
    robjects.globalenv["X"] = X
    sim = rpy2.robjects.language.LangVector.from_string('expression(runifpoint(ex=X))')
    env = ro.envelope(X, fun=fun, nsim=nsim, savefuns=True, simulate=sim, verbose=False)
    res = ro.global_envelope_test(env)
    return (ro.attr(res,'p')[0],res.rx2('r'),res.rx2('lo'),res.rx2('hi'),res.rx2('central'),res.rx2('obs'))

(p,r,lo,hi,central,obs) = envelopeTest(x,fun=ro.Lest)
print(f'How likely x is to be random: {p:g}')
ax5 = plt.subplot(3,2,5)
ax5.plot(r,lo, 'k--')
ax5.plot(r,hi, 'k--')
ax5.plot(r,central, 'k-')
ax5.plot(r,obs, 'k-')
ax5.autoscale(enable=True, tight=True)
ax5.set_title(f'p={p:g}')

(p,r,lo,hi,central,obs) = envelopeTest(z,fun=ro.Lest)
print(f'How likely z is to be random: {p:g}')
ax6 = plt.subplot(3,2,6)
ax6.plot(r,lo, 'k--')
ax6.plot(r,hi, 'k--')
ax6.plot(r,central, 'k-')
ax6.plot(r,obs, 'k-')
ax6.autoscale(enable=True, tight=True)
ax6.set_title(f'p={p:g}')

plt.tight_layout()
plt.show()
