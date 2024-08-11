import numpy as np
x = np.loadtxt("x.csv",delimiter=",", dtype=float)
z = np.loadtxt("z.csv",delimiter=",", dtype=float)

import matplotlib.pyplot as plt
plt.plot(x[:,0],x[:,1],'k.')
plt.title('x')
plt.show()

plt.plot(z[:,0],z[:,1],'k.')
plt.title('z')
plt.show()

import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import matplotlib.pyplot as plt
import numpy as np
from rpy2.robjects import FloatVector
base = importr('base')
spatstat = importr('spatstat')
ro = robjects.r

def sampleK(x, n=100, bounds = FloatVector([0,1])):
    t = [i/float(n) for i in range(n+1)]
    ppx = ro.ppp(FloatVector(x[:,0]), FloatVector(x[:,1]), window = ro.owin(bounds,bounds))
    fx = base.as_function(ro.Lest(ppx))
    fxt = [fx(i)[0] for i in t] 
    return (t,fxt)

(t,fxt) = sampleK(x)
(_,fzt) = sampleK(z)
plt.plot(t,fxt,'r-')
plt.plot(t,fzt,'b-')
plt.show()

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
print(f'How likely x is to be random: {p}')
plt.plot(r,lo, 'k--')
plt.plot(r,hi, 'k--')
plt.plot(r,central, 'k-')
plt.plot(r,obs, 'k-')
plt.axis('tight')
plt.title(f'x, p={p}')
plt.show()

(p,r,lo,hi,central,obs) = envelopeTest(z,fun=ro.Lest)
print(f'How likely z is to be random: {p}')
plt.plot(r,lo, 'k--')
plt.plot(r,hi, 'k--')
plt.plot(r,central, 'k-')
plt.plot(r,obs, 'k-')
plt.axis('tight')
plt.title(f'x, p={p}')
plt.show()
