"""Independent analytic integral checks for the R PSF EE formula (not R execution)."""
import math,unittest
import numpy as np

def ee(fwhm,beta,q,diam,nang=720):
 alpha=(fwhm/2)/math.sqrt(2**(1/beta)-1)
 t=2*math.pi*((np.arange(nang)+.5)/nang)
 v=np.cos(t)**2+(np.sin(t)/q)**2
 return 1-np.mean((1/v)*(1+v*(diam/2/alpha)**2)**(1-beta))/q
class TestAnalyticEE(unittest.TestCase):
 def test_circular_closed_form(self):
  for d in (21,31,75):
   alpha=(3/2)/math.sqrt(2**(1/2.1)-1)
   expected=1-(1+(d/2/alpha)**2)**(1-2.1)
   self.assertAlmostEqual(ee(3,2.1,1,d),expected,places=12)
 def test_elliptical_monotonic(self):
  for q in (.35,.7,1.):
   self.assertTrue(0<ee(3,2.1,q,21)<ee(3,2.1,q,51)<1)
if __name__=='__main__':unittest.main()
