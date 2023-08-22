import numpy as np
import numpy.fft as nf
def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth
def SmoothySpec(a,nums=None):
   b=a.copy()
   if nums is None: nums=2*len(b)//3
   for i in range(nums):
      b[i+1::][:-1] = 0.25*b[i::][:-2]+0.5*b[i+1::][:-1]+0.25*b[i::][2:]
   return b

def SlidingFSpec(a_in,dt,winkind='tukey',nums=10,step=100):
   """
      Code to compute power spectrum using simple
      minded Fourier transform on windowed data.
   """
   import TurbAn.Analysis.Simulations.AnalysisFunctions as af
   nt=len(a_in) # Length of input timeseries
   ll = nt - (nums-1)*step # Length of time-series to work with
   lenspec=ll//2 #
   if lenspec%2:
      lenspec=lenspec-1
      aa=a_in[2:].copy()
   else:
      aa=a_in.copy()
   spec=np.zeros(lenspec+1)
   for i in range(nums):
      a=aa[i*step:i*step+ll]*af.windowff(ll,kind=winkind)
      fa=nf.rfft(a)/len(a)
      print(i,len(spec),len(fa))
      spec+=abs(fa)**2 
   fq=nf.rfftfreq(len(a),d=dt)
   return fq,spec/nums

