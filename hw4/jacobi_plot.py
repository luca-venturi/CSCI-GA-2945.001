from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# plot strong scaling 

p = [64,256,1024,4096]

real_timings = [11.991474, 2.270839, 1.084813, 0.409620]

perfect_timings = [real_timings[0] / (4**j) for j in range(4)]

plt.plot(np.log2(p)*0.5,np.log2(real_timings)*0.5,np.log2(p)*0.5,np.log2(perfect_timings)*0.5)
plt.xlabel('$log_4(p)$')
plt.ylabel('$log_4(t)$')
plt.show()

# plot weak scaling 

p = [1,4,16,64,256,1024,4096]

real_timings = [19.936430, 20.289741, 20.513429, 21.221098, 21.253908, 21.748358, 22.681867]

perfect_timings = [real_timings[0] for j in range(7)]

plt.plot(np.log2(p)*0.5,real_timings,np.log2(p)*0.5,perfect_timings)
plt.ylim(0,30)
plt.xlabel('$log_4(p)$')
plt.ylabel('$t$')
plt.show()
