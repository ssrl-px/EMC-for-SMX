import matplotlib.pyplot as plt
import numpy as np
import sys

ave_bg = np.fromfile(sys.argv[1], dtype=np.double)
plt.plot(range(len(ave_bg)), ave_bg)
plt.show()
