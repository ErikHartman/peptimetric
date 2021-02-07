import tkinter as tk
import matplotlib.pyplot as plt
import numpy as np

objects = ('Python', 'C++', 'Java', 'Perl', 'Scala', 'Lisp')
y_pos = np.arange(len(objects))
perfromance = [8,10,12,3,7,1]

plt.bar(y_pos,perfromance)
plt.xticks(y_pos, objects)

plt.show()

