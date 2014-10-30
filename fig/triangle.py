import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

plt.rc('text', usetex=True)
fig, ax = plt.subplots()

Path = mpath.Path
# path_data = [
#     (Path.MOVETO, (0, 0)),
#     (Path.LINETO, (1, 0)),
#     (Path.LINETO, (0, 1)),
#     ]
# codes, verts = zip(*path_data)
A = (0, 0)
B = (1, 0)
C = (0.9, 0.8)
verts = (A, B, C, A)
codes = (Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY)
path = mpath.Path(verts, codes)
patch = mpatches.PathPatch(path, facecolor='r', alpha=0.2)
ax.add_patch(patch)

ax.text(A[0]-.05, A[1]-.08, r'$(x_0, y_0)$', fontsize=15)
ax.text(B[0]+.00, B[1]-.08, r'$(x_1, y_1)$', fontsize=15)
ax.text(C[0]+.00, C[1]+.02, r'$(x_2, y_2)$', fontsize=15)

# plot control points and connecting lines
x, y = zip(*path.vertices)
line, = ax.plot(x, y, 'ko-')

#ax.grid()
#ax.axis('equal')
ax.axis('off')
ax.axis([-0.1, 1.1, -0.1, 1.1])

plt.show()

