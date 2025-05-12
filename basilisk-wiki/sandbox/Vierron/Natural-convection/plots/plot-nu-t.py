import matplotlib.pyplot as plt
import pandas as pd

"""data = pd.read_csv('/home/administrator/cav2dSP/Ra=4e5/cav2d.asc',skiprows = [0],sep=' ', engine='python',index_col=False)"""

data = pd.read_csv('/home/administrator/cav2ds/out',skiprows = [0],sep=' ', engine='python',index_col=False)
data['[4]nutop'] = data['[4]nutop'].convert_objects(convert_numeric=True)
data['[5]nubot'] = data['[5]nubot'].convert_objects(convert_numeric=True)
data['[11]t'] = data['[11]t'].convert_objects(convert_numeric=True)
plt.plot(data['[11]t'].values,data['[4]nutop'].values,label='Nu_top')
plt.plot(data['[11]t'].values,data['[5]nubot'].values,label='Nu_bottom')
plt.title('Nu=f(t)')
plt.xlabel('t')
plt.ylabel('Nu')
plt.grid
plt.legend()
plt.show()