import pandas as pd
import matplotlib.pyplot as plt
data=pd.read_excel('data.xlsx',index=False,dtype={'Ra': float})
print(data)
plt.xscale('log')
plt.yscale('log')
plt.scatter(data['Ra'].values, data['Nu_top'].values, label='Nu_top')
plt.scatter(data['Ra'].values, data['Nu_bottom'].values, label='Nu_bottom')
plt.grid
plt.title("Nusselt en fonction de Ra")
plt.xlabel('Ra')
plt.ylabel('Nu')
plt.legend()
plt.show()
