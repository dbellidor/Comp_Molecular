from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import random
num_clusters=0
tam_imgX=0
tam_imgY=0
with open('out_30.txt', 'r') as reader:
	num_clusters = reader.readline()
	tam_imgX = reader.readline()
	tam_imgY = reader.readline()
	#print(tam_img)
num_clusters=int(num_clusters)
tam_imgX=int(tam_imgX)
tam_imgY=int(tam_imgY)
print(num_clusters)
print(tam_imgX)

file = open("out_30.txt","r")
file1=open("centroid_30.txt","r")
x=[]
y=[]
lbl_cluster=[]
r=[]
g=[]
b=[]
rc=[]
gc=[]
bc=[]
xc=[]
yc=[]
#color cluster


for j,line in enumerate(file1):
	fields=line.split(",")
	rc.append(int(fields[0]))
	gc.append(int(fields[1]))
	bc.append(int(fields[2]))
	xc.append(int(fields[3]))
	yc.append(int(fields[4]))
file1.close()

color_temp = np.zeros(shape=(num_clusters,3))
"""
print(rc[0],gc[0],bc[0])
print(rc[1],gc[1],bc[1])
print(rc[2],gc[2],bc[2])
print(rc[3],gc[3],bc[3])
"""

for i in range (len(rc)):
	for  j in range (i+1):
		color_temp[i][0]=rc[j]
		color_temp[i][1]=gc[j]
		color_temp[i][2]=bc[j]


color_temp = color_temp.tolist()
print(color_temp[0])
############################################


#color pixel
for i,line in enumerate(file):
	if (i>2):
		 fields = line.split(",")
		 r.append(int(fields[0]))
		 g.append(int(fields[1]))
		 b.append(int(fields[2]))
		 x.append(int(fields[3]))
		 y.append(int(fields[4]))
		 lbl_cluster.append(int(fields[5]))

		 
file.close()


w, h = tam_imgX, tam_imgY
data = np.zeros((h, w, 3), dtype=np.uint8)
for j in range (tam_imgY):
	for i in range (tam_imgX):
		get_color=lbl_cluster[(i+(j*tam_imgX))]
		data[j][i]=color_temp[get_color-1]
img = Image.fromarray(data, 'RGB')
img.save('img_'+str(num_clusters)+'.png')
img.show()



color_cluster=[]
for i in range (num_clusters):
	color_cluster.append(np.random.choice(range(256), size=3))


#print(color_cluster)



#Voronio

w, h = tam_imgX, tam_imgY
data = np.zeros((h, w, 3), dtype=np.uint8)
for j in range (tam_imgY):
	for i in range (tam_imgX):
		get_color=lbl_cluster[(i+(j*tam_imgX))]
		data[j][i]=color_cluster[get_color-1]
img = Image.fromarray(data, 'RGB')
img.save('voronoi_'+str(num_clusters)+'.png')
img.show()






