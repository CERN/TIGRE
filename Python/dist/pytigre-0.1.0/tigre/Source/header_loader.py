import os 

dirname = os.path.dirname(__file__)
lst = []
for filename in os.listdir(dirname):
    if filename.endswith('hpp'):
        lst.append(os.path.join('tigre','Source',filename))
print(lst)