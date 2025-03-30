#!/usr/bin/env python
# coding: utf-8

# In[50]:


import matplotlib.pyplot as plt
import sys
from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv

def G(y, t, m1, l1, m2, l2): 

    # geschwindigkeiten
    a1d, a2d = y[0], y[1]

    # positionen
    a1, a2 = y[2], y[3]

    if not np.isfinite(a1) or not np.isfinite(a2):
        #print(a1, a2)
        raise ValueError("Ungültige Eingabewerte: a1 oder a2 sind nicht endlich")

    m11, m12 = (m1+m2)*l1, m2*l2*cos(a1-a2)
    m21, m22 = l1*cos(a1-a2), l2
    m = np.array([[m11, m12],[m21, m22]])

    f1 = -m2*l2*a2d*a2d*sin(a1-a2) - (m1+m2)*g*sin(a1)
    f2 = l1*a1d*a1d*sin(a1-a2) - g*sin(a2)
    f = np.array([f1, f2])

    accel = inv(m).dot(f)

    return np.array([accel[0], accel[1], a1d, a2d])


def RK4_step(y, t, dt, m1, l1, m2, l2):

	k1 = G(y,t, m1, l1, m2, l2) # len(m) x 4
	k2 = G(y+0.5*k1*dt, t+0.5*dt, m1, l1, m2, l2)
	k3 = G(y+0.5*k2*dt, t+0.5*dt, m1, l1, m2, l2)
	k4 = G(y+k3*dt, t+dt, m1, l1, m2, l2)

	return dt * (k1 + 2*k2 + 2*k3 + k4) /6



'''
def RK4_step(y, t, dt, m1, l1, m2, l2):

	k = G(y,t, m1, l1, m2, l2)

	return dt * k
'''



# parameters

#massen
#m1_values = np.arange(1,10,1)
#m2 = 1.0

#längen
#l1_values = np.arange(1,10,1)
#l1 = 1.0
#l2 = 1.0

#anfangswinkel
#a1, a2 = pi/4, -1.0


a1, a2 = 0.0, 0.0
g = 9.81


#prev_point = None
t = 0.0
delta_t = 0.02


# In[51]:




def time_evolution(t_stop,m1, l1, m2, l2, v1, v2):

    y = np.array([v1, v2, 0.0, 0.0])

    t=0

    t_simulation = np.array([])
    y_simulation = np.array([])

    while t<t_stop:

        try: 
            t_simulation = np.append(t_simulation,t)
            y_simulation = np.append(y_simulation,y)

            if not np.isfinite(y[2]) or not np.isfinite(y[3]):
                print(t)
                raise ValueError("a1 oder a2 sind unendlich geworden")

            
            y = y + RK4_step(y,t,delta_t,m1,l1,m2,l2)

        except ValueError as e:
            #print(f"Überspringe Iteration: {t} aufgrund von Fehler: {e}")
            t += delta_t
            continue

        t += delta_t

        

    y_simulation = y_simulation.reshape(len(t_simulation),4)

    return y_simulation, t_simulation

'''
y_simulation = time_evolution(5,m1[0],l1[0],m2[0],l2[0],1,1)[0]

t_simulation = np.arange(0,5+delta_t,delta_t)



# In[28]:


#winkel masse 1
winkel1 = y_simulation.reshape(len(t_simulation),4)[:,2] 

#winkel masse 2
winkel2 = y_simulation.reshape(len(t_simulation),4)[:,3]


# In[ ]:





# In[29]:
'''


def update(a1, a2, l1, l2):

    try:

        if not np.isfinite(a1) or not np.isfinite(a2):
            print(a1, a2)
            raise ValueError("Ungültige Eingabewerte: a1 oder a2 sind nicht endlich")

        x1 = l1* sin(a1)
        y1 = l1 * cos(a1)
        x2 = x1 + l2 * sin(a2)
        y2 = y1 + l2 * cos(a2)
    
    except ValueError as e:
        print(f"Überspringe Iteration: {a1} {a2} aufgrund von Fehler: {e}")

    return (x1, y1), (x2, y2)

'''

positions = [update(winkel1[i],winkel2[i],l1[0],l2[0]) for i in range(len(t_simulation))]

positions[0][1]


# In[30]:


x_positions_masse1 = [positions[i][0][0] for i in range(len(t_simulation))]
y_positions_masse1 = [positions[i][0][1] for i in range(len(t_simulation))]

x_positions_masse2 = [positions[i][1][0] for i in range(len(t_simulation))]
y_positions_masse2 = [positions[i][1][1] for i in range(len(t_simulation))]

plt.scatter(x_positions_masse1,y_positions_masse1,s=0.3)
plt.scatter(x_positions_masse2,y_positions_masse2, s=0.3)

plt.xlabel('x')
plt.ylabel('y')

plt.gca().invert_yaxis()  # Y-Achse umdrehen


# In[31]:
'''

'''
#m1_values = np.arange(1,10,0.1)
print('N_m', len(m1_values))
#l1_values = np.arange(1,10,1)
print('N_l', len(l1_values))

array = [[time_evolution(0.1,m1,l1) for m1 in m1_values] for l1 in l1_values]


print(len(array)) # len(array) = len(m1)
print(len(array[0]))
print(len(array[0][0]))
print(len(array[0][0][0]))

#array hat die dimension len(m) x len(l) x len(time) x 4
'''


# In[32]:


import pandas as pd


# In[103]:


def trajectory_images(m1_v,l1_v, m2_v, l2_v, e_v, v1_v, v2_v, dpi, labels_csv,t_stop,bilderzeigen=False):
        
        labels = []

        ratio_l = l1_v/l2_v

            
        ratio_l = (ratio_l - min(ratio_l))/(max(ratio_l)-min(ratio_l))

        ratio_m = m1_v/m2_v
        ratio_m = (ratio_m - min(ratio_m))/(max(ratio_m)-min(ratio_m))
    
        for i in range(len(m1_v)):
            
            m1 = m1_v[i]
            l1 = l1_v[i]
            #print(l1)
            m2 = m2_v[i]
            l2 = l2_v[i]
            #print(l2)
            e = e_v[i]
            v1 = v1_v[i]

            v2 = v2_v[i]

            r_m = ratio_m[i]
            r_l = ratio_l[i]


            try:
                y_simulation = time_evolution(t_stop,m1,l1,m2,l2, v1, v2)[0]
                t_simulation = time_evolution(t_stop,m1,l1,m2,l2, v1, v2)[1]
            except ValueError as e:
                print(m1,m2,l1,l2)
                print(f"Überspringe Iteration: {i} aufgrund von Fehler: {e}")
                continue

            #winkel masse 1
            winkel1 = y_simulation.reshape(len(t_simulation),4)[:,2] 

            #winkel masse 2
            winkel2 = y_simulation.reshape(len(t_simulation),4)[:,3]

            #print(l2)

            positions = [update(winkel1[j],winkel2[j],l1,l2) for j in range(len(t_simulation))]

            x_positions_masse1 = [positions[j][0][0] for j in range(len(t_simulation))]
            y_positions_masse1 = [positions[j][0][1] for j in range(len(t_simulation))]

            x_positions_masse2 = [positions[j][1][0] for j in range(len(t_simulation))]
            y_positions_masse2 = [positions[j][1][1] for j in range(len(t_simulation))]


            l1_n = (l1 - min(l1_v))/(max(l1_v)-min(l1_v))
            l2_n = (l2 - min(l2_v))/(max(l2_v)-min(l2_v))
            m1_n = (m1 - min(m1_v))/(max(m1_v)-min(m1_v))
            e_n = (e - min(e_v))/(max(e_v)-min(e_v))
            #v1_n = (v1 - min(v1))/(max(v1)-min(v1))
            m2_n = (m2 - min(m2_v))/(max(m2_v)-min(m2_v))
            #v2_n = (v2 - min(v2))/(max(v2)-min(v2))+1

            

            

            

            #print(m1_n/m2_n)
            

            filename = f'/Users/maraludwig/Documents/Master/Deep Learning/project/images/t_m{np.round(r_m,2)}_l{np.round(r_l,2)}_e{np.round(e_n,2)}.png'

            labels.append([f't_m{np.round(r_m,2)}_l{np.round(r_l,2)}_e{np.round(e_n,2)}.png', np.round(r_m,2), np.round(r_l,2), np.round(e_n,2)])
            

            
            plt.figure(figsize=(4, 4),dpi=dpi)

            plt.plot(x_positions_masse1,y_positions_masse1)
            plt.plot(x_positions_masse2,y_positions_masse2)

            plt.xlabel('x')
            plt.ylabel('y')

            L = (l1+l2)*1.05

            plt.ylim(-L,L)
            plt.xlim(-L,L)

            plt.gca().invert_yaxis()  # Y-Achse umdrehen



            plt.savefig(filename)

            plt.close()

        labels_df = pd.DataFrame(labels, columns=['image_name','mass_ratio', 'length_ratio', 'energy'])
        labels_df.to_csv(labels_csv, index=False)



# In[113]:


def sample_generator(num_samples, l_sum, v_max):


    # gauss distributed values of log(m2/m1)

    M = 1

    log_mass_ratio = np.random.normal(loc=0, scale=1, size=num_samples)

    m_ratio = np.exp(log_mass_ratio)
    m1 = M / (1 + m_ratio)
    m2 = M - m1


    # gauss distributed values of log(l2/l1)
    log_length_ratio = np.random.normal(loc=0, scale=1, size=num_samples)

    # aus verhälntnis und gesamtlänge l_sum die länge l1 und l2 berechnen:
    l_ratio = np.exp(log_length_ratio)
    l1 = l_sum / (1 + l_ratio)
    l2 = l_sum - l1

    v1 = np.random.uniform(low=-v_max, high=v_max, size=num_samples)
    np.random.shuffle(v1)
    v2 = np.random.uniform(low=-v_max, high=v_max, size=num_samples)
    np.random.shuffle(v2)

    # gesamtenergie (ohne potentielle energie)
    e_ges =1/2 * m1 * v1**2 + 1/2 * m2 * v2**2 + m1*g*l2
    return m1, m2, l1, l2, v1, v2, e_ges


'''
#testsamples generieren
m1, m2, l1, l2, v1, v2, e_ges = sample_generator(5, 1, 1)
'''



# In[35]:


import numpy as np



# In[114]:


#m1_values = np.array([1,2])
#l1_values = np.array([1,2])


labels_csv = '/Users/maraludwig/Documents/Master/Deep Learning/project/labels.csv'

'''
trajectory_images(m1,l1,m2,l2,e_ges,v1,v2,labels_csv,t_stop=5)'
'''


# In[115]:


import os
import pandas as pd
from PIL import Image
import torch
from torch.utils.data import Dataset, DataLoader, random_split
from torchvision import transforms

# Beispiel Dataset-Klasse für deine eigenen Bilder und Labels
class Dataset(Dataset):
    def __init__(self, image_dir, labels_csv, transform=None):
        """
        :param image_dir: Pfad zu den Bildern
        :param labels_csv: CSV-Datei mit den Labels (Masse, Längenverhältnis)
        :param transform: Optional: Transformationen für die Bilder
        """
        self.image_dir = image_dir
        self.labels_df = pd.read_csv(labels_csv)  # CSV mit den Labels
        self.transform = transform

    def __len__(self):
        return len(self.labels_df)

    def __getitem__(self, idx):
        try:
            # Holen des Bildpfads und der Labels
            img_name = os.path.join(self.image_dir, self.labels_df.iloc[idx, 0])  # Angenommen, die Bilddatei steht in der ersten Spalte
            image = Image.open(img_name).convert('RGB')

            # Holen der Labels (z.B. Masse, Längenverhältnis)
            mass = self.labels_df.iloc[idx, 1]  # Masse ist in der zweiten Spalte
            length_ratio = self.labels_df.iloc[idx, 2]  # Längenverhältnis in der dritten Spalte
            energy = self.labels_df.iloc[idx, 3]  # Energie in der vierten Spalte

            # Labels in ein Tensor umwandeln
            labels = torch.tensor([mass, length_ratio, energy], dtype=torch.float32)
            #labels = torch.tensor([mass, length_ratio], dtype=torch.float32)

            # Wenn Transformationen definiert sind, anwenden
            if self.transform:
                image = self.transform(image)

            return image, labels
        except Exception as e:
            print(f"Error loading data at index {idx}: {e}")
            return None, None

    def get_labels(self):
        return self.labels_df

# Beispielhafte Transformationen für Bilder


transform = transforms.Compose([
    transforms.Resize((400, 400)),
    transforms.ToTensor(),
])



# Dataset und DataLoader initialisieren
image_dir = '/Users/maraludwig/Documents/Master/Deep Learning/project/images'
labels_csv = '/Users/maraludwig/Documents/Master/Deep Learning/project/labels.csv'
dataset = Dataset(image_dir=image_dir, labels_csv=labels_csv,transform=transform)


# Dataset in Trainings- und Validierungsdaten aufteilen
train_size = int(0.8 * len(dataset))
val_size = len(dataset) - train_size
train_dataset, val_dataset = random_split(dataset, [train_size, val_size])

# DataLoader erstellen
train_loader = DataLoader(train_dataset, batch_size=32, shuffle=True, num_workers=0)
val_loader = DataLoader(val_dataset, batch_size=32, shuffle=False, num_workers=0)

dataset.get_labels()


def show_first_example(dataset,index):
#dataset.__getitem__(1)

    print('Original data:')
    plt.imshow(dataset[index][0].permute(1,2,0).cpu().numpy())
    plt.axis('off')
    print(dataset[0][1])



import torchvision
import torch.nn as nn

# ResNet18 für Regression anpassen
model = torchvision.models.resnet18(weights=torchvision.models.ResNet18_Weights.IMAGENET1K_V1)
model.fc = nn.Linear(model.fc.in_features, 3)  # drei Ausgaben (Massenverhältnis, Längenverhältnis, energy)
#model = model.cuda()

# Verlustfunktion für Regression
loss_fn = torch.nn.MSELoss()


# In[109]:


import time

import pickle

loss_array = []

def train(args, model, device, train_loader, optimizer, epoch, loss_fn):

    model.train()
    start_time = time.time()

    for batch_idx, (data, target) in enumerate(train_loader):
        

        data, target = data.to(device), target.to(device)

        optimizer.zero_grad()
        output = model(data)

        loss = loss_fn(output, target)  # MSE Loss für Regression
        loss_array.append(loss.item())
        loss.backward()
        optimizer.step()

        if batch_idx % args.log_interval == 0:
            print(f'Train Epoch: {epoch} [{batch_idx * len(data)}/{len(train_loader.dataset)} ({100. * batch_idx / len(train_loader):.0f}%)]\tLoss: {loss.item():.6f}')
            if args.dry_run:
                break

    #print('loss_array',loss_array)

    with open("/Users/maraludwig/Documents/Master/Deep Learning/project/train_loss.pkl", "wb") as f:
        pickle.dump(loss_array, f)
    
    print("--- Epoch time: %s seconds ---" % (time.time() - start_time))


t_loss_array = []

def test(model, device, test_loader, loss_fn):
    model.eval()
    test_loss = 0
    with torch.no_grad():
        for data, target in test_loader:
            data, target = data.to(device), target.to(device)
            output = model(data)

            test_loss += loss_fn(output, target).item()

            t_loss_array.append(test_loss)

    test_loss /= len(test_loader)

    #t_loss_array.append(test_loss)

    with open("/Users/maraludwig/Documents/Master/Deep Learning/project/test_loss.pkl", "wb") as f:
        pickle.dump(t_loss_array, f)

    print(f'\nTest set: Average loss: {test_loss:.4f}\n')

optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)


# In[42]:


import types


# In[118]:

'''
num_epochs = 20
args = types.SimpleNamespace(dry_run=False, log_interval=16)

device = "cpu"
for epoch in range(num_epochs):
    train(args, model, device, train_loader, optimizer, epoch, loss_fn)
    test(model, device, val_loader, loss_fn)  # Or with a test loader, if available'
'''


# In[44]:


# Funktion zum Anzeigen eines Beispiels
def show_example(model, dataset, index, device):
    model.eval()
    image, labels = dataset[index]
    if image is None or labels is None:
        print(f"Warning: Found None data at index {index}")
        return

    image = image.unsqueeze(0).to(device)  # Füge eine Batch-Dimension hinzu
    with torch.no_grad():
        output = model(image)

    # Konvertiere das Bild zurück zu einem numpy-Array für die Anzeige
    image = image.squeeze(0).cpu().numpy().transpose((1, 2, 0))

    # Rescale the image to [0, 1] range for display
    #image = (image - image.min()) / (image.max() - image.min())

    # Zeige das Bild und die Vorhersagen
    plt.imshow(image)
    plt.title(f"True: {labels.cpu().numpy()}, Predicted: {output.cpu().numpy()}")
    plt.axis('off')
    plt.show()

# Beispiel anzeigen
'''
show_example(model, val_dataset, index=1, device=device)'
'''

